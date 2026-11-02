#!/usr/bin/env python3
"""Protein-sharing network and gene-content clustering for all 241 genomes.

CD-HIT at 50% identity / 50% coverage defines protein families. Prophages are
clustered by Jaccard similarity of protein-family content (average linkage,
distance <= 0.5). Host relatedness is measured by Mash whole-genome distance,
not MLST alone.
"""
import json
import os
import re
from collections import Counter

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import fisher_exact, mannwhitneyu
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

BASE = "/home/bioinfo/cdprophage"
RESULT = os.path.join(BASE, "scaleup_241", "results")


def parse_clstr(path):
    fam = {}
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">Cluster"):
                current = int(line.split()[1])
            elif re.match(r"^\d+\t", line):
                if ", >" not in line:
                    continue
                header = line.split(", >")[1].split("...")[0]
                fam[header] = "GF%05d" % (current + 1)
    return fam


def mantel(m1, m2, permutations=999, seed=42):
    k = m1.shape[0]
    iu = np.triu_indices_from(m1, k=1)
    a, b = m1[iu], m2[iu]
    r_obs = float(np.corrcoef(a, b)[0, 1])
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(permutations):
        perm = rng.permutation(k)
        rp = float(np.corrcoef(a, m2[np.ix_(perm, perm)][iu])[0, 1])
        if rp >= r_obs:
            count += 1
    return r_obs, (count + 1) / (permutations + 1)


def jaccard_matrix(X):
    sizes = X.sum(axis=1)
    inter = X @ X.T
    union = sizes[:, None] + sizes[None, :] - inter
    J = np.divide(inter, union, out=np.zeros_like(inter, dtype=float), where=union > 0)
    np.fill_diagonal(J, 1.0)
    return J


def main():
    details = pd.read_csv(os.path.join(RESULT, "prophage_details.csv"),
                          dtype={"ST": str, "mlst_clade": str})
    cds = pd.read_csv(os.path.join(RESULT, "prophage_cds_annotation.csv"),
                      dtype={"ST": str, "mlst_clade": str})
    cargo = pd.read_csv(os.path.join(RESULT, "prophage_cargo_summary.csv"),
                        dtype={"ST": str, "mlst_clade": str})
    fam = parse_clstr(os.path.join(RESULT, "prophage_proteins_nr.fasta.clstr"))

    key_to_cds = dict(zip(cds["protein_key"], cds.index))
    rows = []
    for key, gf in fam.items():
        if key not in key_to_cds:
            continue
        c = cds.loc[key_to_cds[key]]
        rows.append({
            "protein_key": key,
            "prophage_id": c["prophage_id"],
            "assembly_accession": c["assembly_accession"],
            "ST": c["ST"],
            "mlst_clade": c["mlst_clade"],
            "gene_family": gf,
            "product": c["product"],
            "phage_core": int(c["phage_core"]),
            "is_defence": int(c["is_defence"]),
            "is_virulence": int(c["is_virulence"]),
            "is_amr": int(c["is_amr"]),
            "is_mge": int(c["is_mge"]),
        })
    prots = pd.DataFrame(rows)

    pids = sorted(details["prophage_id"])
    gfs = sorted(prots["gene_family"].unique())
    pid_index = {p: i for i, p in enumerate(pids)}
    gf_index = {g: i for i, g in enumerate(gfs)}
    X = np.zeros((len(pids), len(gfs)), dtype=np.int8)
    for _, r in prots.iterrows():
        X[pid_index[r["prophage_id"]], gf_index[r["gene_family"]]] = 1

    J = jaccard_matrix(X)
    D = 1.0 - J
    np.fill_diagonal(D, 0.0)
    Z = linkage(squareform(D, checks=False), method="average")
    gc_labels = fcluster(Z, t=0.5, criterion="distance")

    meta = details.set_index("prophage_id")
    gf_rows = []
    for gf, sub in prots.groupby("gene_family"):
        genes = [g for g in sub["product"].fillna("").tolist() if g]
        rep = Counter(genes).most_common(1)[0][0] if genes else ""
        cats = []
        if sub["is_defence"].any():
            cats.append("defence")
        if sub["is_virulence"].any():
            cats.append("virulence")
        if sub["is_amr"].any():
            cats.append("AMR")
        if sub["is_mge"].any():
            cats.append("MGE")
        if sub["phage_core"].any():
            cats.append("phage_core")
        gf_rows.append({
            "gene_family": gf,
            "n_proteins": len(sub),
            "n_prophages": sub["prophage_id"].nunique(),
            "n_genomes": sub["assembly_accession"].nunique(),
            "n_clades": sub["mlst_clade"].nunique(),
            "n_nucleotide_clusters": sub["prophage_id"].map(meta["cluster_id"]).nunique(),
            "representative_product": rep,
            "categories": ";".join(cats) if cats else "other",
            "has_defence": int(sub["is_defence"].any()),
            "has_virulence": int(sub["is_virulence"].any()),
            "has_amr": int(sub["is_amr"].any()),
            "has_mge": int(sub["is_mge"].any()),
            "has_phage_core": int(sub["phage_core"].any()),
        })
    gf_df = pd.DataFrame(gf_rows).sort_values(
        ["n_prophages", "n_clades", "n_genomes"], ascending=False)
    gf_df.to_csv(os.path.join(RESULT, "protein_gene_families.csv"), index=False)

    per = details[["prophage_id", "assembly_accession", "ST", "mlst_clade",
                   "cluster_id", "length_bp", "confidence"]].copy()
    per = per.set_index("prophage_id").loc[pids].reset_index()
    gc_map = {pid: "GC%03d" % lab for pid, lab in zip(pids, gc_labels)}
    per["gc_cluster"] = [gc_map[p] for p in per["prophage_id"]]
    per["n_gene_families"] = X.sum(axis=1)
    core_counts = []
    for p in pids:
        core_counts.append(prots[(prots["prophage_id"] == p) & (prots["phage_core"] == 1)]["gene_family"].nunique())
    per["n_phage_core_families"] = core_counts
    per = per.merge(cargo[["prophage_id", "has_defence_system", "has_virulence",
                            "has_amr", "has_mge"]], on="prophage_id", how="left")
    per.to_csv(os.path.join(RESULT, "prophage_gene_content.csv"), index=False)

    ari = adjusted_rand_score(per["cluster_id"], per["gc_cluster"])
    nmi = normalized_mutual_info_score(per["cluster_id"], per["gc_cluster"])
    gc_sizes = per["gc_cluster"].value_counts()
    gc_singletons = int((gc_sizes == 1).sum())

    clades = per["mlst_clade"].astype(str).values
    iu = np.triu_indices_from(J, k=1)
    jvals = J[iu]
    within = np.array([clades[i] == clades[j] for i, j in zip(*iu)])
    mw = mannwhitneyu(jvals[within], jvals[~within], alternative="greater")

    # Genome-level union avoids pseudoreplication.
    genome_ids = sorted(per["assembly_accession"].unique())
    gidx = {g: i for i, g in enumerate(genome_ids)}
    Xg = np.zeros((len(genome_ids), len(gfs)), dtype=np.int8)
    for _, r in prots.iterrows():
        Xg[gidx[r["assembly_accession"]], gf_index[r["gene_family"]]] = 1
    Jg = jaccard_matrix(Xg)

    mash = pd.read_csv(os.path.join(RESULT, "mash_dist_matrix.csv"), index_col=0)
    mash.index = mash.index.astype(str)
    mash.columns = mash.columns.astype(str)
    gids = [g for g in genome_ids if g in mash.index and g in mash.columns]
    gm = mash.loc[gids, gids].values.astype(float)
    gJ = Jg[[gidx[g] for g in gids]][:, [gidx[g] for g in gids]]
    r_gen, p_gen = mantel(1.0 - gJ, gm, permutations=999, seed=13)

    # Presence/absence of nucleotide clusters at genome level vs Mash.
    Xp = np.zeros((len(genome_ids), len(details["cluster_id"].unique())), dtype=np.int8)
    clim = {c: i for i, c in enumerate(sorted(details["cluster_id"].unique()))}
    for _, r in per.iterrows():
        Xp[gidx[r["assembly_accession"]], clim[r["cluster_id"]]] = 1
    Jp = jaccard_matrix(Xp)
    gpids = [g for g in genome_ids if g in mash.index]
    Jp_sub = Jp[[gidx[g] for g in gpids]][:, [gidx[g] for g in gpids]]
    gm2 = mash.loc[gpids, gpids].values.astype(float)
    r_pres, p_pres = mantel(1.0 - Jp_sub, gm2, permutations=999, seed=17)

    cat_counts = {}
    cat_totals = {}
    for cat, col in (("defence", "has_defence"), ("virulence", "has_virulence"),
                     ("AMR", "has_amr"), ("MGE", "has_mge"), ("phage_core", "has_phage_core")):
        sub = gf_df[gf_df[col] == 1]
        cat_totals[cat] = int(len(sub))
        cat_counts[cat] = int((sub["n_clades"] >= 2).sum())
    cross_def = cat_counts["defence"]
    non_def = cat_totals["defence"] - cross_def
    cross_other = sum(cat_counts[c] for c in ("virulence", "AMR", "MGE", "phage_core"))
    non_other = sum(cat_totals[c] for c in ("virulence", "AMR", "MGE", "phage_core")) - cross_other
    odds, p_cross = fisher_exact([[cross_def, non_def], [cross_other, non_other]])

    summary = {
        "proteins": int(len(prots)),
        "gene_families": int(len(gfs)),
        "prophages": int(len(per)),
        "gene_content_clusters": int(per["gc_cluster"].nunique()),
        "gene_content_singleton_clusters": gc_singletons,
        "gene_content_singleton_rate": round(gc_singletons / per["gc_cluster"].nunique(), 4),
        "ari_vs_nucleotide_clusters": round(float(ari), 3),
        "nmi_vs_nucleotide_clusters": round(float(nmi), 3),
        "jaccard_within_clade_median": round(float(np.median(jvals[within])), 3),
        "jaccard_between_clade_median": round(float(np.median(jvals[~within])), 3),
        "jaccard_within_clade_mean": round(float(np.mean(jvals[within])), 4),
        "jaccard_between_clade_mean": round(float(np.mean(jvals[~within])), 4),
        "fraction_pairs_jaccard_gt_0.2_within": round(float(np.mean(jvals[within] > 0.2)), 4),
        "fraction_pairs_jaccard_gt_0.2_between": round(float(np.mean(jvals[~within] > 0.2)), 4),
        "mannwhitney_within_gt_between_p": float(mw.pvalue),
        "mantel_genome_gene_content_vs_mash": {"r": round(r_gen, 3), "p": p_gen, "n": int(len(gids))},
        "mantel_genome_cluster_presence_vs_mash": {"r": round(r_pres, 3), "p": p_pres, "n": int(len(gpids))},
        "cross_clade_gene_families": cat_counts,
        "gene_families_by_category": cat_totals,
        "defence_cross_clade_enrichment": {
            "odds_ratio": round(float(odds), 3), "p_value": float(p_cross),
            "table": [[cross_def, non_def], [cross_other, non_other]],
        },
        "top_widespread_gene_families": gf_df.head(15)[
            ["gene_family", "n_prophages", "n_genomes", "n_clades", "categories", "representative_product"]
        ].to_dict(orient="records"),
    }
    with open(os.path.join(RESULT, "protein_network_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
