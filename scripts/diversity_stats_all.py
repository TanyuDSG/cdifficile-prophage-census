#!/usr/bin/env python3
"""Diversity statistics for C. difficile prophages across all 241 genomes."""
import json
import math
import os

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform

BASE = "/home/bioinfo/cdprophage"
RESULT = os.path.join(BASE, "scaleup_241", "results")
LAB_DERIVATIVE = "GCF_000953275.1"  # 630Derm, lab derivative of 630 (ST54)

details = pd.read_csv(os.path.join(RESULT, "prophage_details.csv"), dtype={"ST": str, "mlst_clade": str})
genomes = pd.read_csv(os.path.join(RESULT, "genome_prophage_summary_raw.csv"), dtype={"ST": str, "mlst_clade": str})
mash = pd.read_csv(os.path.join(RESULT, "mash_dist_matrix.csv"), index_col=0)
mash.index = mash.index.astype(str)
mash.columns = mash.columns.astype(str)


def analyse(details, genomes, tag=""):
    all_accs = list(genomes["assembly_accession"])
    clusters = sorted(details["cluster_id"].dropna().unique())
    presence = pd.DataFrame(0, index=all_accs, columns=clusters)
    for _, r in details.iterrows():
        presence.loc[r["assembly_accession"], r["cluster_id"]] = 1

    n_clusters = len(clusters)
    cluster_genomes = presence.sum(axis=0)
    cluster_sts = details.groupby("cluster_id")["ST"].nunique()
    cluster_clades = details.groupby("cluster_id")["mlst_clade"].nunique()
    singleton_clusters = int((cluster_genomes == 1).sum())
    sharing = {
        "one_genome": int((cluster_genomes == 1).sum()),
        "two_to_three": int(((cluster_genomes >= 2) & (cluster_genomes <= 3)).sum()),
        "four_to_ten": int(((cluster_genomes >= 4) & (cluster_genomes <= 10)).sum()),
        "eleven_to_fifty": int(((cluster_genomes >= 11) & (cluster_genomes <= 50)).sum()),
        "more_than_fifty": int((cluster_genomes > 50).sum()),
    }

    p = cluster_genomes / cluster_genomes.sum()
    shannon = float(-(p * np.log(p)).sum())
    simpson = float(1 - (p ** 2).sum())
    evenness = float(shannon / math.log(n_clusters)) if n_clusters > 1 else 0.0

    genome_metrics = []
    for acc in all_accs:
        sub = details[details["assembly_accession"] == acc]
        counts = sub["cluster_id"].value_counts()
        sh = float(-((counts / counts.sum()) * np.log(counts / counts.sum())).sum()) if len(counts) > 1 else 0.0
        genome_metrics.append({
            "assembly_accession": acc,
            "n_clusters": sub["cluster_id"].nunique(),
            "n_high": int((sub["confidence"] == "high").sum()),
            "n_complete_like": int(sub["complete_like"].sum()),
            "n_typeable": int(sub["phage_typeable"].sum()),
            "within_genome_shannon": round(sh, 3),
        })
    gm = pd.DataFrame(genome_metrics)
    g = genomes.merge(gm, on="assembly_accession", how="left")
    g["carriage"] = (g["n_prophages"] > 0).astype(int)
    g["high_carriage"] = (g["n_high"] > 0).astype(int)
    g["complete_carriage"] = (g["n_complete_like"] > 0).astype(int)

    def desc(series):
        s = pd.to_numeric(series, errors="coerce").dropna()
        return {
            "n": int(len(s)), "mean": round(float(s.mean()), 3),
            "sd": round(float(s.std(ddof=1)), 3) if len(s) > 1 else 0,
            "median": round(float(s.median()), 3), "min": float(s.min()), "max": float(s.max()),
            "q1": round(float(s.quantile(0.25)), 3), "q3": round(float(s.quantile(0.75)), 3),
        }

    summary = {
        "dataset": tag,
        "genomes": int(len(g)),
        "STs": int(g["ST"].nunique()),
        "clades": int(g["mlst_clade"].nunique()),
        "prophages": int(len(details)),
        "carriage_genomes": int(g["carriage"].sum()),
        "carriage_rate": round(float(g["carriage"].mean()), 4),
        "high_confidence_carriage": int(g["high_carriage"].sum()),
        "high_confidence_carriage_rate": round(float(g["high_carriage"].mean()), 4),
        "complete_like_carriage": int(g["complete_carriage"].sum()),
        "prophages_per_genome": desc(g["n_prophages"]),
        "prophage_percent": desc(g["prophage_percent"]),
        "prophage_length_bp": desc(details["length_bp"]),
        "prophage_gc_percent": desc(details["gc_percent"]),
        "high_confidence_n": int((details["confidence"] == "high").sum()),
        "medium_confidence_n": int((details["confidence"] == "medium").sum()),
        "low_confidence_n": int((details["confidence"] == "low").sum()),
        "complete_like_n": int(details["complete_like"].sum()),
        "typeable_n": int(details["phage_typeable"].sum()),
        "typeable_rate": round(float(details["phage_typeable"].mean()), 4),
        "clusters": n_clusters,
        "singleton_clusters": singleton_clusters,
        "singleton_cluster_rate": round(singleton_clusters / n_clusters, 4) if n_clusters else 0,
        "sharing": sharing,
        "shannon_cluster": round(shannon, 3),
        "simpson_cluster": round(simpson, 3),
        "evenness_cluster": round(evenness, 3),
        "clusters_per_genome": desc(g["n_clusters"]),
        "clusters_in_multiple_clades": int((cluster_clades >= 2).sum()),
        "clusters_in_multiple_STs": int((cluster_sts >= 2).sum()),
        "core_clusters_95pct": int((cluster_genomes >= 0.95 * len(g)).sum()),
    }

    # clade summary
    clade_rows = []
    for clade, sub in g.groupby("mlst_clade"):
        accs = set(sub["assembly_accession"])
        dsub = details[details["assembly_accession"].isin(accs)]
        cl_clusters = set(dsub["cluster_id"])
        other_clusters = set(details[~details["assembly_accession"].isin(accs)]["cluster_id"])
        top = dsub["cluster_id"].value_counts().head(3)
        clade_rows.append({
            "mlst_clade": clade, "n_genomes": len(sub), "n_STs": sub["ST"].nunique(),
            "n_prophages": len(dsub), "mean_prophages": round(len(dsub) / len(sub), 3),
            "carriage_rate": round(float((sub["n_prophages"] > 0).mean()), 3),
            "high_carriage_rate": round(float((sub["n_high"] > 0).mean()), 3),
            "mean_prophage_percent": round(float(sub["prophage_percent"].mean()), 3),
            "mean_clusters": round(float(sub["n_clusters"].mean()), 3),
            "clusters_present": len(cl_clusters),
            "unique_clusters": len(cl_clusters - other_clusters),
            "shared_clusters": len(cl_clusters & other_clusters),
            "top_clusters": ";".join(f"{c}({n})" for c, n in top.items()),
        })
    clade_df = pd.DataFrame(clade_rows).sort_values("mlst_clade")
    clade_df.to_csv(os.path.join(RESULT, f"clade_summary{tag}.csv"), index=False)

    # ST summary (>=5 genomes)
    st_rows = []
    for st, sub in g.groupby("ST"):
        if len(sub) < 5:
            continue
        accs = set(sub["assembly_accession"])
        dsub = details[details["assembly_accession"].isin(accs)]
        st_rows.append({
            "ST": st, "mlst_clade": sub["mlst_clade"].iloc[0], "n_genomes": len(sub),
            "n_prophages": len(dsub), "mean_prophages": round(len(dsub) / len(sub), 3),
            "n_clusters": dsub["cluster_id"].nunique(),
            "mean_prophage_percent": round(float(sub["prophage_percent"].mean()), 3),
        })
    pd.DataFrame(st_rows).sort_values(["n_genomes", "ST"], ascending=[False, True]).to_csv(
        os.path.join(RESULT, f"st_summary{tag}.csv"), index=False)

    # statistical tests
    tests = []
    clade_groups = [s["n_prophages"].values for _, s in g.groupby("mlst_clade") if len(s) >= 10]
    if len(clade_groups) >= 2:
        h, pval = stats.kruskal(*clade_groups)
        tests.append({"test": "Kruskal-Wallis", "metric": "prophages_per_genome",
                      "groups": "MLST clade (n>=10)", "statistic": round(float(h), 3),
                      "p_value": float(pval), "note": ""})
    clade_groups_pct = [s["prophage_percent"].values for _, s in g.groupby("mlst_clade") if len(s) >= 10]
    if len(clade_groups_pct) >= 2:
        h, pval = stats.kruskal(*clade_groups_pct)
        tests.append({"test": "Kruskal-Wallis", "metric": "prophage_bp_percent",
                      "groups": "MLST clade (n>=10)", "statistic": round(float(h), 3),
                      "p_value": float(pval), "note": ""})

    tab = pd.crosstab(g["mlst_clade"] == "1", g["high_carriage"])
    if tab.shape == (2, 2):
        odds, pval = stats.fisher_exact(tab.values)
        tests.append({"test": "Fisher exact", "metric": "high_confidence_carriage",
                      "groups": "clade 1 vs other clades", "statistic": round(float(odds), 3),
                      "p_value": float(pval), "note": f"table={tab.values.tolist()}"})

    for x, y in [("genome_size", "n_prophages"), ("genome_size", "prophage_percent"),
                 ("genome_gc_percent", "n_prophages"), ("n_prophages", "prophage_percent")]:
        rho, pval = stats.spearmanr(g[x], g[y])
        tests.append({"test": "Spearman", "metric": f"{x} vs {y}", "groups": "all genomes",
                      "statistic": round(float(rho), 3), "p_value": float(pval), "note": ""})

    # Mantel: cluster presence vs Mash distance (genomes with Mash data)
    accs_m = [a for a in all_accs if a in mash.index]
    X = presence.loc[accs_m].values.astype(bool)
    n = len(accs_m)
    jac = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            union = np.logical_or(X[i], X[j]).sum()
            inter = np.logical_and(X[i], X[j]).sum()
            jac[i, j] = jac[j, i] = 0.0 if union == 0 else 1 - inter / union

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

    md = mash.loc[accs_m, accs_m].values.astype(float)
    r_all, p_all = mantel(jac, md)
    tests.append({"test": "Mantel", "metric": "prophage cluster Jaccard vs Mash genome distance",
                  "groups": f"all genomes (n={n})", "statistic": round(r_all, 3),
                  "p_value": p_all, "note": "999 permutations, one-sided"})

    clade1 = [a for a in accs_m if g.set_index("assembly_accession").loc[a, "mlst_clade"] == "1"]
    if len(clade1) >= 5:
        idx1 = [accs_m.index(a) for a in clade1]
        r_c1, p_c1 = mantel(jac[np.ix_(idx1, idx1)], md[np.ix_(idx1, idx1)], seed=43)
        tests.append({"test": "Mantel", "metric": "prophage cluster Jaccard vs Mash genome distance",
                      "groups": f"clade 1 only (n={len(clade1)})", "statistic": round(r_c1, 3),
                      "p_value": p_c1, "note": "999 permutations, one-sided"})
    else:
        r_c1, p_c1 = None, None

    ps = np.array([t["p_value"] for t in tests])
    order = np.argsort(ps)
    q = np.empty_like(ps)
    m = len(ps)
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        idx = order[rank]
        val = min(prev, ps[idx] * m / (rank + 1))
        q[idx] = val
        prev = val
    for t, qv in zip(tests, q):
        t["q_value_BH"] = round(float(qv), 4)
    pd.DataFrame(tests).to_csv(os.path.join(RESULT, f"statistical_tests{tag}.csv"), index=False)

    summary["mantel_all"] = {"r": round(r_all, 3), "p": p_all, "n": n}
    summary["mantel_clade1"] = {"r": round(r_c1, 3) if r_c1 is not None else None,
                                "p": p_c1, "n": len(clade1)}
    summary["top_clusters"] = [
        {"cluster_id": c, "n_genomes": int(n_g), "n_clades": int(cluster_clades[c])}
        for c, n_g in cluster_genomes.sort_values(ascending=False).head(15).items()
    ]

    # accumulation curve
    rng = np.random.default_rng(7)
    curves = []
    cluster_arr = np.array(clusters)
    for _ in range(500):
        order = rng.permutation(len(all_accs))
        seen = set()
        curve = [0]
        for idx in order:
            seen |= set(cluster_arr[X[idx]])
            curve.append(len(seen))
        curves.append(curve)
    curve_mean = np.mean(curves, axis=0)
    acc_df = pd.DataFrame({"n_genomes": np.arange(len(curve_mean)), "clusters_mean": curve_mean,
                           "clusters_p05": np.percentile(curves, 5, axis=0),
                           "clusters_p95": np.percentile(curves, 95, axis=0)})
    acc_df.to_csv(os.path.join(RESULT, f"cluster_accumulation{tag}.csv"), index=False)
    summary["rarefaction_57_genomes"] = round(float(curve_mean[57]), 2)
    summary["rarefaction_100_genomes"] = round(float(curve_mean[100]), 2)
    if len(curve_mean) > 57:
        summary["percent_clusters_at_57"] = round(float(curve_mean[57] / curve_mean[-1] * 100), 1)

    # reference phage types
    type_rows = []
    for phage, sub in details[details["phage_typeable"] == 1].groupby("best_phage"):
        type_rows.append({"best_phage": phage, "n_prophages": len(sub),
                          "n_genomes": sub["assembly_accession"].nunique(),
                          "n_clusters": sub["cluster_id"].nunique(),
                          "clades": ";".join(sorted(sub["mlst_clade"].unique())),
                          "median_pident": round(float(sub["best_pident"].median()), 2),
                          "median_qcov": round(float(sub["best_qcov"].median()), 3)})
    pd.DataFrame(type_rows).sort_values("n_prophages", ascending=False).to_csv(
        os.path.join(RESULT, f"reference_phage_types{tag}.csv"), index=False)

    summary["reference_phages_detected"] = int(len(type_rows))
    return summary, g, presence


main, main_genomes, presence = analyse(details, genomes, tag="")
main_genomes.to_csv(os.path.join(RESULT, "genome_summary_analysis.csv"), index=False)

sens_details = details[details["assembly_accession"] != LAB_DERIVATIVE]
sens_genomes = genomes[genomes["assembly_accession"] != LAB_DERIVATIVE]
sens, _, _ = analyse(sens_details, sens_genomes, tag="_no630Derm")

main["sensitivity_without_630Derm"] = {k: sens[k] for k in
                                       ("genomes", "prophages", "clusters", "singleton_clusters",
                                        "carriage_rate", "prophages_per_genome", "prophage_percent")}
with open(os.path.join(RESULT, "diversity_summary.json"), "w") as fh:
    json.dump(main, fh, indent=2)
print(json.dumps(main, indent=2))
