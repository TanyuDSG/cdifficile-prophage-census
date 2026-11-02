#!/usr/bin/env python3
"""Post-hoc analyses requested by the second-round review.

Outputs
  mantel_partial.json       Mantel and partial Mantel (control year, country)
  positive_control_630.json phiCD630-1/-2 recovery and cluster merging
  reference_similarity.json module-level vs genome-level similarity
  defence_overlap.json      keyword cargo vs DefenseFinder systems
"""
import json
import os

import numpy as np
import pandas as pd

BASE = "/home/bioinfo/cdprophage"
RES = os.path.join(BASE, "scaleup_241", "results")

details = pd.read_csv(os.path.join(RES, "prophage_details.csv"),
                      dtype={"ST": str, "mlst_clade": str})
meta = pd.read_csv(os.path.join(RES, "genome_metadata_full.csv"),
                   dtype={"ST": str, "mlst_clade": str})
mash = pd.read_csv(os.path.join(RES, "mash_dist_matrix.csv"), index_col=0)
mash.index = mash.index.astype(str)
mash.columns = mash.columns.astype(str)

accs = sorted(details["assembly_accession"].unique())
clusters = sorted(details["cluster_id"].dropna().unique())
presence = pd.DataFrame(0, index=accs, columns=clusters)
for _, r in details.iterrows():
    presence.loc[r["assembly_accession"], r["cluster_id"]] = 1
X = presence.values.astype(bool)

md = meta.set_index("assembly_accession")


def jaccard_matrix(mat):
    n = mat.shape[0]
    out = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            union = np.logical_or(mat[i], mat[j]).sum()
            inter = np.logical_and(mat[i], mat[j]).sum()
            out[i, j] = out[j, i] = 0.0 if union == 0 else 1 - inter / union
    return out


def mantel(a, b, permutations=999, seed=42):
    iu = np.triu_indices_from(a, k=1)
    x, y = a[iu], b[iu]
    r = float(np.corrcoef(x, y)[0, 1])
    rng = np.random.default_rng(seed)
    cnt = 0
    for _ in range(permutations):
        p = rng.permutation(a.shape[0])
        rp = float(np.corrcoef(x, b[np.ix_(p, p)][iu])[0, 1])
        if rp >= r:
            cnt += 1
    return r, (cnt + 1) / (permutations + 1)


def partial_mantel(a, b, c, permutations=999, seed=43):
    """r(a,b | c) with permutation of a."""
    iu = np.triu_indices_from(a, k=1)
    x, y, z = a[iu], b[iu], c[iu]
    r_ab = float(np.corrcoef(x, y)[0, 1])
    r_ac = float(np.corrcoef(x, z)[0, 1])
    r_bc = float(np.corrcoef(y, z)[0, 1])
    denom = np.sqrt(max((1 - r_ac ** 2) * (1 - r_bc ** 2), 1e-12))
    r = (r_ab - r_ac * r_bc) / denom
    rng = np.random.default_rng(seed)
    cnt = 0
    n = a.shape[0]
    for _ in range(permutations):
        p = rng.permutation(n)
        xp = a[np.ix_(p, p)][iu]
        r_ab_p = float(np.corrcoef(xp, y)[0, 1])
        r_ac_p = float(np.corrcoef(xp, z)[0, 1])
        rp = (r_ab_p - r_ac_p * r_bc) / denom
        if rp >= r:
            cnt += 1
    return r, (cnt + 1) / (permutations + 1), r_ab, r_ac, r_bc


def year_matrix(idx):
    y = md.loc[idx, "collection_year"].astype(float).values
    return np.abs(y[:, None] - y[None, :])


def country_matrix(idx):
    c = md.loc[idx, "country"].fillna("").values
    return (c[:, None] != c[None, :]).astype(float)


out = {}

# ---- Mantel and partial Mantel ----
avail = [a for a in accs if a in mash.index]
jac_all = jaccard_matrix(presence.loc[avail].values)
mas_all = mash.loc[avail, avail].values.astype(float)
r_all, p_all = mantel(jac_all, mas_all)
out["mantel_all"] = {"r": round(r_all, 3), "p": round(p_all, 4), "n": len(avail),
                     "note": "999 permutations; p=0.001 is the resolution floor"}

sub = [a for a in avail if pd.notna(md.loc[a, "collection_year"]) and str(md.loc[a, "country"]).strip() not in ("", "nan", "NA")]
if len(sub) >= 20:
    idx = [avail.index(a) for a in sub]
    jac = jac_all[np.ix_(idx, idx)]
    mas = mas_all[np.ix_(idx, idx)]
    yr = year_matrix(sub)
    co = country_matrix(sub)
    r_year, p_year = mantel(jac, yr, seed=44)
    r_country, p_country = mantel(jac, co, seed=45)
    r_partial, p_partial, r_ab, r_ac, r_bc = partial_mantel(jac, mas, yr, seed=46)
    out["mantel_year"] = {"r": round(r_year, 3), "p": round(p_year, 4), "n": len(sub)}
    out["mantel_country"] = {"r": round(r_country, 3), "p": round(p_country, 4), "n": len(sub)}
    out["partial_mantel_year"] = {"r": round(r_partial, 3), "p": round(p_partial, 4),
                                  "n": len(sub), "r_ab": round(r_ab, 3),
                                  "r_ac": round(r_ac, 3), "r_bc": round(r_bc, 3)}
    # both controls jointly: residualise year and country by regression
    def residualise(mat, controls):
        iu = np.triu_indices_from(mat, k=1)
        y = mat[iu]
        A = np.column_stack([c[iu] for c in controls] + [np.ones(len(y))])
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        res = y - A @ beta
        outm = np.zeros_like(mat)
        outm[iu] = res
        outm[(iu[1], iu[0])] = res
        return outm
    jr = residualise(jac, [yr, co])
    mr = residualise(mas, [yr, co])
    r2, p2 = mantel(jr, mr, seed=47)
    out["partial_mantel_year_country"] = {"r": round(r2, 3), "p": round(p2, 4), "n": len(sub)}
else:
    out["partial_mantel_year"] = {"r": None, "p": None, "n": len(sub)}

with open(os.path.join(RES, "mantel_partial.json"), "w") as fh:
    json.dump(out, fh, indent=2)

# ---- 630 positive control ----
pc = details[details["assembly_accession"] == "GCF_000009205.2"].copy()
pc = pc.sort_values("length_bp", ascending=False)
out630 = {
    "accession": "GCF_000009205.2",
    "strain": "630",
    "n_regions": int(len(pc)),
    "regions": [
        {"prophage_id": r["prophage_id"], "length_bp": int(r["length_bp"]),
         "cds_count": int(r["cds_count"]), "confidence": r["confidence"],
         "cluster_id": r["cluster_id"], "complete_like": int(r["complete_like"])}
        for _, r in pc.iterrows()
    ],
}
high = pc[pc["confidence"] == "high"].sort_values("length_bp", ascending=False)
if len(high) >= 2:
    out630["two_largest_high"] = {
        "first": {"prophage_id": high.iloc[0]["prophage_id"], "length_bp": int(high.iloc[0]["length_bp"]),
                  "cds_count": int(high.iloc[0]["cds_count"]), "cluster_id": high.iloc[0]["cluster_id"]},
        "second": {"prophage_id": high.iloc[1]["prophage_id"], "length_bp": int(high.iloc[1]["length_bp"]),
                   "cds_count": int(high.iloc[1]["cds_count"]), "cluster_id": high.iloc[1]["cluster_id"]},
        "same_cluster": bool(high.iloc[0]["cluster_id"] == high.iloc[1]["cluster_id"]),
    }
with open(os.path.join(RES, "positive_control_630.json"), "w") as fh:
    json.dump(out630, fh, indent=2)

# ---- reference similarity: module vs genome level ----
hits = details[details["best_pident"].notna()].copy()
hits["best_pident"] = pd.to_numeric(hits["best_pident"], errors="coerce")
hits["best_qcov"] = pd.to_numeric(hits["best_qcov"], errors="coerce")
refsim = {
    "regions": int(len(details)),
    "regions_with_hit": int(len(hits)),
    "regions_with_hit_pct": round(100 * len(hits) / len(details), 1),
    "pident_ge90": int((hits["best_pident"] >= 90).sum()),
    "pident_ge90_pct_of_hits": round(100 * float((hits["best_pident"] >= 90).mean()), 1),
    "median_pident": round(float(hits["best_pident"].median()), 1),
    "min_pident": round(float(hits["best_pident"].min()), 1),
    "max_pident": round(float(hits["best_pident"].max()), 1),
    "median_qcov_pct": round(100 * float(hits["best_qcov"].median()), 1),
    "qcov_ge50": int((hits["best_qcov"] >= 0.5).sum()),
    "qcov_ge80": int((hits["best_qcov"] >= 0.8).sum()),
    "typeable_regions": int((details["phage_typeable"] == 1).sum()),
    "typeable_clusters": int(details[details["phage_typeable"] == 1]["cluster_id"].nunique()),
}
with open(os.path.join(RES, "reference_similarity.json"), "w") as fh:
    json.dump(refsim, fh, indent=2)

# ---- defence keyword vs DefenseFinder ----
cargo_path = os.path.join(RES, "prophage_cargo_summary.csv")
if os.path.exists(cargo_path):
    cargo = pd.read_csv(cargo_path, dtype={"ST": str, "mlst_clade": str})
    high = cargo[cargo["confidence"] == "high"]
    dfo = {
        "prophages": int(len(cargo)),
        "keyword_defence": int(cargo["has_defence"].sum()),
        "df_complete_any": int((cargo["n_df_systems"] > 0).sum()),
        "df_any_hit": int((cargo["has_any_df_hit"] > 0).sum()),
        "both": int(((cargo["has_defence"] == 1) & (cargo["n_df_systems"] > 0)).sum()),
        "keyword_only": int(((cargo["has_defence"] == 1) & (cargo["n_df_systems"] == 0)).sum()),
        "df_only": int(((cargo["has_defence"] == 0) & (cargo["n_df_systems"] > 0)).sum()),
        "high_n": int(len(high)),
        "high_keyword_defence": int(high["has_defence"].sum()),
        "high_df_complete_any": int((high["n_df_systems"] > 0).sum()),
        "high_df_any_hit": int((high["has_any_df_hit"] > 0).sum()),
    }
    with open(os.path.join(RES, "defence_overlap.json"), "w") as fh:
        json.dump(dfo, fh, indent=2)

print(json.dumps({"mantel": out, "refsim": refsim, "630": out630}, indent=2))
