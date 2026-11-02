#!/usr/bin/env python3
"""Collect every number quoted in the manuscript into one JSON file."""
import json
import os

import numpy as np
import pandas as pd

BASE = "/home/bioinfo/cdprophage"
RES = os.path.join(BASE, "scaleup_241", "results")


def j(name):
    with open(os.path.join(RES, name)) as fh:
        return json.load(fh)


ds = j("diversity_summary.json")
cs = j("cargo_summary.json")
df = j("defensefinder_summary.json")
pn = j("protein_network_summary.json")
mp = j("mantel_partial.json")
refsim = j("reference_similarity.json")
pc630 = j("positive_control_630.json")
dfo = j("defence_overlap.json")

clades = pd.read_csv(os.path.join(RES, "clade_summary.csv"), dtype={"mlst_clade": str})
tests = pd.read_csv(os.path.join(RES, "statistical_tests.csv"))
meta = pd.read_csv(os.path.join(RES, "genome_metadata_full.csv"), dtype={"ST": str, "mlst_clade": str})
cargo = pd.read_csv(os.path.join(RES, "prophage_cargo_summary.csv"), dtype={"ST": str, "mlst_clade": str})
details = pd.read_csv(os.path.join(RES, "prophage_details.csv"), dtype={"ST": str, "mlst_clade": str})

V = {}
V.update(ds)
V["clade_sizes"] = ", ".join(str(int(x)) for x in clades.sort_values("mlst_clade")["n_genomes"])
V["clade_table"] = clades.to_dict(orient="records")

V["metadata"] = {
    "countries": int((meta["country"].fillna("") != "").sum()),
    "collection_year": int(meta["collection_year"].notna().sum()),
    "human": int((meta["host"] == "Homo sapiens").sum()),
    "dog": int((meta["host"] == "Dog").sum()),
    "top_countries": meta["country"].replace("", np.nan).dropna().value_counts().head(6).to_dict(),
}
V["mantel"] = mp
V["reference_similarity"] = refsim
V["positive_control_630"] = pc630
V["defence_overlap"] = dfo
V["defensefinder"] = df
V["protein_network"] = pn
V["cargo"] = cs

# cargo fractions as percentages
en = cs["enrichment_vs_chromosome"]
V["cargo_table"] = {
    k: {
        "prophage_pct": round(100 * en[k]["prophage_cds_fraction_mean"], 2),
        "chromosome_pct": round(100 * en[k]["chromosome_cds_fraction_mean"], 2),
        "q": en[k]["q_value_bh"],
        "p": en[k]["p_value"],
        "n": en[k]["n_genomes"],
    } for k in en
}
V["cargo_counts"] = {
    "cds_total": cs["cds_total"],
    "defence_cds": cs["defence_cds"],
    "virulence_cds": cs["virulence_cds"],
    "amr_cds": cs["amr_cds"],
    "mge_cds": cs["mge_cds"],
}
V["cargo_rates"] = {
    "defence_all": round(100 * cs["defence_rate_all"], 1),
    "defence_high": round(100 * cs["defence_rate_high"], 1),
    "virulence_all": round(100 * cs["virulence_rate_all"], 1),
    "virulence_high": round(100 * cs["virulence_rate_high"], 1),
    "amr_all": round(100 * cs["amr_rate_all"], 1),
    "amr_high": round(100 * cs["amr_rate_high"], 1),
    "mge_all": round(100 * cs["mge_rate_all"], 1),
    "mge_high": round(100 * cs["mge_rate_high"], 1),
}
V["defence_system_rate"] = {
    "df_complete_all": round(100 * df["prophage_complete_system_rate"], 1),
    "df_complete_high": round(100 * float((cargo[cargo["confidence"] == "high"]["n_df_systems"] > 0).mean()), 1),
    "df_any_all": round(100 * float(cargo["has_any_df_hit"].mean()), 1),
    "df_any_high": round(100 * float(cargo[cargo["confidence"] == "high"]["has_any_df_hit"].mean()), 1),
}
V["tests"] = tests.to_dict(orient="records")

hi = details[details["confidence"] == "high"]
clade_per_cluster = details.groupby("cluster_id")["mlst_clade"].nunique()
multi_clade_ids = set(clade_per_cluster[clade_per_cluster >= 2].index)
V["multi_clade_clusters"] = int(len(multi_clade_ids))
V["region_fraction_multi_clade"] = round(float(details["cluster_id"].isin(multi_clade_ids).mean()), 4)

V["high_confidence_subset"] = {
    "regions": int(len(hi)),
    "clusters": int(hi["cluster_id"].nunique()),
    "genomes": int(hi["assembly_accession"].nunique()),
    "singleton_clusters": int((hi.groupby("cluster_id")["assembly_accession"].nunique() == 1).sum()),
    "mean_per_genome": round(len(hi) / max(hi["assembly_accession"].nunique(), 1), 2),
}

try:
    V["catalogue"] = j("catalogue_curation.json")
except FileNotFoundError:
    V["catalogue"] = {}

try:
    dsys = pd.read_csv(os.path.join(RES, "defensefinder_systems.csv"))
    rm = dsys[dsys["system_name"].astype(str).str.startswith("RM")]
    V["rm_systems"] = {
        "calls": int(len(rm)),
        "complete_calls": int(rm["complete"].sum()),
        "regions": int(rm["prophage_id"].nunique()),
        "region_ids": sorted(rm["prophage_id"].unique().tolist()),
    }
except FileNotFoundError:
    V["rm_systems"] = {}

try:
    arc = pd.read_csv(os.path.join(RES, "all_regions_raw", "cluster_summary.csv"))
    V["all_region_sensitivity"] = {
        "regions": int(arc["n_prophages"].sum()),
        "clusters": int(len(arc)),
        "largest_cluster": int(arc["n_prophages"].max()),
        "largest_cluster_fraction": round(float(arc["n_prophages"].max() / arc["n_prophages"].sum()), 3),
    }
except FileNotFoundError:
    V["all_region_sensitivity"] = {}

with open(os.path.join(RES, "manuscript_numbers.json"), "w") as fh:
    json.dump(V, fh, indent=2, default=str)

print(json.dumps(V, indent=2, default=str))
