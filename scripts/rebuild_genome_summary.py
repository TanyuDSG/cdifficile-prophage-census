#!/usr/bin/env python3
"""Recompute per-genome prophage totals from the curated catalogue.

annotate_prophages_all.py writes genome totals for every PhiSpy call; after
curation only packaging-positive regions belong to the catalogue, so the
per-genome totals have to be rebuilt from prophage_details.csv.
"""
import os

import pandas as pd

BASE = "/home/bioinfo/cdprophage"
RES = os.path.join(BASE, "scaleup_241", "results")

details = pd.read_csv(os.path.join(RES, "prophage_details.csv"), dtype={"ST": str, "mlst_clade": str})
base = pd.read_csv(os.path.join(RES, "genome_prophage_summary_raw.csv"),
                   dtype={"ST": str, "mlst_clade": str})

agg = details.groupby("assembly_accession").agg(
    n_prophages=("prophage_id", "nunique"),
    prophage_bp=("length_bp", "sum"),
)
out = base.drop(columns=["n_prophages", "prophage_bp", "prophage_percent"]).merge(
    agg, on="assembly_accession", how="left")
out["n_prophages"] = out["n_prophages"].fillna(0).astype(int)
out["prophage_bp"] = out["prophage_bp"].fillna(0).astype(int)
out["prophage_percent"] = (100 * out["prophage_bp"] / out["genome_size"]).round(3)
out = out.sort_values(["mlst_clade", "ST", "assembly_accession"])
out.to_csv(os.path.join(RES, "genome_prophage_summary_raw.csv"), index=False)
print(out[["n_prophages", "prophage_percent"]].describe().to_string())
print("genomes", len(out), "regions", int(out["n_prophages"].sum()),
      "carriage", round(float((out["n_prophages"] > 0).mean()), 4))
