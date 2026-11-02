#!/usr/bin/env python3
"""Summarise BUSCO v5.4.7 genome-mode results across the 241 genomes."""
import glob
import json
import os

import pandas as pd

BASE = "/home/bioinfo/cdprophage"
OUTDIR = os.path.join(BASE, "scaleup_241", "busco_all")
RES = os.path.join(BASE, "scaleup_241", "results")

rows = []
bad = []
for path in glob.glob(os.path.join(OUTDIR, "*", "short_summary*.json")):
    try:
        with open(path) as fh:
            d = json.load(fh)
    except (json.JSONDecodeError, OSError):
        bad.append(path)
        continue
    r = d["results"]
    acc = os.path.basename(os.path.dirname(path))
    rows.append({
        "assembly_accession": acc,
        "complete_pct": r["Complete"],
        "single_copy_pct": r["Single copy"],
        "duplicated_pct": r["Multi copy"],
        "fragmented_pct": r["Fragmented"],
        "missing_pct": r["Missing"],
        "n_markers": r["n_markers"],
        "lineage": d["lineage_dataset"]["name"],
        "busco_version": d["versions"]["busco"],
    })

df = pd.DataFrame(rows).sort_values("assembly_accession")
out_csv = os.path.join(RES, "busco_summary.csv")
df.to_csv(out_csv, index=False)

compl = df["complete_pct"]
summary = {
    "n_genomes": int(len(df)),
    "lineage": sorted(df["lineage"].unique().tolist()),
    "busco_version": sorted(df["busco_version"].unique().tolist()),
    "complete_min": float(compl.min()),
    "complete_median": float(compl.median()),
    "complete_max": float(compl.max()),
    "n_below_95": int((compl < 95.0).sum()),
    "n_below_99": int((compl < 99.0).sum()),
    "n_markers": int(df["n_markers"].iloc[0]) if len(df) else 0,
}
with open(os.path.join(RES, "busco_summary.json"), "w") as fh:
    json.dump(summary, fh, indent=2)
print(json.dumps(summary, indent=2))
if bad:
    print("UNPARSABLE SUMMARY FILES:")
    for path in bad:
        print(" ", path)
