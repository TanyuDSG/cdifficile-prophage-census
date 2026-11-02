#!/usr/bin/env python3
"""Merge BioSample geography/date/source into the genome metadata table."""
import glob
import json
import os

import pandas as pd

BASE = "/home/bioinfo/cdprophage"
META = os.path.join(BASE, "genome_metadata.csv")
OUT = os.path.join(BASE, "scaleup_241", "results", "genome_metadata_full.csv")


def get(d, *path, default=""):
    cur = d
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    if isinstance(cur, list):
        cur = ";".join(str(x) for x in cur)
    return "" if cur is None else str(cur)


rows = []
for path in sorted(glob.glob(os.path.join(BASE, "genomes", "metadata", "*.json"))):
    acc = os.path.basename(path).replace(".json", "")
    with open(path, encoding="utf-8") as fh:
        obj = json.load(fh)
    bs = obj.get("assemblyInfo", {}).get("biosample", {})
    rows.append({
        "assembly_accession": acc,
        "geo_loc_name": get(bs, "geoLocName"),
        "collection_date": get(bs, "collectionDate"),
        "host": get(bs, "host"),
        "isolation_source": get(bs, "isolationSource"),
        "biosample": get(bs, "accession"),
    })
extra = pd.DataFrame(rows)

meta = pd.read_csv(META, dtype={"ST": str, "mlst_clade": str})
meta["assembly_accession"] = meta["assembly_accession"].astype(str)
merged = meta.drop(columns=[c for c in ("country", "collection_date", "host", "isolation_source")
                            if c in meta.columns]).merge(
    extra, on="assembly_accession", how="left")
merged["country"] = merged["geo_loc_name"].fillna("").str.split(":").str[0].str.strip()
merged["year"] = pd.to_datetime(merged["release_date"], errors="coerce").dt.year
merged["collection_year"] = pd.to_numeric(
    merged["collection_date"].astype(str).str.extract(r"(\d{4})")[0], errors="coerce")
merged.to_csv(OUT, index=False)

print("rows", len(merged))
print("with country", int((merged["country"] != "").sum()))
print("with collection year", int(merged["collection_year"].notna().sum()))
print(merged["country"].replace("", "NA").value_counts().head(12).to_string())
print(merged["collection_year"].describe().to_string())
