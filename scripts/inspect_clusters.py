#!/usr/bin/env python3
import pandas as pd

p = "/home/bioinfo/cdprophage/scaleup_241/results/prophage_details.csv"
d = pd.read_csv(p, dtype={"ST": str, "mlst_clade": str})
print(d.columns.tolist())
for c in ["terminase_n", "capsid_n", "tail_n", "lysin_n", "integrase_n"]:
    print(c, int((d[c] > 0).sum()))
d["structural"] = d[["terminase_n", "capsid_n", "tail_n", "lysin_n"]].sum(axis=1) > 0
d["any_hall"] = d[["terminase_n", "capsid_n", "tail_n", "lysin_n", "integrase_n"]].sum(axis=1) > 0
print("structural hallmark", int(d.structural.sum()), "any hallmark", int(d.any_hall.sum()),
      "none", int((~d.any_hall).sum()))
print(d.groupby("confidence").agg(n=("prophage_id", "size"), structural=("structural", "sum"),
                                  med_len=("length_bp", "median")).to_string())
print("low no structural", int(((d.confidence == "low") & (~d.structural)).sum()))
print(d[~d.structural].length_bp.describe())
print(d[d.confidence == "low"].hallmark_categories.value_counts().head(10))
