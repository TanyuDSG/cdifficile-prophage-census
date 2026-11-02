#!/usr/bin/env python3
"""Characterise the 216-member structural-hallmark cluster."""
from collections import Counter, defaultdict

import pandas as pd

B = "/home/bioinfo/cdprophage/scaleup_241/results/"
det = pd.read_csv(B + "prophage_details.csv", dtype={"ST": str, "mlst_clade": str})
det["structural"] = det[["terminase_n", "capsid_n", "tail_n", "lysin_n"]].sum(axis=1) > 0
blast = pd.read_csv(B + "prophage_self_blast.tsv", sep="\t", header=None,
                    names=["q", "s", "pid", "alen", "qlen", "slen", "evalue", "bits"])
blast = blast[blast.q != blast.s]
blast["cq"] = blast.alen / blast.qlen
blast["cs"] = blast.alen / blast.slen
nodes = set(det[det.structural].prophage_id)
d = blast[(blast.pid >= 70) & blast.q.isin(nodes) & blast.s.isin(nodes)]
d = d[(d.cq >= 0.5) & (d.cs >= 0.5)]
parent = {}

def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

for q, s in zip(d.q, d.s):
    rq, rs = find(q), find(s)
    if rq != rs:
        parent[rq] = rs
groups = defaultdict(list)
for x in nodes:
    groups[find(x)].append(x)
big = max(groups.values(), key=len)
sub = det.set_index("prophage_id").loc[big]
print("cluster size", len(sub))
print(sub[["length_bp", "cds_count", "gc_percent", "gc_delta", "terminase_n", "capsid_n",
           "tail_n", "lysin_n", "integrase_n", "best_pident", "best_qcov"]].describe().round(2).to_string())
print("\nhallmark categories")
print(sub.hallmark_categories.value_counts().head(10).to_string())
print("\nbest phage")
print(sub.best_phage.value_counts().head(10).to_string())
print("\nclades", sub.mlst_clade.value_counts().to_dict())
print("\nstrain 630 regions")
print(det[det.assembly_accession == "GCF_000009205.2"][
    ["prophage_id", "start", "stop", "length_bp", "confidence", "hallmark_categories",
     "cluster_id", "best_phage"]].to_string())
print("\nbig cluster example ids", big[:5])
