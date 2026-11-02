#!/usr/bin/env python3
import pandas as pd
from collections import defaultdict

B = "/home/bioinfo/cdprophage/scaleup_241/results/"
det = pd.read_csv(B + "prophage_details.csv", dtype={"ST": str, "mlst_clade": str})
det["packaging"] = (det.terminase_n > 0) | (det.capsid_n > 0)
blast = pd.read_csv(B + "prophage_self_blast.tsv", sep="\t", header=None,
                    names=["q", "s", "pid", "alen", "qlen", "slen", "evalue", "bits"])
blast = blast[blast.q != blast.s]
blast["cq"] = blast.alen / blast.qlen
blast["cs"] = blast.alen / blast.slen
nodes = set(det[det.packaging].prophage_id)
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
of = {}
for i, (k, v) in enumerate(sorted(groups.items(), key=lambda kv: -len(kv[1])), 1):
    for x in v:
        of[x] = f"PC{i:04d}"
p630 = det[det.assembly_accession == "GCF_000009205.2"].copy()
p630["new_cluster"] = p630.prophage_id.map(of)
print(p630[["prophage_id", "length_bp", "cds_count", "confidence", "terminase_n", "capsid_n",
            "tail_n", "lysin_n", "cluster_id", "new_cluster"]].to_string())
print("same cluster pp1/pp5:", of.get("GCF_000009205.2_pp1") == of.get("GCF_000009205.2_pp5"))
print("cluster sizes", len(of), max(len(v) for v in groups.values()))
