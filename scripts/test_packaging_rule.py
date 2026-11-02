#!/usr/bin/env python3
"""Evaluate a packaging-module definition of prophage regions."""
from collections import defaultdict

import pandas as pd

B = "/home/bioinfo/cdprophage/scaleup_241/results/"
det = pd.read_csv(B + "prophage_details.csv", dtype={"ST": str, "mlst_clade": str})
det["packaging"] = (det.terminase_n > 0) | (det.capsid_n > 0)
print("terminase", int((det.terminase_n > 0).sum()), "capsid", int((det.capsid_n > 0).sum()),
      "packaging", int(det.packaging.sum()), "fraction", round(float(det.packaging.mean()), 3))
print(det[det.packaging].confidence.value_counts().to_dict())
print("median length packaging", int(det[det.packaging].length_bp.median()))
print("per genome", round(float(det.packaging.sum()) / det.assembly_accession.nunique(), 2))
print("genomes with packaging", det[det.packaging].assembly_accession.nunique())

blast = pd.read_csv(B + "prophage_self_blast.tsv", sep="\t", header=None,
                    names=["q", "s", "pid", "alen", "qlen", "slen", "evalue", "bits"])
blast = blast[blast.q != blast.s]
blast["cq"] = blast.alen / blast.qlen
blast["cs"] = blast.alen / blast.slen
for both in (True, False):
    nodes = set(det[det.packaging].prophage_id)
    d = blast[(blast.pid >= 70) & blast.q.isin(nodes) & blast.s.isin(nodes)]
    d = d[(d.cq >= 0.5) & (d.cs >= 0.5)] if both else d[(d.cq >= 0.5) | (d.cs >= 0.5)]
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
    sizes = sorted((len(v) for v in groups.values()), reverse=True)
    big = max(groups.values(), key=len)
    sub = det.set_index("prophage_id").loc[big]
    print(("reciprocal" if both else "shorter"), "clusters", len(groups), "top", sizes[:8],
          "big high frac", round(float((sub.confidence == "high").mean()), 2),
          "big med len", int(sub.length_bp.median()), "big clades", sub.mlst_clade.nunique())
