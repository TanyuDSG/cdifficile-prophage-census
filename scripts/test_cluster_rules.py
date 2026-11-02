#!/usr/bin/env python3
"""Compare clustering rules on the prophage catalogue."""
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

length = dict(zip(det.prophage_id, det.length_bp))
confidence = dict(zip(det.prophage_id, det.confidence))
structural = dict(zip(det.prophage_id, det.structural))
clade = dict(zip(det.prophage_id, det.mlst_clade))


def run(nodes, min_id=70, both=True, min_len=0):
    nodes = set(nodes)
    d = blast[(blast.pid >= min_id) & blast.q.isin(nodes) & blast.s.isin(nodes)].copy()
    d = d[(d.cq >= 0.5) & (d.cs >= 0.5)] if both else d[(d.cq >= 0.5) | (d.cs >= 0.5)]
    if min_len:
        d = d[[length[q] >= min_len and length[s] >= min_len for q, s in zip(d.q, d.s)]]
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
    for x in parent:
        groups[find(x)].append(x)
    for x in nodes - set(parent):
        groups[x].append(x)
    sizes = sorted((len(v) for v in groups.values()), reverse=True)
    top = max(groups.values(), key=len)
    return len(groups), sizes[:6], len(top), \
        round(sum(confidence[x] == "high" for x in top) / len(top), 2), \
        round(sum(structural[x] for x in top) / len(top), 2)


print("ALL regions, shorter coverage:", run(det.prophage_id, 70, both=False))
print("ALL regions, reciprocal coverage:", run(det.prophage_id, 70, both=True))
nodes = list(det[det.structural].prophage_id)
print("structural+, reciprocal:", run(nodes, 70, both=True))
print("structural+, shorter:", run(nodes, 70, both=False))
print("structural+, reciprocal, >=8kb:", run(nodes, 70, both=True, min_len=8000))
print("structural+, reciprocal, >=10kb:", run(nodes, 70, both=True, min_len=10000))
hi = list(det[det.confidence == "high"].prophage_id)
print("high only, reciprocal:", run(hi, 70, both=True))
print("high only, shorter:", run(hi, 70, both=False))
