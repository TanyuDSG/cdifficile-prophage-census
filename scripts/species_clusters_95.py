#!/usr/bin/env python3
"""Species-level demarcation of primary-catalogue prophages.

Applies a >=95% nucleotide identity / >=50% reciprocal coverage criterion to
the all-vs-all BLASTn output already generated for the 70% genus-level
clustering.  This is the operational equivalent of the VIRIDIC species
threshold used by Holtappels et al. (2026); it is not a VIRIDIC run.
"""
import collections
import json
import os

import pandas as pd

BASE = "/home/bioinfo/cdprophage"
RES = os.path.join(BASE, "scaleup_241", "results")
BLAST = os.path.join(RES, "prophage_primary_blast.tsv")
DETAILS = os.path.join(RES, "prophage_details_all.csv")
OUT_CSV = os.path.join(RES, "species_clusters_95.csv")
OUT_JSON = os.path.join(RES, "species_clusters_95_summary.json")

MIN_IDENT = 95.0
MIN_COV = 0.5

# best HSP per ordered pair, reciprocal-coverage filtered
best = {}
with open(BLAST) as fh:
    for line in fh:
        q, s, pident, alen, qlen, slen, _evalue, _bits = line.split()[:8]
        if q == s:
            continue
        pident = float(pident)
        alen, qlen, slen = int(alen), int(qlen), int(slen)
        if alen / qlen < MIN_COV or alen / slen < MIN_COV:
            continue
        score = alen * pident
        if (q, s) not in best or score > best[(q, s)][0]:
            best[(q, s)] = (score, pident, alen, qlen, slen)

parent = {}


def find(x):
    parent.setdefault(x, x)
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(a, b):
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[rb] = ra


for (q, s), (_score, pident, _alen, _qlen, _slen) in best.items():
    if pident >= MIN_IDENT:
        union(q, s)

# include every catalogue region, even if it had no qualifying pair
details = pd.read_csv(DETAILS)
catalogue = details.loc[details["in_primary"] == 1, ["prophage_id", "length_bp"]].copy()
for pid in catalogue["prophage_id"]:
    find(pid)

clusters = collections.defaultdict(list)
for pid in catalogue["prophage_id"]:
    clusters[find(pid)].append(pid)

# representative = longest region in cluster
length = dict(zip(catalogue["prophage_id"], catalogue["length_bp"]))
rows = []
for members in clusters.values():
    rep = max(members, key=lambda p: length.get(p, 0))
    cid = f"SC{len(rows) + 1:03d}"
    for pid in members:
        rows.append({"prophage_id": pid, "species_cluster_95": cid,
                     "cluster_size": len(members), "representative": rep})
rows.sort(key=lambda r: (r["species_cluster_95"], r["prophage_id"]))
out = pd.DataFrame(rows)
out.to_csv(OUT_CSV, index=False)

sizes = sorted((len(m) for m in clusters.values()), reverse=True)
summary = {
    "criterion": ">=95% nucleotide identity over >=50% of both sequences",
    "n_regions": int(len(catalogue)),
    "n_species_clusters": int(len(clusters)),
    "largest_cluster": int(sizes[0]) if sizes else 0,
    "singletons": int(sum(1 for s in sizes if s == 1)),
    "clusters_ge5": int(sum(1 for s in sizes if s >= 5)),
    "n_pairs_ge95": int(sum(1 for v in best.values() if v[1] >= MIN_IDENT)),
    "n_pairs_qualifying": int(len(best)),
    "size_distribution": sizes,
    "note": "BLASTn operational criterion; not a VIRIDIC run.",
}
with open(OUT_JSON, "w") as fh:
    json.dump(summary, fh, indent=2)
print(json.dumps({k: v for k, v in summary.items() if k != "size_distribution"}, indent=2))
