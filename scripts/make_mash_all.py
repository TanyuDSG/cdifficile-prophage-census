#!/usr/bin/env python3
"""Mash whole-genome distances and NJ tree for all 241 genomes."""
import os
import subprocess

import numpy as np
import pandas as pd
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo

BASE = "/home/bioinfo/cdprophage"
RESULT = os.path.join(BASE, "scaleup_241", "results")
MASH = os.path.join(BASE, "tools", "root", "usr", "bin", "mash")
FASTA = os.path.join(BASE, "genomes", "fasta")
OUT = os.path.join(BASE, "scaleup_241", "mash")
os.makedirs(OUT, exist_ok=True)
os.environ["LD_LIBRARY_PATH"] = (
    os.path.join(BASE, "tools", "mashlibs", "root", "usr", "lib", "x86_64-linux-gnu")
    + ":" + os.path.join(BASE, "tools", "root", "usr", "lib", "x86_64-linux-gnu")
)

accs = sorted(pd.read_csv(os.path.join(BASE, "genome_metadata.csv"))["assembly_accession"].astype(str))
fastas = [os.path.join(FASTA, a + ".fna") for a in accs]
missing = [f for f in fastas if not os.path.exists(f)]
if missing:
    raise SystemExit("missing fasta: %s" % missing[:5])

sketch = os.path.join(OUT, "genomes.msh")
subprocess.run([MASH, "sketch", "-p", "8", "-k", "21", "-s", "10000", "-o", sketch] + fastas, check=True)

dist_tsv = os.path.join(OUT, "mash_dist.tsv")
with open(dist_tsv, "w") as fh:
    subprocess.run([MASH, "dist", "-p", "8", sketch, sketch], check=True, stdout=fh)

rows = []
with open(dist_tsv) as fh:
    for line in fh:
        ref, qry, d, pval, shared = line.split()[:5]
        ref = os.path.basename(ref).replace(".fna", "")
        qry = os.path.basename(qry).replace(".fna", "")
        rows.append((ref, qry, float(d), float(pval), int(str(shared).split("/")[0])))
long = pd.DataFrame(rows, columns=["ref", "query", "distance", "p_value", "shared_hashes"])
long.to_csv(os.path.join(OUT, "mash_dist_long.csv"), index=False)

mat = long.pivot(index="ref", columns="query", values="distance").reindex(index=accs, columns=accs)
mat.to_csv(os.path.join(RESULT, "mash_dist_matrix.csv"))

# NJ tree
lower = []
for i, a in enumerate(accs):
    row = [0.0 if j == i else float(mat.loc[a, accs[j]]) for j in range(i + 1)]
    lower.append(row)
dm = DistanceMatrix(accs, lower)
tree = DistanceTreeConstructor().nj(dm)
try:
    tree.root_at_midpoint()
except Exception as exc:
    print("midpoint rooting failed:", exc)
Phylo.write(tree, os.path.join(RESULT, "mash_nj_tree.nwk"), "newick")

print(f"genomes={len(accs)} distances={len(long)} median={long['distance'].median():.4f}")
