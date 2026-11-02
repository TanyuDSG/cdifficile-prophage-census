#!/usr/bin/env python3
"""Cluster primary prophage sequences by BLAST identity/coverage (single linkage).

Criterion: >=70% nucleotide identity and >=50% coverage of both sequences.
Requiring reciprocal coverage stops short fragments from bridging otherwise
distinct elements; a shorter-sequence rule collapsed two-thirds of the raw
PhiSpy predictions into one cluster.
"""
import csv
import os
import subprocess
from collections import defaultdict

from Bio import SeqIO

BASE = "/home/bioinfo/cdprophage"
RESULT = os.path.join(BASE, "scaleup_241", "results")
TOOLS = os.path.join(BASE, "tools", "root")
BLASTN = os.path.join(TOOLS, "usr", "bin", "blastn")
MAKEBLASTDB = os.path.join(TOOLS, "usr", "bin", "makeblastdb")
os.environ["LD_LIBRARY_PATH"] = os.pathsep.join([
    os.path.join(TOOLS, "usr", "lib", "ncbi-blast+"),
    os.path.join(TOOLS, "usr", "lib", "x86_64-linux-gnu"),
])

FASTA = os.path.join(RESULT, "prophage_sequences.fasta")
DB = os.path.join(RESULT, "prophage_primary_db")
HITS = os.path.join(RESULT, "prophage_primary_blast.tsv")


def cluster_fasta(fasta, db_prefix, hits_path, min_pident=70.0, min_cov=0.5):
    subprocess.run([MAKEBLASTDB, "-in", fasta, "-dbtype", "nucl", "-out", db_prefix],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run([BLASTN, "-query", fasta, "-db", db_prefix, "-outfmt",
                    "6 qseqid sseqid pident length qlen slen evalue bitscore",
                    "-evalue", "1e-10", "-max_target_seqs", "10000",
                    "-num_threads", "8", "-out", hits_path],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse(fasta, "fasta")}
    parent = {acc: acc for acc in lengths}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    edges = []
    with open(hits_path) as fh:
        for line in fh:
            q, s, pident, alen, qlen, slen = line.split()[:6]
            if q == s:
                continue
            pident = float(pident); alen = int(alen); qlen = int(qlen); slen = int(slen)
            cov_q = alen / qlen if qlen else 0
            cov_s = alen / slen if slen else 0
            cov = min(cov_q, cov_s)
            if pident >= min_pident and cov_q >= min_cov and cov_s >= min_cov:
                union(q, s)
                edges.append((q, s, pident, cov))
    groups = defaultdict(list)
    for acc in lengths:
        groups[find(acc)].append(acc)
    clusters = []
    for members in groups.values():
        members.sort(key=lambda x: lengths[x], reverse=True)
        clusters.append(members)
    clusters.sort(key=lambda m: (-len(m), -lengths[m[0]], m[0]))
    cluster_of = {}
    for i, members in enumerate(clusters, 1):
        for acc in members:
            cluster_of[acc] = f"PC{i:04d}"
    return clusters, cluster_of, edges


clusters, cluster_of, edges = cluster_fasta(FASTA, DB, HITS)

with open(os.path.join(RESULT, "prophage_clusters.csv"), "w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(["cluster_id", "cluster_size", "representative", "members"])
    for i, members in enumerate(clusters, 1):
        writer.writerow([f"PC{i:04d}", len(members), members[0], ";".join(members)])

with open(os.path.join(RESULT, "prophage_details_primary_raw.csv")) as fh:
    rows = list(csv.DictReader(fh))
for row in rows:
    row["cluster_id"] = cluster_of.get(row["prophage_id"], "")
with open(os.path.join(RESULT, "prophage_details.csv"), "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

by_cluster = defaultdict(list)
for row in rows:
    by_cluster[row["cluster_id"]].append(row)
with open(os.path.join(RESULT, "cluster_summary.csv"), "w", newline="") as fh:
    fields = ["cluster_id", "n_prophages", "n_genomes", "n_STs", "n_clades",
              "median_length_bp", "min_length_bp", "max_length_bp",
              "mean_gc_percent", "n_high", "n_complete_like", "n_typeable",
              "clades", "STs", "best_phage", "representative"]
    writer = csv.DictWriter(fh, fieldnames=fields)
    writer.writeheader()
    for cid in sorted(by_cluster):
        items = by_cluster[cid]
        lengths = sorted(int(r["length_bp"]) for r in items)
        med = lengths[len(lengths) // 2] if len(lengths) % 2 else (lengths[len(lengths) // 2 - 1] + lengths[len(lengths) // 2]) / 2
        writer.writerow({
            "cluster_id": cid,
            "n_prophages": len(items),
            "n_genomes": len({r["assembly_accession"] for r in items}),
            "n_STs": len({r["ST"] for r in items}),
            "n_clades": len({r["mlst_clade"] for r in items}),
            "median_length_bp": int(med),
            "min_length_bp": min(lengths),
            "max_length_bp": max(lengths),
            "mean_gc_percent": round(sum(float(r["gc_percent"]) for r in items) / len(items), 2),
            "n_high": sum(1 for r in items if r["confidence"] == "high"),
            "n_complete_like": sum(int(r["complete_like"]) for r in items),
            "n_typeable": sum(int(r["phage_typeable"]) for r in items),
            "clades": ";".join(sorted({r["mlst_clade"] for r in items})),
            "STs": ";".join(sorted({r["ST"] for r in items}, key=lambda x: int(x) if x.isdigit() else 0)),
            "best_phage": items[0]["best_phage"],
            "representative": items[0]["prophage_id"],
        })

singletons = sum(1 for m in clusters if len(m) == 1)
print(f"clusters={len(clusters)} prophages={len(rows)} singletons={singletons} edges={len(edges)}")
