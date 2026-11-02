#!/usr/bin/env python3
"""Define the primary prophage catalogue.

A region enters the primary catalogue when it encodes a packaging module
(terminase or capsid). PhiSpy also reports tail-only contractile elements and
integrase-only islands; these lack a packaging module and are retained in
prophage_details_all.csv but excluded from the primary catalogue.
"""
import csv
import json
import os

from Bio import SeqIO

BASE = "/home/bioinfo/cdprophage"
RESULT = os.path.join(BASE, "scaleup_241", "results")

with open(os.path.join(RESULT, "prophage_details_raw.csv")) as fh:
    rows = list(csv.DictReader(fh))

for r in rows:
    r["packaging_module"] = int(int(r["terminase_n"]) > 0 or int(r["capsid_n"]) > 0)
    r["in_primary"] = r["packaging_module"]

primary = [r for r in rows if r["in_primary"] == 1]
primary_ids = {r["prophage_id"] for r in primary}

with open(os.path.join(RESULT, "prophage_details_all.csv"), "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

with open(os.path.join(RESULT, "prophage_details_primary_raw.csv"), "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=[k for k in rows[0] if k not in ("packaging_module", "in_primary")])
    writer.writeheader()
    for r in primary:
        writer.writerow({k: v for k, v in r.items() if k not in ("packaging_module", "in_primary")})

n_seq = 0
with open(os.path.join(RESULT, "prophage_sequences_all.fasta")) as fh, \
        open(os.path.join(RESULT, "prophage_sequences.fasta"), "w") as out:
    for rec in SeqIO.parse(fh, "fasta"):
        if rec.id in primary_ids:
            SeqIO.write(rec, out, "fasta")
            n_seq += 1

summary = {
    "regions_phispy": len(rows),
    "regions_primary": len(primary),
    "regions_excluded": len(rows) - len(primary),
    "excluded_with_terminase": sum(1 for r in rows if r["in_primary"] == 0 and int(r["terminase_n"]) > 0),
    "excluded_with_capsid": sum(1 for r in rows if r["in_primary"] == 0 and int(r["capsid_n"]) > 0),
    "excluded_tail_only": sum(1 for r in rows if r["in_primary"] == 0 and int(r["tail_n"]) > 0),
    "excluded_no_structural": sum(1 for r in rows if r["in_primary"] == 0
                                  and int(r["terminase_n"]) == 0 and int(r["capsid_n"]) == 0
                                  and int(r["tail_n"]) == 0 and int(r["lysin_n"]) == 0),
    "excluded_lysin_only": sum(1 for r in rows if r["in_primary"] == 0
                               and int(r["terminase_n"]) == 0 and int(r["capsid_n"]) == 0
                               and int(r["tail_n"]) == 0 and int(r["lysin_n"]) > 0),
    "primary_high": sum(1 for r in primary if r["confidence"] == "high"),
    "primary_medium": sum(1 for r in primary if r["confidence"] == "medium"),
    "primary_low": sum(1 for r in primary if r["confidence"] == "low"),
    "primary_genomes": len({r["assembly_accession"] for r in primary}),
    "sequences_written": n_seq,
}
with open(os.path.join(RESULT, "catalogue_curation.json"), "w") as fh:
    json.dump(summary, fh, indent=2)
print(json.dumps(summary, indent=2))
