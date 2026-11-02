#!/usr/bin/env python3
"""Extract PhiSpy prophages for all 241 complete C. difficile genomes.

Inputs
  genome_metadata.csv                    241 accessions with ST/clade
  genomes/fasta/{acc}.fna                genome nucleotide sequences
  genomes/gbff/{acc}.gbff                annotated GenBank (CDS products)
  phispy/{acc}/prophage_coordinates.tsv  PhiSpy 5.0.10 predictions
  refs/phages/phage_refs                 complete C. difficile phage genomes
Outputs (scaleup_241/results/)
  prophage_details_raw.csv
  prophage_sequences_all.fasta
  genome_prophage_summary_raw.csv
  prophage_annotation_qc.json
"""
import csv
import json
import os
import re
import subprocess
from collections import Counter, defaultdict

from Bio import SeqIO

BASE = "/home/bioinfo/cdprophage"
RESULT_DIR = os.path.join(BASE, "scaleup_241", "results")
FASTA_DIR = os.path.join(BASE, "genomes", "fasta")
GBFF_DIR = os.path.join(BASE, "genomes", "gbff")
PHISPY_DIR = os.path.join(BASE, "phispy")
PHAGE_DB = os.path.join(BASE, "refs", "phages", "phage_refs")
os.makedirs(RESULT_DIR, exist_ok=True)

TOOLS = os.path.join(BASE, "tools", "root")
BLASTN = os.path.join(TOOLS, "usr", "bin", "blastn")
os.environ["LD_LIBRARY_PATH"] = os.pathsep.join([
    os.path.join(TOOLS, "usr", "lib", "ncbi-blast+"),
    os.path.join(TOOLS, "usr", "lib", "x86_64-linux-gnu"),
])

HALLMARK = {
    "terminase": re.compile(r"terminase|\bter[lLsS]\b|dna packaging", re.I),
    "capsid": re.compile(r"capsid|portal protein|head protein|major head|prohead", re.I),
    "tail": re.compile(r"\btail\b|tail fiber|tail fibre|tape measure|baseplate|tail sheath|tail tube|neck protein|gp[0-9]+.*tail", re.I),
    "lysin": re.compile(r"lysin|endolysin|holin|amidase|muramidase|autolysin|spanin", re.I),
    "integrase": re.compile(r"integrase|recombinase|\bxerc?\b|\bxerd?\b|excisionase|resolvase", re.I),
    "regulation": re.compile(r"repressor|antirepressor|transcriptional regulator|\bcro\b|excisionase", re.I),
    "replication": re.compile(r"replicat|helicase|primase|dna polymerase", re.I),
}


def read_genome(acc):
    records = {}
    order = []
    with open(os.path.join(FASTA_DIR, acc + ".fna")) as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            records[rec.id] = rec
            order.append(rec.id)
    return records, order


def read_gbff(acc):
    path = os.path.join(GBFF_DIR, acc + ".gbff")
    cds_by_contig = defaultdict(list)
    if not os.path.exists(path):
        return cds_by_contig
    for rec in SeqIO.parse(path, "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            qual = feat.qualifiers
            product = " ".join(qual.get("product", []))
            gene = " ".join(qual.get("gene", []))
            note = " ".join(qual.get("note", []))
            locus = " ".join(qual.get("locus_tag", []))
            text = " ".join([product, gene, note, locus])
            cds_by_contig[rec.id].append({
                "start": int(feat.location.start),
                "end": int(feat.location.end),
                "strand": feat.location.strand,
                "product": product,
                "gene": gene,
                "note": note,
                "locus_tag": locus,
                "text": text,
            })
    return cds_by_contig


def read_coords(acc):
    path = os.path.join(PHISPY_DIR, acc, "prophage_coordinates.tsv")
    if not os.path.exists(path):
        return []
    out = []
    with open(path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            try:
                start = int(fields[2])
                stop = int(fields[3])
            except ValueError:
                continue
            if stop < start:
                start, stop = stop, start
            out.append({"pp": fields[0], "contig": fields[1],
                        "start": start, "stop": stop,
                        "attL": fields[8] if len(fields) > 8 else "",
                        "attR": fields[9] if len(fields) > 9 else ""})
    return out


def classify(cds_list):
    counts = Counter()
    for cds in cds_list:
        for name, pattern in HALLMARK.items():
            if pattern.search(cds["text"]):
                counts[name] += 1
    return counts


def confidence(length, counts):
    cats = [name for name in HALLMARK if counts.get(name, 0) > 0]
    if length >= 15000 and counts["terminase"] and counts["capsid"] and (counts["tail"] or counts["lysin"]):
        return "high", cats
    if length >= 10000 and len(cats) >= 2 and (counts["terminase"] or counts["capsid"]):
        return "medium", cats
    return "low", cats


def gc_fraction(seq):
    seq = str(seq).upper()
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def main():
    with open(os.path.join(BASE, "genome_metadata.csv")) as fh:
        meta = list(csv.DictReader(fh))

    fasta_out = os.path.join(RESULT_DIR, "prophage_sequences_all.fasta")
    rows = []
    qc = {"genomes_total": len(meta), "genomes_with_gbff": 0, "genomes_with_phispy": 0,
          "genomes_with_prophage": 0, "missing_fasta": [], "missing_gbff": [],
          "missing_phispy": [], "unmatched_contigs": []}
    genome_rows = []
    handle = open(fasta_out, "w")
    for rep in meta:
        acc = rep["assembly_accession"]
        genome_records, order = read_genome(acc)
        cds_by_contig = read_gbff(acc)
        if cds_by_contig:
            qc["genomes_with_gbff"] += 1
        else:
            qc["missing_gbff"].append(acc)
        coords = read_coords(acc)
        if coords:
            qc["genomes_with_phispy"] += 1
        else:
            qc["missing_phispy"].append(acc)
        genome_len = sum(len(rec.seq) for rec in genome_records.values())
        genome_gc = gc_fraction("".join(str(rec.seq) for rec in genome_records.values()))
        pp_bp = 0
        for i, pp in enumerate(coords, 1):
            contig = pp["contig"]
            rec = genome_records.get(contig)
            if rec is None:
                for key in genome_records:
                    if key.split(".")[0] == contig.split(".")[0]:
                        rec = genome_records[key]
                        contig = key
                        break
            if rec is None:
                qc["unmatched_contigs"].append(f"{acc}:{pp['contig']}")
                continue
            start, stop = pp["start"], pp["stop"]
            start = max(0, start)
            stop = min(len(rec.seq), stop)
            seq = rec.seq[start:stop]
            length = len(seq)
            cds_list = [
                c for c in cds_by_contig.get(contig, [])
                if c["end"] > start and c["start"] < stop
            ]
            counts = classify(cds_list)
            level, cats = confidence(length, counts)
            complete_like = bool(
                length >= 20000 and counts["terminase"] and counts["capsid"]
                and counts["tail"] and counts["lysin"]
            )
            name = f"{acc}_pp{i}"
            handle.write(f">{name} {acc} {contig}:{start}-{stop} len={length} gc={gc_fraction(seq):.3f}\n")
            for j in range(0, len(seq), 80):
                handle.write(str(seq[j:j + 80]) + "\n")
            pp_bp += length
            rows.append({
                "prophage_id": name,
                "assembly_accession": acc,
                "strain": rep.get("strain", ""),
                "ST": rep.get("ST", ""),
                "mlst_clade": rep.get("mlst_clade", ""),
                "contig": contig,
                "start": start,
                "stop": stop,
                "length_bp": length,
                "gc_percent": round(gc_fraction(seq) * 100, 2),
                "genome_gc_percent": round(genome_gc * 100, 2),
                "gc_delta": round((gc_fraction(seq) - genome_gc) * 100, 2),
                "cds_count": len(cds_list),
                "cds_per_kb": round(len(cds_list) / (length / 1000), 2) if length else 0,
                "confidence": level,
                "complete_like": int(complete_like),
                "hallmark_categories": ";".join(sorted(cats)),
                "terminase_n": counts["terminase"],
                "capsid_n": counts["capsid"],
                "tail_n": counts["tail"],
                "lysin_n": counts["lysin"],
                "integrase_n": counts["integrase"],
                "regulation_n": counts["regulation"],
                "replication_n": counts["replication"],
                "attL": pp["attL"],
                "attR": pp["attR"],
                "best_phage": "",
                "best_pident": "",
                "best_qcov": "",
                "phage_typeable": 0,
            })
        if coords:
            qc["genomes_with_prophage"] += 1
        genome_rows.append({
            "assembly_accession": acc, "strain": rep.get("strain", ""),
            "ST": rep.get("ST", ""), "mlst_clade": rep.get("mlst_clade", ""),
            "genome_size": genome_len, "genome_gc_percent": round(genome_gc * 100, 2),
            "n_prophages": len([r for r in rows if r["assembly_accession"] == acc]),
            "prophage_bp": pp_bp,
            "prophage_percent": round(pp_bp / genome_len * 100, 3) if genome_len else 0,
        })
    handle.close()

    # BLAST all prophage sequences against the complete phage reference set.
    blast_tsv = os.path.join(RESULT_DIR, "prophage_vs_phage_refs.tsv")
    cmd = [BLASTN, "-query", fasta_out, "-db", PHAGE_DB,
           "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
           "-max_target_seqs", "5", "-evalue", "1e-10", "-num_threads", "8",
           "-out", blast_tsv]
    subprocess.run(cmd, check=True)
    best = {}
    with open(blast_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            qid, sid, pident, alen, qlen = parts[0], parts[1], float(parts[2]), int(parts[3]), int(parts[4])
            qcov = alen / qlen if qlen else 0
            score = pident * qcov
            if qid not in best or score > best[qid][0]:
                best[qid] = (score, sid, pident, qcov)
    for row in rows:
        hit = best.get(row["prophage_id"])
        if hit:
            _, sid, pident, qcov = hit
            row["best_phage"] = sid
            row["best_pident"] = round(pident, 2)
            row["best_qcov"] = round(qcov, 3)
            row["phage_typeable"] = int(pident >= 70 and qcov >= 0.5)

    with open(os.path.join(RESULT_DIR, "prophage_details_raw.csv"), "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else [])
        writer.writeheader()
        writer.writerows(rows)

    with open(os.path.join(RESULT_DIR, "genome_prophage_summary_raw.csv"), "w", newline="") as fh:
        fields = ["assembly_accession", "strain", "ST", "mlst_clade", "genome_size",
                  "genome_gc_percent", "n_prophages", "prophage_bp", "prophage_percent"]
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(genome_rows)

    qc["prophages_total"] = len(rows)
    qc["confidence_counts"] = dict(Counter(r["confidence"] for r in rows))
    qc["complete_like"] = sum(r["complete_like"] for r in rows)
    qc["typeable"] = sum(r["phage_typeable"] for r in rows)
    qc["genomes_without_prophage"] = [g["assembly_accession"] for g in genome_rows if g["n_prophages"] == 0]
    with open(os.path.join(RESULT_DIR, "prophage_annotation_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=2)
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
