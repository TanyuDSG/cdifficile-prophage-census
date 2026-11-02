#!/usr/bin/env python3
"""Annotate prophage CDS and cargo categories for the 241-genome dataset.

Keyword annotation is transparent and reproducible from the deposited GenBank
records. The chromosome baseline explicitly excludes every CDS overlapping a
PhiSpy prophage region, removing the circularity of the pilot analysis.

Outputs (scaleup_241/results/)
  prophage_cds_annotation.csv
  prophage_proteins.fasta
  prophage_cargo_summary.csv
  genome_cargo_counts.csv
  cargo_summary.json
"""
import json
import os
import re
from collections import Counter, defaultdict

import pandas as pd
from Bio import SeqIO
from scipy import stats

BASE = "/home/bioinfo/cdprophage"
GBFF = os.path.join(BASE, "genomes", "gbff")
RESULT = os.path.join(BASE, "scaleup_241", "results")

DEFENCE = {
    "CRISPR-Cas": re.compile(
        r"\bcas[0-9]+[a-z]?\b|\bcrispr\b|\bcascade\b|\bcsn2\b|\bcsm[0-9]\b|"
        r"\bcmr[0-9]\b|\bcas10\b|\bcas9\b|\bcsx[0-9]+\b", re.I),
    "Restriction-Modification": re.compile(
        r"restriction endonuclease|restriction-modification|restriction modification|"
        r"\bhsd[msr]\b|\becoki\b|\bmrr\b|type i restriction|type ii restriction|"
        r"type iii restriction|\bdna methyltransferase\b|\bdam methylase\b|\bdcm\b", re.I),
    "Toxin-Antitoxin": re.compile(
        r"toxin-antitoxin|toxin antitoxin|\bantitoxin\b|\brel[eE]\b|\bmaz[ef]\b|"
        r"\bpem[ik]\b|\bvap[bc]\b|\bhig[ab]\b|\byaf[oq]\b|"
        r"\bhic[ab]\b|\bmqs[ar]\b|\brnl[ab]\b|\bsym[er]\b|\bdoc\b|\bfic\b", re.I),
    "Abortive infection": re.compile(
        r"abortive infection|\babi[abcdefghijklmnopqrstuvwxyz]\b|\brex[ab]\b|"
        r"\brlo[c]\b|\bprr[c]\b|\blit\b|\bpif\b", re.I),
    "Argonaute": re.compile(r"argonaute|\bago\b", re.I),
    "CBASS/cyclic-oligo": re.compile(
        r"\bcbass\b|\bcgas\b|\bcd-ntase\b|\bcdntase\b|\bcap[0-9]\b|"
        r"cyclic oligonucleotide|phosphodiesterase", re.I),
    "Other defence": re.compile(
        r"antiphage|anti-phage|bacteriophage resistance|phage resistance|"
        r"\bthoeris\b|\btir\b|\bsirt\b|gasdermin|\bnin[bc]\b|\bavcad\b|"
        r"\bdsra\b|\bdsrb\b|phage shock", re.I),
}
TA_TOXIN = re.compile(
    r"toxin-antitoxin system toxin|toxin antitoxin system toxin|"
    r"\brel[eE]\b|\bmaz[ef]\b|\bpem[ik]\b|\bvap[bc]\b|\bhigB\b|"
    r"\byaf[oq]\b|\bhicA\b|\bmqsR\b|\brnlA\b|\bsymE\b|"
    r"\bdoc\b|\bfic\b", re.I)
TA_ANTITOXIN = re.compile(
    r"toxin-antitoxin system antitoxin|toxin antitoxin system antitoxin|"
    r"family antitoxin|system antitoxin|^antitoxin$|"
    r"\brel[bB]\b|\bmaz[eE]\b|\bpem[iI]\b|\bvap[bB]\b|"
    r"\bhigA\b|\byaf[nN]\b|\bhicB\b|\bmqsA\b|\brnlB\b|"
    r"\bsymR\b", re.I)
RM_ENDO = re.compile(
    r"restriction endonuclease|\bhsd[rR]\b|\bmrr\b|type i restriction|"
    r"type ii restriction|type iii restriction", re.I)
RM_METH = re.compile(
    r"dna methyltransferase|\bhsd[mM]\b|\bdam methylase\b|\bdcm\b", re.I)
CRISPR_CAS = re.compile(
    r"\bcas[0-9]+[a-z]?\b|\bcrispr-associated\b|\bcascade\b|\bcsn2\b|"
    r"\bcsm[0-9]\b|\bcmr[0-9]\b|\bcas10\b|\bcas9\b|\bcsx[0-9]+\b", re.I)

VIRULENCE = {
    "TcdA/TcdB toxin": re.compile(
        r"\btcd[ab]\b|\btoxin a\b|\btoxin b\b|\btoxin-a\b|\btoxin-b\b|"
        r"cytotoxin a|cytotoxin b", re.I),
    "Binary toxin CDT": re.compile(
        r"\bcdt[ab]\b|binary toxin|clostridial binary toxin", re.I),
    "Adhesion/S-layer": re.compile(
        r"adhesin|collagen-binding|collagen binding|fibronectin-binding|"
        r"fibronectin binding|lpxtg|sortase|s-layer|surface layer|"
        r"cell wall anchor|internalin|fimbri|pilin|flagell", re.I),
    "Capsule/spore/germination": re.compile(
        r"capsul|spore coat|spore protein|germination|exosporium|"
        r"sporulation", re.I),
    "Secreted enzyme/hemolysin": re.compile(
        r"hemolysin|haemolysin|phospholipase|collagenase|hyaluronidase|"
        r"extracellular protease|serine protease|cysteine protease|"
        r"metalloprotease|secreted protease", re.I),
    "Immune evasion": re.compile(
        r"immune evasion|phagocytosis|complement|superantigen|"
        r"immunoglobulin-binding|ige-binding", re.I),
}
AMR = {
    "Beta-lactam": re.compile(
        r"beta-lactamase|β-lactamase|\bbla\b|\boxa\b|\bndm\b|\bkpc\b|\bmcr\b", re.I),
    "Aminoglycoside": re.compile(
        r"aminoglycoside|\baac\(|\baph\(|\bant\(|\bstr[ab]\b", re.I),
    "Tetracycline": re.compile(r"tetracycline|\btet[amqosw]\b", re.I),
    "Macrolide/Lincosamide": re.compile(
        r"macrolide|\berm[abc]?\b|\bmph[ab]?\b|\bmef[ae]?\b|\blnu[ab]?\b|"
        r"\bmsr[ab]?\b", re.I),
    "Chloramphenicol": re.compile(r"chloramphenicol|\bcat[a]?\b", re.I),
    "Glycopeptide": re.compile(r"vancomycin|\bvan[ab]\b", re.I),
    "Sulfonamide/Trimethoprim": re.compile(
        r"sulfonamide|\bsul[12]\b|trimethoprim|\bdfr[ab]?\b", re.I),
    "Fluoroquinolone": re.compile(r"fluoroquinolone|\bqnr\b|quinolone resistance", re.I),
    "Multidrug efflux": re.compile(
        r"multidrug|multi-drug|efflux pump|efflux transporter|acrab|mex[ab]", re.I),
}
MGE = {
    "Transposase/IS": re.compile(
        r"transposase|insertion sequence|\bis[0-9]+\b|\btn[0-9]+\b|"
        r"transposon", re.I),
    "Conjugation/Mobilization": re.compile(
        r"conjugal transfer|conjugation|relaxase|mobilization|mobilase|"
        r"\btrae?\b|\btrad\b|\btrb[abcdefg]\b|\bvirb[0-9]+\b", re.I),
    "Plasmid replication": re.compile(
        r"plasmid replication|plasmid replicon|\brep[ab]\b|rolling-circle", re.I),
}
PHAGE_CORE = re.compile(
    r"terminase|capsid|portal protein|major head|prohead|\btail\b|tail fiber|"
    r"tail fibre|tape measure|baseplate|tail sheath|tail tube|neck protein|"
    r"lysin|endolysin|holin|amidase|muramidase|spanin", re.I)


def text_of(feat):
    q = feat.qualifiers
    return " ".join(" ".join(q.get(k, [])) for k in ("product", "gene", "note", "locus_tag"))


def category_hits(text, patterns):
    return [name for name, pat in patterns.items() if pat.search(text)]


def main():
    details = pd.read_csv(os.path.join(RESULT, "prophage_details.csv"),
                          dtype={"ST": str, "mlst_clade": str})
    regions = defaultdict(list)
    for _, r in details.iterrows():
        regions[r["assembly_accession"]].append(r)

    cds_rows = []
    genome_rows = []
    handle = open(os.path.join(RESULT, "prophage_proteins.fasta"), "w")

    for acc, sub in regions.items():
        gbff = os.path.join(GBFF, acc + ".gbff")
        genome_cds = []
        chrom_cds = []
        genome_crispr = 0
        chrom_crispr = 0
        for rec in SeqIO.parse(gbff, "genbank"):
            contig = rec.id
            for feat in rec.features:
                if feat.type == "repeat_region":
                    txt = text_of(feat)
                    if re.search(r"crispr", txt, re.I):
                        genome_crispr += 1
                        overlap = [r for r in sub if r["contig"] == contig
                                   and feat.location.end > r["start"] and feat.location.start < r["stop"]]
                        if not overlap:
                            chrom_crispr += 1
                    continue
                if feat.type != "CDS":
                    continue
                txt = text_of(feat)
                prot = " ".join(feat.qualifiers.get("translation", []))
                start, end = int(feat.location.start), int(feat.location.end)
                genome_cds.append((feat, txt, prot))
                assigned = None
                for r in sub:
                    if r["contig"] != contig:
                        continue
                    if end <= r["start"] or start >= r["stop"]:
                        continue
                    assigned = r
                    break
                if assigned is None:
                    chrom_cds.append((feat, txt, prot))
                    continue
                r = assigned
                defence = category_hits(txt, DEFENCE)
                virulence = category_hits(txt, VIRULENCE)
                amr = category_hits(txt, AMR)
                mge = category_hits(txt, MGE)
                core = bool(PHAGE_CORE.search(txt))
                locus = " ".join(feat.qualifiers.get("locus_tag", []))
                protein_key = f"{r['prophage_id']}|{locus if locus else start}"
                row = {
                    "prophage_id": r["prophage_id"],
                    "assembly_accession": acc,
                    "ST": r["ST"],
                    "mlst_clade": r["mlst_clade"],
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": feat.location.strand,
                    "locus_tag": locus,
                    "protein_key": protein_key,
                    "gene": " ".join(feat.qualifiers.get("gene", [])),
                    "product": " ".join(feat.qualifiers.get("product", [])),
                    "protein_length": len(prot.rstrip("*")),
                    "phage_core": int(core),
                    "defence_categories": ";".join(defence),
                    "virulence_categories": ";".join(virulence),
                    "amr_categories": ";".join(amr),
                    "mge_categories": ";".join(mge),
                    "is_defence": int(bool(defence)),
                    "is_virulence": int(bool(virulence)),
                    "is_amr": int(bool(amr)),
                    "is_mge": int(bool(mge)),
                    "ta_toxin": int(bool(TA_TOXIN.search(txt))),
                    "ta_antitoxin": int(bool(TA_ANTITOXIN.search(txt))),
                    "rm_endonuclease": int(bool(RM_ENDO.search(txt))),
                    "rm_methyltransferase": int(bool(RM_METH.search(txt))),
                    "crispr_cas": int(bool(CRISPR_CAS.search(txt))),
                }
                cds_rows.append(row)
                if prot and not prot.startswith("MISSING"):
                    handle.write(f">{protein_key}\n{prot}\n")

        def count(cds_list, patterns):
            return sum(1 for _, t, _ in cds_list if category_hits(t, patterns))

        genome_rows.append({
            "assembly_accession": acc,
            "genome_cds": len(genome_cds),
            "chromosome_cds": len(chrom_cds),
            "genome_defence_cds": count(genome_cds, DEFENCE),
            "genome_virulence_cds": count(genome_cds, VIRULENCE),
            "genome_amr_cds": count(genome_cds, AMR),
            "genome_mge_cds": count(genome_cds, MGE),
            "chromosome_defence_cds": count(chrom_cds, DEFENCE),
            "chromosome_virulence_cds": count(chrom_cds, VIRULENCE),
            "chromosome_amr_cds": count(chrom_cds, AMR),
            "chromosome_mge_cds": count(chrom_cds, MGE),
            "genome_crispr_arrays": genome_crispr,
            "chromosome_crispr_arrays": chrom_crispr,
        })
    handle.close()

    cds = pd.DataFrame(cds_rows)
    cds.to_csv(os.path.join(RESULT, "prophage_cds_annotation.csv"), index=False)

    per = []
    for pid, sub in cds.groupby("prophage_id"):
        row = details[details["prophage_id"] == pid].iloc[0]
        per.append({
            "prophage_id": pid,
            "assembly_accession": row["assembly_accession"],
            "ST": row["ST"],
            "mlst_clade": row["mlst_clade"],
            "cluster_id": row["cluster_id"],
            "length_bp": row["length_bp"],
            "gc_percent": row["gc_percent"],
            "gc_delta": row["gc_delta"],
            "confidence": row["confidence"],
            "complete_like": row["complete_like"],
            "n_cds": len(sub),
            "n_defence": int(sub["is_defence"].sum()),
            "n_virulence": int(sub["is_virulence"].sum()),
            "n_amr": int(sub["is_amr"].sum()),
            "n_mge": int(sub["is_mge"].sum()),
            "has_defence": int(sub["is_defence"].any()),
            "has_ta_system": int(sub["ta_toxin"].any() and sub["ta_antitoxin"].any()),
            "has_rm_system": int(sub["rm_endonuclease"].any() and sub["rm_methyltransferase"].any()),
            "has_crispr_cas": int(sub["crispr_cas"].any()),
            "has_virulence": int(sub["is_virulence"].any()),
            "has_amr": int(sub["is_amr"].any()),
            "has_mge": int(sub["is_mge"].any()),
            "defence_categories": ";".join(sorted(set(
                c for s in sub["defence_categories"] for c in s.split(";") if c))),
            "virulence_categories": ";".join(sorted(set(
                c for s in sub["virulence_categories"] for c in s.split(";") if c))),
            "amr_categories": ";".join(sorted(set(
                c for s in sub["amr_categories"] for c in s.split(";") if c))),
            "mge_categories": ";".join(sorted(set(
                c for s in sub["mge_categories"] for c in s.split(";") if c))),
        })
    per = pd.DataFrame(per)
    per["n_crispr_arrays"] = 0
    per["has_crispr_array"] = 0
    per["has_crispr_system"] = 0
    per["has_defence_system"] = (
        (per["has_ta_system"] == 1) | (per["has_rm_system"] == 1) | (per["has_crispr_cas"] == 1)
    ).astype(int)
    per.to_csv(os.path.join(RESULT, "prophage_cargo_summary.csv"), index=False)

    genome_df = pd.DataFrame(genome_rows)
    pp_counts = cds.groupby("assembly_accession").agg(
        prophage_cds=("prophage_id", "size"),
        prophage_defence_cds=("is_defence", "sum"),
        prophage_virulence_cds=("is_virulence", "sum"),
        prophage_amr_cds=("is_amr", "sum"),
        prophage_mge_cds=("is_mge", "sum"),
    ).reset_index()
    genome_df = genome_df.merge(pp_counts, on="assembly_accession", how="left").fillna(0)
    genome_df.to_csv(os.path.join(RESULT, "genome_cargo_counts.csv"), index=False)

    high = per[per["confidence"] == "high"]
    summary = {
        "prophages": int(len(per)),
        "high_confidence_prophages": int(len(high)),
        "cds_total": int(len(cds)),
        "defence_cds": int(cds["is_defence"].sum()),
        "virulence_cds": int(cds["is_virulence"].sum()),
        "amr_cds": int(cds["is_amr"].sum()),
        "mge_cds": int(cds["is_mge"].sum()),
        "prophages_with_defence": int(per["has_defence"].sum()),
        "prophages_with_virulence": int(per["has_virulence"].sum()),
        "prophages_with_amr": int(per["has_amr"].sum()),
        "prophages_with_mge": int(per["has_mge"].sum()),
        "defence_rate_all": round(float(per["has_defence"].mean()), 4),
        "virulence_rate_all": round(float(per["has_virulence"].mean()), 4),
        "amr_rate_all": round(float(per["has_amr"].mean()), 4),
        "mge_rate_all": round(float(per["has_mge"].mean()), 4),
        "defence_rate_high": round(float(high["has_defence"].mean()), 4),
        "virulence_rate_high": round(float(high["has_virulence"].mean()), 4),
        "amr_rate_high": round(float(high["has_amr"].mean()), 4),
        "mge_rate_high": round(float(high["has_mge"].mean()), 4),
        "defence_category_region_counts": dict(Counter(
            c for s in per["defence_categories"] for c in s.split(";") if c)),
        "defence_category_counts": dict(Counter(
            c for s in per["defence_categories"] for c in s.split(";") if c)),
        "category_counts_unit": "regions; categories are non-exclusive",
        "virulence_category_counts": dict(Counter(
            c for s in per["virulence_categories"] for c in s.split(";") if c)),
        "amr_category_counts": dict(Counter(
            c for s in per["amr_categories"] for c in s.split(";") if c)),
        "mge_category_counts": dict(Counter(
            c for s in per["mge_categories"] for c in s.split(";") if c)),
        "genome_cds_mean": round(float(genome_df["genome_cds"].mean()), 2),
        "chromosome_cds_mean": round(float(genome_df["chromosome_cds"].mean()), 2),
        "genome_defence_cds_mean": round(float(genome_df["genome_defence_cds"].mean()), 2),
        "chromosome_defence_cds_mean": round(float(genome_df["chromosome_defence_cds"].mean()), 2),
        "genome_virulence_cds_mean": round(float(genome_df["genome_virulence_cds"].mean()), 2),
        "chromosome_virulence_cds_mean": round(float(genome_df["chromosome_virulence_cds"].mean()), 2),
        "genome_amr_cds_mean": round(float(genome_df["genome_amr_cds"].mean()), 2),
        "chromosome_amr_cds_mean": round(float(genome_df["chromosome_amr_cds"].mean()), 2),
        "genome_mge_cds_mean": round(float(genome_df["genome_mge_cds"].mean()), 2),
        "chromosome_mge_cds_mean": round(float(genome_df["chromosome_mge_cds"].mean()), 2),
    }
    enrich = {}
    for cat in ("defence", "virulence", "amr", "mge"):
        pp = genome_df["prophage_%s_cds" % cat] / genome_df["prophage_cds"].replace(0, float("nan"))
        bg = genome_df["chromosome_%s_cds" % cat] / genome_df["chromosome_cds"]
        ok = pp.notna()
        try:
            stat, pval = stats.wilcoxon(pp[ok], bg[ok], alternative="two-sided", zero_method="wilcox")
        except ValueError:
            stat, pval = float("nan"), 1.0
        enrich[cat] = {
            "prophage_cds_fraction_mean": round(float(pp[ok].mean()), 5),
            "chromosome_cds_fraction_mean": round(float(bg[ok].mean()), 5),
            "median_difference": round(float((pp[ok] - bg[ok]).median()), 5),
            "wilcoxon_stat": round(float(stat), 2) if stat == stat else None,
            "p_value": float(pval),
            "n_genomes": int(ok.sum()),
        }
    ordered = sorted(enrich, key=lambda c: enrich[c]["p_value"])
    n_tests = len(ordered)
    prev_q = 1.0
    for rank, cat in reversed(list(enumerate(ordered, 1))):
        q = min(prev_q, enrich[cat]["p_value"] * n_tests / rank)
        enrich[cat]["q_value_bh"] = float(q)
        prev_q = q
    summary["enrichment_vs_chromosome"] = enrich
    with open(os.path.join(RESULT, "cargo_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
