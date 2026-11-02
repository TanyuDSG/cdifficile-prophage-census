#!/usr/bin/env python3
"""Assemble Supplementary Tables S1-S8 into a single workbook."""
import importlib.util
import json
import os

import pandas as pd
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter

BASE = "/home/bioinfo/cdprophage"
RES = os.path.join(BASE, "scaleup_241", "results")
OUT = os.path.join(BASE, "scaleup_241", "Supplementary_Tables_S1-S8.xlsx")


def read(name, **kw):
    path = os.path.join(RES, name)
    if not os.path.exists(path):
        return None
    return pd.read_csv(path, **kw)


def load_json(name):
    path = os.path.join(RES, name)
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


# S1 genomes
g = read("genome_summary_analysis.csv", dtype={"ST": str, "mlst_clade": str})
meta = read("genome_metadata_full.csv", dtype={"ST": str, "mlst_clade": str})
if g is not None:
    cols = ["assembly_accession", "strain", "ST", "mlst_clade", "genome_size",
            "genome_gc_percent", "country", "collection_year", "host", "isolation_source",
            "n_prophages", "prophage_bp", "prophage_percent", "n_high",
            "n_complete_like", "n_typeable", "n_clusters"]
    if meta is not None:
        keep = [c for c in cols if c in meta.columns]
        s1 = meta[keep].copy()
        extra = g[["assembly_accession", "n_prophages", "prophage_bp", "prophage_percent",
                   "n_high", "n_complete_like", "n_typeable", "n_clusters"]]
        s1 = s1.drop(columns=[c for c in extra.columns if c != "assembly_accession" and c in s1.columns])
        s1 = s1.merge(extra, on="assembly_accession", how="left")
    else:
        s1 = g[[c for c in cols if c in g.columns]].copy()
    busco = read("busco_summary.csv")
    if busco is not None:
        s1 = s1.merge(busco[["assembly_accession", "complete_pct", "missing_pct"]],
                      on="assembly_accession", how="left")
    s1 = s1.rename(columns={
        "assembly_accession": "Assembly accession", "strain": "Strain", "ST": "ST",
        "mlst_clade": "MLST clade", "genome_size": "Genome size (bp)",
        "genome_gc_percent": "GC (%)", "country": "Country",
        "collection_year": "Collection year", "host": "Host",
        "isolation_source": "Isolation source", "n_prophages": "Prophage regions",
        "prophage_bp": "Prophage DNA (bp)", "prophage_percent": "Prophage DNA (%)",
        "n_high": "High-confidence regions", "n_complete_like": "Complete-like regions",
        "n_typeable": "Reference-typeable regions", "n_clusters": "Distinct clusters",
        "complete_pct": "BUSCO complete (%)", "missing_pct": "BUSCO missing (%)"})
    s1 = s1.sort_values(["MLST clade", "ST", "Strain"], key=lambda c: c.astype(str))

# S2 clusters
s2 = read("cluster_summary.csv")
if s2 is not None:
    details = read("prophage_details.csv", dtype={"ST": str, "mlst_clade": str})
    if details is not None:
        rep = details.sort_values("length_bp", ascending=False).drop_duplicates("cluster_id")
        rep = rep.set_index("cluster_id")
        s2["Representative"] = s2["cluster_id"].map(rep["prophage_id"])
        s2["Representative accession"] = s2["cluster_id"].map(rep["assembly_accession"])
        s2["Length (bp)"] = s2["cluster_id"].map(rep["length_bp"])
        s2["GC (%)"] = s2["cluster_id"].map(rep["gc_percent"])
        s2["Confidence"] = s2["cluster_id"].map(rep["confidence"])
        s2["Closest reference phage"] = s2["cluster_id"].map(rep["best_phage"]).fillna("")
        conf_counts = details.groupby("cluster_id")["confidence"].apply(
            lambda s: ";".join(f"{k}:{v}" for k, v in s.value_counts().items()))
        s2["Confidence mix"] = s2["cluster_id"].map(conf_counts)
    s2 = s2.rename(columns={
        "cluster_id": "Cluster", "n_prophages": "Regions", "n_genomes": "Genomes",
        "n_clades": "MLST clades", "n_STs": "STs",
        "median_length_bp": "Median length (bp)", "min_length_bp": "Min length (bp)",
        "max_length_bp": "Max length (bp)", "mean_gc_percent": "Mean GC (%)",
        "n_high": "High-confidence regions", "n_complete_like": "Complete-like regions",
        "n_typeable": "Typeable regions"})

# S2b species-level clusters at 95% identity (robustness check)
s2b = read("species_clusters_95.csv")
if s2b is not None:
    s2b = s2b.rename(columns={
        "prophage_id": "Region", "species_cluster_95": "Species-level cluster (95%)",
        "cluster_size": "Regions in cluster", "representative": "Representative region"})
    s2b = s2b.sort_values(["Species-level cluster (95%)", "Region"])

# S3 protein families
s3 = read("protein_gene_families.csv")

# S4a cargo per prophage
s4a = read("prophage_cargo_summary.csv", dtype={"ST": str, "mlst_clade": str})
if s4a is not None:
    s4a = s4a.rename(columns={
        "prophage_id": "Prophage", "assembly_accession": "Assembly accession",
        "strain": "Strain", "ST": "ST", "mlst_clade": "MLST clade",
        "length_bp": "Length (bp)", "confidence": "Confidence",
        "cds_count": "CDS", "n_defence": "Defence CDS", "n_virulence": "Virulence CDS",
        "n_amr": "AMR CDS", "n_mge": "MGE CDS", "n_df_systems": "Complete defence systems",
        "n_df_partial": "Partial defence systems", "df_system_types": "Defence system types"})

# S4b keyword dictionary
spec = importlib.util.spec_from_file_location(
    "cargo", os.path.join(BASE, "scaleup_241", "scripts", "annotate_cargo_all.py"))
cargo_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cargo_mod)
rows = []
for cat, d in (("Defence", cargo_mod.DEFENCE), ("Virulence", cargo_mod.VIRULENCE),
               ("AMR", cargo_mod.AMR), ("MGE", cargo_mod.MGE)):
    for sub, pattern in d.items():
        rows.append({"Category": cat, "Subcategory": sub, "Regular expression": pattern.pattern})
rows.append({"Category": "Phage core", "Subcategory": "phage_core",
             "Regular expression": cargo_mod.PHAGE_CORE.pattern})
s4b = pd.DataFrame(rows)

# S5a tests
s5a = read("statistical_tests.csv")
mp = load_json("mantel_partial.json")
extra_tests = []
if "mantel_year" in mp:
    extra_tests.append({"test": "Mantel", "metric": "prophage cluster Jaccard vs collection year distance",
                        "groups": f"genomes with metadata (n={mp['mantel_year']['n']})",
                        "statistic": mp["mantel_year"]["r"], "p_value": mp["mantel_year"]["p"],
                        "q_value_BH": "", "note": "999 permutations"})
    extra_tests.append({"test": "Mantel", "metric": "prophage cluster Jaccard vs country distance",
                        "groups": f"genomes with metadata (n={mp['mantel_country']['n']})",
                        "statistic": mp["mantel_country"]["r"], "p_value": mp["mantel_country"]["p"],
                        "q_value_BH": "", "note": "999 permutations"})
    extra_tests.append({"test": "Partial Mantel", "metric": "prophage cluster Jaccard vs Mash distance | collection year",
                        "groups": f"genomes with metadata (n={mp['partial_mantel_year']['n']})",
                        "statistic": mp["partial_mantel_year"]["r"], "p_value": mp["partial_mantel_year"]["p"],
                        "q_value_BH": "", "note": "999 permutations; p=0.001 is resolution floor"})
    if "partial_mantel_year_country" in mp:
        extra_tests.append({"test": "Partial Mantel", "metric": "prophage cluster Jaccard vs Mash distance | year + country",
                            "groups": f"genomes with metadata (n={mp['partial_mantel_year_country']['n']})",
                            "statistic": mp["partial_mantel_year_country"]["r"],
                            "p_value": mp["partial_mantel_year_country"]["p"],
                            "q_value_BH": "", "note": "999 permutations; p=0.001 is resolution floor"})
if s5a is not None and extra_tests:
    s5a = pd.concat([s5a, pd.DataFrame(extra_tests)], ignore_index=True)

# S5b paired cargo comparisons
cs = load_json("cargo_summary.json")
rows = []
for cat, d in cs.get("enrichment_vs_chromosome", {}).items():
    rows.append({"Category": cat, "Prophage CDS fraction": d["prophage_cds_fraction_mean"],
                 "Chromosome CDS fraction": d["chromosome_cds_fraction_mean"],
                 "Median paired difference": d["median_difference"],
                 "Wilcoxon statistic": d["wilcoxon_stat"], "P value": d["p_value"],
                 "q value (BH)": d["q_value_bh"], "Genomes": d["n_genomes"]})
s5b = pd.DataFrame(rows)

# S6a sensitivity
sens = cs and load_json("diversity_summary.json").get("sensitivity_without_630Derm", {})
ds = load_json("diversity_summary.json")
sens_rows = [{"Analysis": "Primary (all complete genomes)", "Genomes": ds.get("genomes"),
              "Prophage regions": ds.get("prophages"), "Clusters": ds.get("clusters"),
              "Singleton clusters": ds.get("singleton_clusters"),
              "Carriage rate": ds.get("carriage_rate"),
              "Mean regions per genome": ds.get("prophages_per_genome", {}).get("mean"),
              "Prophage DNA (%)": ds.get("prophage_percent", {}).get("mean")}]
if sens:
    sens_rows.append({"Analysis": "Excluding 630Derm (lab derivative)", "Genomes": sens.get("genomes"),
                      "Prophage regions": sens.get("prophages"), "Clusters": sens.get("clusters"),
                      "Singleton clusters": sens.get("singleton_clusters"),
                      "Carriage rate": sens.get("carriage_rate"),
                      "Mean regions per genome": sens.get("prophages_per_genome"),
                      "Prophage DNA (%)": sens.get("prophage_percent")})
if details is not None:
    hi = details[details["confidence"] == "high"]
    sens_rows.append({"Analysis": "High-confidence regions only",
                      "Genomes": hi["assembly_accession"].nunique(),
                      "Prophage regions": len(hi), "Clusters": hi["cluster_id"].nunique(),
                      "Singleton clusters": int((hi.groupby("cluster_id")["assembly_accession"].nunique() == 1).sum()),
                      "Carriage rate": "", "Mean regions per genome": round(len(hi) / max(hi["assembly_accession"].nunique(), 1), 2),
                      "Prophage DNA (%)": ""})
s6a = pd.DataFrame(sens_rows)

# S6b QC
qc = load_json("prophage_annotation_qc.json")
s6b = pd.DataFrame([{"Metric": k, "Value": json.dumps(v) if isinstance(v, (list, dict)) else v}
                    for k, v in qc.items()])

# S7 reference phages
s7 = read("reference_phage_types.csv")

# S8 excluded candidate regions
allreg = read("prophage_details_all.csv", dtype={"ST": str, "mlst_clade": str})
s8 = None
if allreg is not None:
    ex = allreg[allreg["in_primary"] == 0].copy()

    def reason(r):
        if int(r["terminase_n"]) > 0 or int(r["capsid_n"]) > 0:
            return "packaging module (retained)"
        if int(r["tail_n"]) > 0:
            return "tail-only contractile element"
        if int(r["lysin_n"]) > 0:
            return "lysin only"
        return "no structural hallmark"

    ex["Exclusion reason"] = ex.apply(reason, axis=1)
    s8 = ex[["prophage_id", "assembly_accession", "ST", "mlst_clade", "contig",
             "start", "stop", "length_bp", "confidence", "terminase_n", "capsid_n",
             "tail_n", "lysin_n", "Exclusion reason"]].rename(columns={
        "prophage_id": "Region", "assembly_accession": "Assembly accession",
        "ST": "ST", "mlst_clade": "MLST clade", "length_bp": "Length (bp)",
        "confidence": "Confidence", "terminase_n": "Terminase genes",
        "capsid_n": "Capsid genes", "tail_n": "Tail genes", "lysin_n": "Lysin genes"})
    s8 = s8.sort_values(["Exclusion reason", "Length (bp)"], ascending=[True, False])

sheets = [
    ("S1_Genomes", s1), ("S2_Clusters", s2), ("S2b_Species95", s2b),
    ("S3_ProteinFamilies", s3),
    ("S4a_Cargo", s4a), ("S4b_Dictionary", s4b), ("S5a_Tests", s5a),
    ("S5b_PairedCargo", s5b), ("S6a_Sensitivity", s6a), ("S6b_QC", s6b),
    ("S7_ReferencePhages", s7), ("S8_ExcludedRegions", s8),
]
readme = pd.DataFrame({
    "Sheet": [s[0] for s in sheets],
    "Content": [
        "Genome panel, metadata, BUSCO v5.4.7 completeness and prophage carriage (n = 241)",
        "Prophage nucleotide clusters (70% identity, 50% coverage, single linkage)",
        "Species-level clusters at 95% identity over 50% of both sequences (robustness check)",
        "Protein gene families (CD-HIT 50% identity, 50% coverage) and clade distribution",
        "Per-prophage cargo assignments and DefenseFinder systems",
        "Keyword dictionary used for cargo annotation",
        "All statistical tests, test statistics and adjusted P values",
        "Paired prophage-versus-chromosome cargo comparisons",
        "Sensitivity analyses (630Derm excluded; high-confidence subset)",
        "PhiSpy annotation quality-control summary",
        "Complete reference phage and phage-like sequences used for typeability",
        "Candidate regions excluded from the primary catalogue (no packaging module)",
    ],
})

with pd.ExcelWriter(OUT, engine="openpyxl") as writer:
    readme.to_excel(writer, sheet_name="README", index=False)
    for name, df in sheets:
        if df is None:
            pd.DataFrame({"Note": ["Not generated"]}).to_excel(writer, sheet_name=name, index=False)
        else:
            df.to_excel(writer, sheet_name=name, index=False)

# formatting
import openpyxl
wb = openpyxl.load_workbook(OUT)
for ws in wb.worksheets:
    ws.freeze_panes = "A2"
    for cell in ws[1]:
        cell.font = Font(bold=True, color="FFFFFF")
        cell.fill = PatternFill("solid", fgColor="1F4E79")
        cell.alignment = Alignment(vertical="center", wrap_text=True)
    for col in ws.columns:
        width = max(len(str(c.value)) if c.value is not None else 0 for c in col)
        ws.column_dimensions[get_column_letter(col[0].column)].width = min(max(10, width + 2), 60)
wb.save(OUT)
print("written", OUT, os.path.getsize(OUT))
