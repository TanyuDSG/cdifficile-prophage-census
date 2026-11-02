# *Clostridioides difficile* prophage census

Data and code accompanying the manuscript **"Species-wide prophage diversity and
defence-cargo enrichment in *Clostridioides difficile*: a census of 241 complete
genomes"**, submitted to *Microbial Genomics*.

| Quantity | Value |
|---|---|
| Complete RefSeq genomes analysed | 241 |
| Predicted prophage regions (all candidates) | 1,380 |
| Regions in the primary catalogue | 399 |
| Nucleotide clusters | 130 |
| Singleton clusters | 76 (58.5%) |
| Protein gene families | 1,447 |
| DefenceFinder system calls | 18 (16 complete, 2 incomplete) |

The repository holds the derived tables that reproduce every number, table and
figure of the manuscript, the prophage nucleotide sequences, the figure source
files and the analysis pipeline. It does **not** contain the manuscript or the
genome assemblies themselves: all 241 assemblies are public RefSeq records and
are listed by accession in Table S1 of the manuscript.

## Contents

```
data/       derived tables (14 files, ~19 MB)
scripts/    analysis and figure code (28 files)
figures/    Figs 1-6 and S1-S2, 600 dpi PNG + vector PDF
README.md   this file
requirements.txt
LICENSE         MIT (code)
LICENSE-DATA    CC BY 4.0 (data)
```

## Data dictionary

| File | Rows | Description |
|---|---|---|
| `genome_metadata_full.csv` | 241 | One row per genome: accession, strain, ST, MLST clade, genome size, GC, release date, BioSample/BioProject, country, year, host, isolation source |
| `genome_summary_analysis.csv` | 241 | Per-genome prophage summary: region count, prophage bp and percentage, cluster count, confidence tiers, within-genome Shannon index, carriage flags |
| `busco_summary.csv` | 241 | BUSCO completeness per genome (clostridia_odb10, 247 markers) |
| `prophage_details.csv` | 399 | Primary catalogue: coordinates, length, GC and GC delta, CDS count, hallmark gene counts, attL/attR, best BLASTn hit, confidence tier, cluster assignment |
| `prophage_details_all.csv` | 1,380 | All predicted regions, including those removed by the packaging-module filter (`packaging_module`, `in_primary`) |
| `prophage_clusters.csv` | 130 | Cluster membership: representative region and full member list |
| `cluster_summary.csv` | 130 | Per-cluster statistics: region and genome counts, STs and clades spanned, length range, GC, confidence tiers, best reference phage |
| `species_clusters_95.csv` | 399 | 95% identity / 50% coverage species-level demarcation for each region |
| `protein_gene_families.csv` | 1,447 | CD-HIT protein families with prevalence across regions, genomes and clades, representative product and functional category flags |
| `prophage_cargo_summary.csv` | 399 | Cargo annotation per region: defence, virulence, AMR and MGE gene counts, defence categories, CRISPR arrays, complete defence systems |
| `defensefinder_systems.csv` | 18 | DefenceFinder/CasFinder/RM system calls with wholeness, gene counts and member proteins |
| `statistical_tests.csv` | 9 | Kruskal-Wallis, Fisher, Spearman and Mantel tests with BH-adjusted q values |
| `manuscript_numbers.json` | - | Single source of truth for every statistic quoted in the manuscript (60 top-level entries) |
| `prophage_sequences.fasta` | 399 | Nucleotide sequences of the primary-catalogue regions (18.4 MB) |

`scripts/README.txt` gives the intended order of the pipeline and the external
program versions.

## Requirements

Python 3.10 or newer:

```
pip install -r requirements.txt
```

External programs used for the primary analysis (not required to work with the
deposited tables):

| Program | Version |
|---|---|
| BUSCO | 5.4.7 (`clostridia_odb10`) |
| PhiSpy | 5.0.10 |
| BLASTn | 2.12.0 |
| CD-HIT | 4.8.1 |
| DefenseFinder | 3.0.0 (MacSyFinder 2.1.4, models 3.1.0) |
| CasFinder | 3.1.0 |
| HMMER | 3.3.2 |
| Mash | 2.3 |

## Running the code

The scripts were run from a single working directory and expect this layout:

```
<BASE>/
  genomes/            GenBank records of the 241 assemblies
  scaleup_241/
    results/          intermediate and derived tables
    figures/          figure output
    manuscript/       document output
```

`<BASE>` is currently the hard-coded path `/home/bioinfo/cdprophage` in the
scripts. To run elsewhere, rewrite it once:

```bash
BASE="$PWD/work"
grep -rl "/home/bioinfo/cdprophage" scripts | xargs sed -i "s#/home/bioinfo/cdprophage#${BASE}#g"
```

Then create the layout and populate `results/` from the deposited tables:

```bash
mkdir -p "$BASE/scaleup_241"/{results,figures,manuscript}
cp data/*.csv data/manuscript_numbers.json "$BASE/scaleup_241/results/"
```

Two JSON files consumed by the manuscript builder — `busco_summary.json` and
`species_clusters_95_summary.json` — are aggregates of `busco_summary.csv` and
`species_clusters_95.csv` and can be regenerated with `summarize_busco.py` and
`species_clusters_95.py`.

Reproducing the analysis from raw genomes requires the external programs above
and the GenBank records listed in Table S1; steps 1-14 of `scripts/README.txt`
cover that path. Working from the deposited tables alone (steps 10-16) is
enough to regenerate all figures, supplementary tables and manuscript
statistics.

## Notes and limitations

- Nucleotide clusters are single-linkage, so the 130 clusters are a lower bound
  on the true number; the two 630 prophages (phiCD630-1 and phiCD630-2) merge
  into one cluster for this reason.
- The primary catalogue excludes regions lacking a packaging module; the
  excluded 981 candidates are retained in `prophage_details_all.csv` for
  transparency.
- Cargo annotation is keyword-based on RefSeq protein products and is
  deliberately conservative; gene-level hits are reported separately from
  complete systems.
- Clade assignment follows the MLST clade scheme of the manuscript; genome
  metadata reflect the RefSeq release at the time of analysis (see
  `release_date`).

## Data availability

Genome assemblies: NCBI RefSeq, accessions listed in Table S1 of the manuscript.
Derived data: this repository. Prophage region sequences:
`data/prophage_sequences.fasta`.

## Citation

If you use these data or code, please cite the manuscript above. Full citation
details will be added on publication.

## Licences

- Code (`scripts/`): MIT, see `LICENSE`.
- Data (`data/`, `figures/`): Creative Commons Attribution 4.0 International,
  see `LICENSE-DATA`.
