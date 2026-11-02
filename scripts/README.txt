Analysis scripts - C. difficile prophage census (241 complete genomes)
=====================================================================

All scripts were run with Python 3.10 on Ubuntu 22.04. External programs:
BUSCO 5.4.7 (clostridia_odb10, 247 markers), PhiSpy 5.0.10, BLASTn 2.12.0,
CD-HIT 4.8.1, DefenseFinder 3.0.0 with MacSyFinder 2.1.4, HMMER 3.3.2,
DefenseFinder models 3.1.0, CasFinder 3.1.0 and Mash 2.3.

Suggested order
---------------
1.  build_metadata_full.py        genome list, MLST, clade and metadata
2.  run_busco_all.sh              BUSCO v5.4.7 genome mode, clostridia_odb10
    summarize_busco.py            per-genome completeness summary (busco_summary)
3.  run_phispy_all.py             PhiSpy prediction for every genome
    run_phispy_missing.sh         restart helper for unfinished genomes
4.  annotate_prophages_all.py     parse PhiSpy output, structural genes, QC
5.  curate_catalogue.py           packaging-competent primary catalogue
6.  cluster_prophages_all.py      reciprocal-coverage nucleotide clusters
    species_clusters_95.py        95% identity / 50% coverage demarcation
7.  run_pipeline_curated.sh       end-to-end orchestration of steps 3-6
8.  run_defensefinder_all.py      DefenseFinder, CasFinder and RM models
9.  annotate_cargo_all.py         keyword cargo and chromosome baseline
10. diversity_stats_all.py        diversity, accumulation, rarefaction
11. posthoc_all.py                Mantel, partial Mantel, Fisher, Wilcoxon
12. protein_network_all.py        CD-HIT families and gene-content clusters
13. make_mash_all.py              Mash distances and NJ ordering
14. collect_numbers.py            single source of numbers for the manuscript
15. make_figures_all.py           Figures 1-6 and S1-S2 (working versions)
    make_figures_all_journal.py   journal versions: Arial 8 pt, 600 dpi, panel tags
16. make_supp_tables.py           Supplementary Tables S1-S8

The manuscript DOCX and the submission documents (cover letter, point-by-point
response, checklist, validation report) were assembled with separate scripts
that are not part of this repository. Every number they contain comes from
manuscript_numbers.json (step 14) and from the tables and figures produced
above, so the analysis reported in the paper is reproducible from this
directory alone.

Inputs
------
The analysis starts from the NCBI RefSeq assembly summary and the GenBank
records of the 241 complete genomes. Genome accessions are listed in Table S1.
Paths in the scripts point to the working directory used for the analysis and
must be adjusted before rerunning.

The data/ folder holds the derived tables needed to reproduce the manuscript
without rerunning the external programs: the genome panel (genome_metadata,
genome_summary_analysis, busco_summary), the region catalogue
(prophage_details, prophage_details_all), cluster tables (cluster_summary,
prophage_clusters, species_clusters_95), gene content
(protein_gene_families), cargo (prophage_cargo_summary), defence systems
(defensefinder_systems), statistics (statistical_tests) and the region
sequences (prophage_sequences.fasta).
