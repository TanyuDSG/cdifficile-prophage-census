#!/usr/bin/env bash
# Curated 241-genome pipeline: packaging-positive primary catalogue,
# reciprocal 70%/50% nucleotide clustering, then cargo/defence/diversity/figures.
set -euo pipefail
BASE=/home/bioinfo/cdprophage
RES=$BASE/scaleup_241/results
LOG=$BASE/scaleup_241/logs
mkdir -p "$RES" "$LOG"
export LD_LIBRARY_PATH="$BASE/tools/root/usr/lib/ncbi-blast+:$BASE/tools/mashlibs/root/usr/lib/x86_64-linux-gnu:$BASE/tools/root/usr/lib/x86_64-linux-gnu"
export PATH="$BASE/tools/root/usr/bin:$PATH"

echo "[$(date +%F' '%T)] curate catalogue"
python3 "$BASE/scaleup_241/scripts/curate_catalogue.py" > "$LOG/curate.log" 2>&1

echo "[$(date +%F' '%T)] cluster primary regions"
python3 "$BASE/scaleup_241/scripts/cluster_prophages_all.py" > "$LOG/cluster_curated.log" 2>&1

echo "[$(date +%F' '%T)] annotate cargo"
python3 "$BASE/scaleup_241/scripts/annotate_cargo_all.py" > "$LOG/annotate_cargo_curated.log" 2>&1

echo "[$(date +%F' '%T)] CD-HIT"
"$BASE/tools/root/usr/bin/cd-hit" -i "$RES/prophage_proteins.fasta" \
  -o "$RES/prophage_proteins_nr.fasta" -c 0.5 -n 3 -aS 0.5 -M 12000 -T 8 -d 0 \
  > "$LOG/cdhit_curated.log" 2>&1

echo "[$(date +%F' '%T)] DefenseFinder"
python3 "$BASE/scaleup_241/scripts/run_defensefinder_all.py" > "$LOG/defensefinder_curated.log" 2>&1

echo "[$(date +%F' '%T)] protein network"
python3 "$BASE/scaleup_241/scripts/protein_network_all.py" > "$LOG/protein_network_curated.log" 2>&1

echo "[$(date +%F' '%T)] diversity stats"
python3 "$BASE/scaleup_241/scripts/diversity_stats_all.py" > "$LOG/diversity_curated.log" 2>&1

echo "[$(date +%F' '%T)] posthoc"
python3 "$BASE/scaleup_241/scripts/posthoc_all.py" > "$LOG/posthoc_curated.log" 2>&1

echo "[$(date +%F' '%T)] figures"
python3 "$BASE/scaleup_241/scripts/make_figures_all.py" > "$LOG/figures_curated.log" 2>&1

echo "[$(date +%F' '%T)] supplementary tables"
python3 "$BASE/scaleup_241/scripts/make_supp_tables.py" > "$LOG/supp_curated.log" 2>&1

echo "[$(date +%F' '%T)] collect numbers"
python3 "$BASE/scaleup_241/scripts/collect_numbers.py" > "$LOG/collect_curated.log" 2>&1

echo "[$(date +%F' '%T)] done"
