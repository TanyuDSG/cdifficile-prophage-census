#!/usr/bin/env bash
set -euo pipefail
BASE=/home/bioinfo/cdprophage
RES=$BASE/scaleup_241/results
LOG=$BASE/scaleup_241/logs
mkdir -p "$RES" "$LOG"
export LD_LIBRARY_PATH="$BASE/tools/root/usr/lib/ncbi-blast+:$BASE/tools/mashlibs/root/usr/lib/x86_64-linux-gnu:$BASE/tools/root/usr/lib/x86_64-linux-gnu"

echo "[$(date +%F' '%T)] waiting for PhiSpy"
while true; do
  n=$(find "$BASE/phispy" -maxdepth 2 -name prophage_coordinates.tsv | wc -l)
  echo "  $n/241"
  [ "$n" -ge 241 ] && break
  sleep 300
done

echo "[$(date +%F' '%T)] annotate prophages"
python3 "$BASE/scaleup_241/scripts/annotate_prophages_all.py" > "$LOG/annotate_prophages.log" 2>&1

echo "[$(date +%F' '%T)] cluster prophages"
python3 "$BASE/scaleup_241/scripts/cluster_prophages_all.py" > "$LOG/cluster_prophages.log" 2>&1

echo "[$(date +%F' '%T)] annotate cargo"
python3 "$BASE/scaleup_241/scripts/annotate_cargo_all.py" > "$LOG/annotate_cargo.log" 2>&1

echo "[$(date +%F' '%T)] CD-HIT"
"$BASE/tools/root/usr/bin/cd-hit" -i "$RES/prophage_proteins.fasta" \
  -o "$RES/prophage_proteins_nr.fasta" -c 0.5 -n 5 -M 12000 -T 8 -d 0 \
  > "$LOG/cdhit.log" 2>&1

echo "[$(date +%F' '%T)] DefenseFinder"
python3 "$BASE/scaleup_241/scripts/run_defensefinder_all.py" > "$LOG/defensefinder.log" 2>&1

echo "[$(date +%F' '%T)] protein network"
python3 "$BASE/scaleup_241/scripts/protein_network_all.py" > "$LOG/protein_network.log" 2>&1

echo "[$(date +%F' '%T)] diversity stats"
python3 "$BASE/scaleup_241/scripts/diversity_stats_all.py" > "$LOG/diversity_stats.log" 2>&1

echo "[$(date +%F' '%T)] done"
