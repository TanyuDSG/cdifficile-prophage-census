#!/usr/bin/env bash
# Run BUSCO v5.4.7 genome mode (clostridia_odb10) on every genome.
set -uo pipefail
export PATH="/home/bioinfo/miniforge3/envs/busco/bin:$PATH"
BASE=/home/bioinfo/cdprophage
OUT=$BASE/scaleup_241/busco_all
LOGS=$BASE/scaleup_241/logs/busco
LINEAGE=$HOME/busco_downloads/lineages/clostridia_odb10
mkdir -p "$OUT" "$LOGS"

run_one() {
    f="$1"
    acc=$(basename "$f" .fna)
    if compgen -G "$OUT/$acc/short_summary*.json" > /dev/null; then
        return 0
    fi
    busco -i "$f" -o "$acc" -m genome -l "$LINEAGE" --offline \
        -c 2 -f --out_path "$OUT" > "$LOGS/$acc.log" 2>&1
    echo "$(date +%T) $acc exit=$?"
}
export -f run_one
export OUT LOGS LINEAGE

find "$BASE/genomes/fasta" -maxdepth 1 -name '*.fna' -print0 \
    | xargs -0 -P 4 -I{} bash -c 'run_one "$@"' _ {}
echo "ALL DONE"
