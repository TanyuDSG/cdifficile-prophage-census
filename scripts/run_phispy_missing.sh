#!/bin/bash
# Rerun PhiSpy on the four genomes with origin-spanning annotations that
# triggered the zero-length ORF bug (patched in PhiSpyModules/makeTest.py).
set -u
BASE=/home/bioinfo/cdprophage
export PATH="$BASE/tools/root/usr/bin:$HOME/.local/bin:$PATH"
export LD_LIBRARY_PATH="$BASE/tools/root/usr/lib/ncbi-blast+:$BASE/tools/root/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH:-}"
ACC=(GCF_037918925.1 GCF_037931695.1 GCF_054790395.1 GCF_060588345.1)

run_one() {
  a="$1"
  out="$BASE/phispy/$a"
  if [ -f "$out/prophage_coordinates.tsv" ]; then echo "skip $a"; return 0; fi
  if [ -f "$out/phispy.log" ]; then mv "$out/phispy.log" "$out/phispy.previous.log"; fi
  phispy -o "$out" --threads 2 --quiet "$BASE/genomes/gbff/$a.gbff" \
    > "$BASE/scaleup_241/logs/phispy_${a}.rerun.log" 2>&1
  echo "$a exit=$? coords=$([ -f "$out/prophage_coordinates.tsv" ] && echo yes || echo no)"
}

for a in "${ACC[@]}"; do run_one "$a" & done
wait
echo DONE
