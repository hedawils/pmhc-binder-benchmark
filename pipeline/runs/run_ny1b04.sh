#!/usr/bin/env bash
# Arm 1c: ny1b04 partially diffused onto our standard NY-ESO-1 target.
# Second self-transfer source, so arm 1b's ceiling has n=2 instead of n=1.
# Deliberately does NOT schedule a shutdown; the box stays up.
set -uo pipefail
BASE="$HOME/pmhc"
cd "$BASE"
echo "=== ARM 1c ny1b04 self-transfer  started $(date -u '+%F %T UTC') ==="

for PT in 10 20 30; do
  echo "--- partial_T=$PT  $(date -u '+%T') ---"
  bash "$BASE/run_partial_transfer.sh" 17 "ny1b04pt${PT}" \
       "$BASE/control2/complex.pdb" "$PT" "$BASE/out_ny1b04" \
    || echo "NY1B04 DIFFUSION FAILED at partial_T=$PT"
done

N=$(ls "$BASE"/out_ny1b04/ny1b04pt*.pdb 2>/dev/null | wc -l)
echo "=== backbones produced: $N  $(date -u '+%T') ==="
if [ "$N" -gt 0 ]; then
  bash "$BASE/run_arm.sh" ny1b04 "$BASE/out_ny1b04" 6 0 'ny1b04pt*.pdb' \
    || echo "NY1B04 ARM FAILED"
else
  echo "no backbones; skipping scoring"
fi
echo "=== ARM 1c COMPLETE  $(date -u '+%F %T UTC') ==="
wc -l "$BASE/ny1b04_metrics.tsv" 2>/dev/null || echo "no metrics file"
