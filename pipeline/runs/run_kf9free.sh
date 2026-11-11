#!/usr/bin/env bash
# Free-diffusion control arm for the KF9 target: the missing baseline that
# turns "borrowing beats free diffusion" from n=1 target into n=2.
#
# 150 backbones x 4 sequences = 600 designs. Sized to LAND, not to maximise n:
# a 10am talk needs results with hours of margin, not minutes. Power is still
# p ~ 0.001 if free diffusion suffers the same difficulty penalty on KF9 that
# the borrowed arm did.
# exactly so the two free arms are directly comparable as well.
# Sized from a measured 3.7 min/backbone, so ~6h generation + ~2h scoring.
#
# Identical settings to the NY-ESO-1 free arm (130-residue binder, full 50
# denoising steps) except the target and its hotspots, so the two free arms
# are themselves comparable.
set -uo pipefail
BASE="$HOME/pmhc"; cd "$BASE"
export PMHC_PEPTIDE=KAYNVTQAF
echo "=== KF9 FREE ARM started $(date -u '+%F %T UTC') ==="
bash "$BASE/run_rfdiffusion_kf9.sh" 150 kf9free || echo "KF9 FREE DIFFUSION FAILED"
N=$(ls "$BASE"/out/kf9free_*.pdb 2>/dev/null | wc -l)
echo "=== backbones: $N  $(date -u '+%T') ==="
if [ "$N" -gt 0 ]; then
  bash "$BASE/run_arm.sh" kf9free "$BASE/out" 4 0 'kf9free_*.pdb' || echo "KF9 FREE ARM FAILED"
else
  echo "no backbones; skipping scoring"
fi
echo "=== KF9 FREE ARM COMPLETE $(date -u '+%F %T UTC') ==="
wc -l "$BASE/kf9free_metrics.tsv" 2>/dev/null || echo "no metrics"
