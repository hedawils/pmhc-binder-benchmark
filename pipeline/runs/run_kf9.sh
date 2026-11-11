#!/usr/bin/env bash
# Cross-locus arm: Liu's privileged scaffolds against SARS-CoV-2 N KF9
# (KAYNVTQAF) on HLA-C*12:02 -- a target that is nobody's cognate, so no
# cohort is scored on home ground and "cross-locus" stops being collinear
# with "the PHOX2B cohort".
#
# partial_T is 10, matching arm 3, so the ONLY change from arm 3 is the
# target. Note the caveat: arm 3's docked poses started above the focus
# cutoff (median 0.354, 138/140 above) and only needed perturbing. KF9's
# start lower (median 0.220, 46/147 above), so T=10 may under-remodel here.
# Matching arm 3 keeps the comparison clean; a T sweep is the follow-up.
set -uo pipefail
BASE="$HOME/pmhc"; cd "$BASE"
export PMHC_PEPTIDE=KAYNVTQAF        # to_two_chain.py's boundary assertion
echo "=== ARM KF9 started $(date -u '+%F %T UTC') ==="
bash "$BASE/run_scaffold_arm_kf9.sh" "$BASE/docked_kf9" "$BASE/out_kf9" 10 \
  || echo "KF9 DIFFUSION FAILED"
N=$(ls "$BASE"/out_kf9/kf9*_*.pdb 2>/dev/null | wc -l)
echo "=== KF9 backbones: $N  $(date -u '+%T') ==="
if [ "$N" -gt 0 ]; then
  bash "$BASE/run_arm.sh" kf9 "$BASE/out_kf9" 4 0 'kf9*_*.pdb' || echo "KF9 ARM FAILED"
else
  echo "no backbones; skipping scoring"
fi
echo "=== ARM KF9 COMPLETE $(date -u '+%F %T UTC') ==="
wc -l "$BASE/kf9_metrics.tsv" 2>/dev/null || echo "no metrics file"
