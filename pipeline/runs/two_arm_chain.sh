#!/usr/bin/env bash
# Unattended overnight orchestration for the two-arm comparison.
#
# Waits for free-design diffusion to finish, builds a matched-length scaffold
# library, runs fold-conditioned diffusion, then takes both arms through
# MPNN + AF2 + metrics and stops the box.
#
# WHY TWO ARMS
#
# The published TCR-mimic work found hits from ~600 fold-conditioned designs
# (100 scaffolds x 6 sequences). The free-design route that produced a
# cryo-EM-validated binder on this target needed ~5,500. If that ~9x
# difference in sample efficiency is real, fold conditioning is the only
# route with a plausible chance of a hit at the scale available here.
#
# Running both, with the same target, hotspots, MPNN settings and AF2
# settings, makes conditioning the single moving variable -- a controlled
# comparison rather than one arm plus a literature citation. Scaffold length
# is matched deliberately (see make_scaffolds.sh) so it is not a second one.
#
# ORDERING
#
# The fold arm goes through AF2 first. Both arms write metrics incrementally,
# but a hit, if there is one, is likelier in the fold arm, so it gets the
# uncontended GPU hours. A partial free arm still yields a usable pass-rate
# estimate; a partial fold arm might simply miss its hits.
#
# The free arm is capped at 200 of 400 backbones. Pass-rate precision barely
# improves past 200, and the ~11 GPU-hours freed buy an entire second arm --
# a categorical gain over a tighter interval on a number that may well be
# zero.
#
# FAILURE ISOLATION
#
# Deliberately no `set -e` around the stages. If the fold arm fails -- it is
# the newer path -- the free arm must still run to completion rather than
# being taken down with it.
set -uo pipefail

BASE="$HOME/pmhc"
LOG="$BASE/two_arm.log"
exec > >(tee -a "$LOG") 2>&1

echo "######################################################################"
echo "# two-arm chain armed at $(date -u)"
echo "######################################################################"

step() { echo; echo ">>> $* :: $(date -u)"; }

# --- wait for free diffusion ------------------------------------------------
step "waiting for free diffusion (tmux session 'pilot')"
while tmux has-session -t pilot 2>/dev/null; do sleep 60; done
FREE_N=$(ls "$BASE"/out/pilot_*.pdb 2>/dev/null | wc -l)
echo "free diffusion finished; backbones: $FREE_N"

# --- fold arm: scaffolds, then conditioned diffusion ------------------------
FOLD_OK=1
step "matched-length scaffold library (50 monomers, 100-130 res)"
bash "$BASE/make_scaffolds.sh" 50 100-130 || { echo "SCAFFOLDS FAILED"; FOLD_OK=0; }

if [ "$FOLD_OK" -eq 1 ]; then
  step "fold-conditioned diffusion (100 backbones)"
  bash "$BASE/run_rfdiffusion_scaffolded.sh" 100 scaffolded \
    || { echo "SCAFFOLDED DIFFUSION FAILED"; FOLD_OK=0; }
fi

FOLD_N=$(ls "$BASE"/out_scaffolded/scaffolded_*.pdb 2>/dev/null | wc -l)
echo "fold-conditioned backbones: $FOLD_N"
[ "$FOLD_N" -gt 0 ] || FOLD_OK=0

# --- AF2: fold arm first, then free ----------------------------------------
# Each arm names its own prefix. The free arm's directory also holds 25
# hstest_*.pdb from an earlier hotspot test, and those sort before pilot_,
# so a bare *.pdb glob would have put them at the front of the capped 200.
if [ "$FOLD_OK" -eq 1 ]; then
  step "ARM fold -- 6 seqs per backbone"
  bash "$BASE/run_arm.sh" fold "$BASE/out_scaffolded" 6 0 'scaffolded_*.pdb' \
    || echo "FOLD ARM FAILED"
else
  echo; echo "skipping fold arm: no backbones"
fi

step "ARM free -- 4 seqs per backbone, capped at 200"
bash "$BASE/run_arm.sh" free "$BASE/out" 4 200 'pilot_*.pdb' \
  || echo "FREE ARM FAILED"

# --- summary ----------------------------------------------------------------
echo
echo "######################################################################"
echo "# two-arm chain done at $(date -u)"
for arm in fold free; do
  f="$BASE/${arm}_metrics.tsv"
  if [ -f "$f" ]; then
    echo "#   $arm: $(( $(wc -l < "$f") - 1 )) designs scored"
  else
    echo "#   $arm: NO METRICS"
  fi
done
echo "######################################################################"

# Stop, not terminate: compute billing ends, EBS and all results survive.
echo "stopping instance in 5 minutes"
sudo shutdown -h +5
