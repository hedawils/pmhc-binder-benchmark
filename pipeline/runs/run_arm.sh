#!/usr/bin/env bash
# MPNN + AF2 + metrics for ONE experimental arm. No shutdown -- the caller
# owns the instance lifecycle, because four arms now share one box.
#
#   usage:  run_arm.sh <arm> <backbone_dir> <seqs> [max_backbones] [pattern]
#
#   run_arm.sh free ~/pmhc/out            4 200 'pilot_*.pdb'
#   run_arm.sh fold ~/pmhc/out_scaffolded 6 0   'scaffolded_*.pdb'
#
# The pattern is not cosmetic. Output directories accumulate backbones from
# earlier runs -- a hotspot test left 25 hstest_*.pdb in the free arm's
# directory -- and a bare *.pdb would pull them in. Worse, the cap takes a
# sorted list, and hstest_ sorts before pilot_, so an uncapped glob would
# have filled an eighth of the free arm with designs made under different
# settings. Each arm names its own prefix.
#
# Generalised out of run_pipeline.sh, which hardcoded the free arm's paths
# and ended in `shutdown`. Everything else is unchanged, including
# relax_cycles=0 -- see run_pipeline.sh for why that makes pass rates a
# lower bound against the published protocol. Both arms carry the same
# deviation, so the between-arm comparison is unaffected.
#
# Nothing is filtered here. Every design gets every metric; thresholds are
# applied in analysis.
#
# Resumable: dl_binder_design keeps per-arm checkpoints, so re-running after
# an interruption picks up where it stopped rather than redoing GPU work.
set -euo pipefail

ARM=${1:?arm name}
SRC=${2:?backbone dir}
SEQS=${3:-4}
MAX=${4:-0}
PATTERN=${5:-*.pdb}

BASE="$HOME/pmhc"
DLBD="$BASE/dl_binder_design"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"

echo "########## ARM: $ARM  (seqs=$SEQS, cap=${MAX:-none}, glob=$PATTERN) ##########"
echo "started $(date -u)"

# --- backbone selection -----------------------------------------------------
# The free arm is capped so its AF2 budget can be shared with the fold arm.
# Diffusion designs are i.i.d. samples, so any deterministic subset is an
# unbiased subsample -- no need to randomise.
#
# sort -V, not sort: plain sort is lexicographic, which would order
# pilot_0, pilot_1, pilot_10, pilot_100 and make the cap select an
# arbitrary-looking prefix. Version sort takes designs 0..N-1 in the order
# they were generated, which is what "the first 200" should mean.
STAGE="$BASE/${ARM}_backbones"
rm -rf "$STAGE"; mkdir -p "$STAGE"
ALL=$(find "$SRC" -maxdepth 1 -name "$PATTERN" -type f | sort -V)
[ -n "$ALL" ] || { echo "no backbones matching $PATTERN in $SRC" >&2; exit 1; }
# -maxdepth 1 also excludes the traj/ subdirectory RFdiffusion writes, whose
# intermediate *_traj.pdb files are not designs.
if [ "$MAX" -gt 0 ]; then
  echo "$ALL" | head -n "$MAX" | xargs -I{} cp {} "$STAGE/"
else
  echo "$ALL" | xargs -I{} cp {} "$STAGE/"
fi
echo "=== [1/4] backbones selected: $(ls "$STAGE"/*.pdb | wc -l) of $(echo "$ALL" | wc -l) ==="

# --- two-chain conversion ---------------------------------------------------
"$MM" run -n af2 python "$BASE/to_two_chain.py" "$STAGE" "$BASE/${ARM}_2chain" | tail -3

echo
echo "=== [2/4] ProteinMPNN, ${SEQS} sequences per backbone ==="
cd "$DLBD/mpnn_fr"
"$MM" run -n mpnn python dl_interface_design.py \
  -pdbdir "$BASE/${ARM}_2chain" \
  -outpdbdir "$BASE/${ARM}_mpnn" \
  -checkpoint_name "$BASE/${ARM}.mpnn.checkpoint" \
  -relax_cycles 0 -seqs_per_struct "$SEQS" 2>&1 | grep -c 'reported success' \
  | xargs echo "  backbones designed:"

echo
echo "=== [3/4] AF2 initial guess ==="
cd "$DLBD/af2_initial_guess"
"$MM" run -n af2 python predict.py \
  -pdbdir "$BASE/${ARM}_mpnn" \
  -outpdbdir "$BASE/${ARM}_af2" \
  -checkpoint_name "$BASE/${ARM}.af2.checkpoint" \
  -scorefilename "$BASE/${ARM}.sc" 2>&1 | grep -c 'reported success' \
  | xargs echo "  predictions scored:"

echo
# Collector runs in the mpnn env: af2 is pinned to biopython 1.79 for
# AlphaFold's SCOPData import, and 1.79 calls np.int, removed in numpy 1.24.
echo "=== [4/4] join scores with peptide_focus ==="
"$MM" run -n mpnn python "$BASE/scripts/collect_metrics.py" \
  "$BASE/${ARM}.sc" "$BASE/${ARM}_af2" "$BASE/${ARM}_metrics.tsv"

echo
echo "ARM $ARM COMPLETE at $(date -u) -- rows: $(wc -l < "$BASE/${ARM}_metrics.tsv")"
