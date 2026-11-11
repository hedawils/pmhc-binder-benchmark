#!/usr/bin/env bash
# Arm 5: partial diffusion from Liu et al.'s released privileged scaffolds.
#
#   usage:  run_scaffold_arm.sh <docked_dir> <out_dir> [partial_T] [max]
#
# WHY THIS ARM IS FIRST IN LINE
#
# The expensive, low-yield step in de novo pMHC binder design is finding a
# backbone that contacts the peptide while avoiding the MHC -- the published
# work says so explicitly, and it is why they built a scaffold set at all.
# Their 152 scaffolds are that step's output, released. Starting from them
# skips it entirely, which is the only version of this experiment that fits
# one GPU.
#
# WHY partial_T IS LOW HERE AND HIGHER IN THE KH46 TRANSFER ARM
#
# partial_T governs how far a design may depart from its input pose, and the
# right value depends on how good the input already is. Measured statically
# with the same scorer that grades every design:
#
#   KH46 placed on MART-1      peptide_focus 0.195   BELOW the 0.236 cutoff
#   Liu scaffolds on NY-ESO-1  median        0.354   138/140 ABOVE it
#
# The KH46 transfer starts below the filter and has to rebuild an interface
# over a cavity, so it sweeps 10/20/30. These start above it, so the pose is
# the asset and the job is to keep it: less noising, not more. partial_T=10
# is also ~2x cheaper per design than 20.
#
# WHAT THIS ARM CANNOT CLAIM
#
# peptide_focus measures buried AREA. These binders' side chains were designed
# for other peptides, so a high static focus means the geometry is
# peptide-centred, NOT that the chemistry complements SLLMWITQV. Whether that
# geometry survives resequencing is the question, and the funnel answers it.
set -uo pipefail

DOCKED=${1:?docked complex dir}
OUT_DIR=${2:?output dir}
PT=${3:-10}
MAX=${4:-0}

BASE="$HOME/pmhc"
mkdir -p "$OUT_DIR"

ALL=$(find "$DOCKED" -maxdepth 1 -name 'scaf*.pdb' -type f | sort -V)
[ -n "$ALL" ] || { echo "no docked complexes in $DOCKED" >&2; exit 1; }
if [ "$MAX" -gt 0 ]; then ALL=$(echo "$ALL" | head -n "$MAX"); fi
N=$(echo "$ALL" | wc -l | tr -d ' ')

echo "=== arm 5: $N docked scaffolds, partial_T=$PT ==="

# One RFdiffusion call per scaffold: partial diffusion needs a contig exactly
# as long as its input, and each scaffold's binder is a different length
# (54-118), so they cannot be batched. Already-finished designs are skipped so
# the arm resumes after a restart.
done_n=0; fail=0
for pdb in $ALL; do
  stem=$(basename "$pdb" .pdb)
  if [ -f "$OUT_DIR/liu${stem}_0.pdb" ]; then
    done_n=$((done_n + 1)); continue
  fi
  bash "$BASE/run_partial_transfer.sh" 1 "liu${stem}" "$pdb" "$PT" "$OUT_DIR" \
    >/dev/null 2>&1 || { echo "  FAILED $stem"; fail=$((fail + 1)); continue; }
  done_n=$((done_n + 1))
  [ $((done_n % 20)) -eq 0 ] && echo "  $done_n/$N at $(date -u +%H:%M)"
done

echo "  diffused $done_n of $N ($fail failed)"
echo "  backbones: $(ls "$OUT_DIR"/liu*_*.pdb 2>/dev/null | wc -l)"
echo "ARM 5 DIFFUSION COMPLETE -> $OUT_DIR"
