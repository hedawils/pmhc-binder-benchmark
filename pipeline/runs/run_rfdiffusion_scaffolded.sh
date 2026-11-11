#!/usr/bin/env bash
# Fold-conditioned (scaffold-guided) binder design against NY-ESO-1/HLA-A*02:01.
#
#   usage:  run_rfdiffusion_scaffolded.sh <num_designs> <output_prefix>
#
# WHY THIS EXISTS
#
# The published TCR-mimic work against this target did not do free binder
# design. It screened 50 RFdiffusion monomer folds, picked a 99-residue
# four-helix bundle, and fold-conditioned that scaffold onto the pMHC.
# Conditioning on a known-good topology is a materially easier problem than
# inventing one, so a low free-design pass rate is not evidence of a
# misconfiguration -- it is a different, harder question.
#
# Hotspots are unchanged: the paper used Met4, Trp5, Thr7, Gln8, which is
# what we derived independently from KH46's contacts and already use.
#
# Scaffolds come from make_scaffolds.sh: unconditional monomers generated
# at 100-130 residues, bracketing the free arm's fixed 130. The shipped PPI
# library is deliberately NOT used -- its folds are 50-65 residues, below
# the range of every binder validated on this target (99 published, 128
# KH46, 145 NY1-B04), so conditioning on it would move length and
# conditioning at the same time and confound the comparison.
#
# Conditioning on KH46's own fold would inflate the pass rate by measuring
# whether RFdiffusion can rebuild a known solution rather than find one,
# so that is not used either.
#
# Note there is no contigmap.contigs: the scaffold library defines the
# binder, and the target is supplied through scaffoldguided.target_path.
set -euo pipefail

#   usage:  run_rfdiffusion_scaffolded.sh <n> <prefix> [scaffold_dir] [out_dir]
#                                         [sampled_insertion]
#
# The scaffold directory is a parameter because two arms use this path with
# different libraries: the fold arm with generic matched-length monomers, and
# the KH46 control with the blueprint of a binder already validated on this
# target -- circular as a RESULT, but a legitimate plumbing check.
#
# SAMPLED_INSERTION, AND WHY THE DEFAULT OF 0 IS WRONG FOR THE FOLD ARM
#
# The scaffoldguided defaults are sampled_insertion: 0, sampled_N: 0,
# sampled_C: 0. mask_loops is True, but that only allows loop CONFORMATION to
# be rebuilt -- it does not change how many residues are in a loop. With
# insertion at 0, every design conditioned on one scaffold comes out at
# exactly that scaffold's length, varying only by diffusion noise. Measured,
# not assumed: all ten KH46-control backbones are exactly 128 residues, the
# blueprint's length.
#
# The published work selected a 99-residue fold and reported backbones of
# 99-129 aa. That 30-residue spread IS insertion, so conditioning without it
# is not their protocol. It also fixes a length confound against the free
# arm: scaffolds generated at 100-130 (median ~115) plus a uniform 0-30
# insertion centre at ~130, the free arm's fixed length, rather than sitting
# ~15 residues below it.
N=${1:-20}
PREFIX=${2:-scaffolded}
BASE="$HOME/pmhc"
SCAFFOLD_DIR=${3:-$BASE/scaffold_folds}
OUT_DIR=${4:-$BASE/out_scaffolded}
INSERTION=${5:-0}
RFD="$BASE/RFdiffusion"
PKG=/usr/local/lib/python3.9/dist-packages
MODELS="$PKG/models"
# See make_scaffolds.sh: the default config names an example PDB that Hydra
# loads regardless of mode, and it lives at a package-relative path.
EXAMPLES="$PKG/examples"

mkdir -p "$OUT_DIR"

# The library must exist first.
if [ ! -d "$SCAFFOLD_DIR" ] || [ -z "$(ls -A "$SCAFFOLD_DIR" 2>/dev/null)" ]; then
  echo "no scaffold library at $SCAFFOLD_DIR; run make_scaffolds.sh first" >&2
  exit 1
fi

# One-time: secondary structure and block adjacency for the TARGET.
# Scaffold-guided mode needs these for the target as well as the binder.
if [ ! -f "$BASE/target_folds/9min_target_a1a2_pep_ss.pt" ]; then
  mkdir -p "$BASE/target_folds"
  docker run --rm --gpus all \
    -v "$RFD/models:${MODELS}" \
    -v "$RFD/examples:${EXAMPLES}" \
    -v "$BASE/targets:/app/targets" \
    -v "$BASE/target_folds:/app/target_folds" \
    -v "$RFD/helper_scripts:/app/helper_scripts" \
    --entrypoint python3.9 rfdiffusion:fixed \
    /app/helper_scripts/make_secstruc_adj.py \
    --input_pdb /app/targets/9min_target_a1a2_pep.pdb \
    --out_dir /app/target_folds
fi

docker run --rm --gpus all \
  -v "$RFD/models:${MODELS}" \
  -v "$RFD/examples:${EXAMPLES}" \
  -v "$BASE/targets:/app/targets" \
  -v "$BASE/target_folds:/app/target_folds" \
  -v "$SCAFFOLD_DIR:/app/scaffold_folds" \
  -v "$OUT_DIR:/app/out" \
  rfdiffusion:fixed \
  inference.output_prefix="/app/out/${PREFIX}" \
  scaffoldguided.scaffoldguided=True \
  scaffoldguided.target_pdb=True \
  scaffoldguided.target_path=/app/targets/9min_target_a1a2_pep.pdb \
  scaffoldguided.target_ss=/app/target_folds/9min_target_a1a2_pep_ss.pt \
  scaffoldguided.target_adj=/app/target_folds/9min_target_a1a2_pep_adj.pt \
  scaffoldguided.scaffold_dir=/app/scaffold_folds \
  scaffoldguided.sampled_insertion="$INSERTION" \
  'ppi.hotspot_res=[C4,C5,C7,C8]' \
  inference.num_designs="$N" \
  denoiser.noise_scale_ca=0 \
  denoiser.noise_scale_frame=0
