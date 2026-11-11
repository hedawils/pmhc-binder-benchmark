#!/usr/bin/env bash
# RFdiffusion binder design against NY-ESO-1 / HLA-A*02:01. Runs on the box.
#
#   usage:  run_rfdiffusion.sh <num_designs> <output_prefix>
#
# Contigs: A1-180 is the HLA alpha1/alpha2 platform, C1-9 is SLLMWITQV,
# and 130-130 asks for a fixed-length 130-residue binder. 130 comes from
# the validated designs rather than a guess -- KH46 is 128 resolved
# residues and NY1-B04 is 145.
#
# Hotspots are peptide positions P4, P5, P7, P8. These are the residues
# that KH46, NY1-B04 and both natural TCRs actually contact (measured in
# scripts/calibrate_peptide_focus.py). Note what is NOT here: no MHC helix
# hotspots. Naming only peptide residues biases diffusion toward engaging
# the peptide, which is the behaviour peptide_focus later measures. P9 is
# excluded deliberately -- it is the F-pocket anchor, buried, and contacted
# by none of the five reference structures.
#
# noise_scale 0 is the recommended setting for binder design; non-zero
# noise measurably lowers in-silico success rates for PPI tasks.
set -euo pipefail

N=${1:-5}
PREFIX=${2:-smoke}

mkdir -p "$HOME/pmhc/out"

# The package is pip-installed into dist-packages and resolves checkpoints
# relative to its own location, not to /app/RFdiffusion. Mount the weights
# where it actually looks rather than where the repo keeps them.
MODELS=/usr/local/lib/python3.9/dist-packages/models

docker run --rm --gpus all \
  -v "$HOME/pmhc/RFdiffusion/models:${MODELS}" \
  -v "$HOME/pmhc/targets:/app/targets" \
  -v "$HOME/pmhc/out:/app/out" \
  rfdiffusion:fixed \
  inference.output_prefix="/app/out/${PREFIX}" \
  inference.input_pdb=/app/targets/c1202_kf9.pdb \
  'contigmap.contigs=[A1-180/0 C1-9/0 130-130]' \
  'ppi.hotspot_res=[C1,C4,C7,C8]' \
  inference.num_designs="$N" \
  denoiser.noise_scale_ca=0 \
  denoiser.noise_scale_frame=0
