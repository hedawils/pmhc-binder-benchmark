#!/usr/bin/env bash
# Partial diffusion of a VALIDATED binder architecture onto a NEW pMHC target.
#
#   usage:  run_partial_transfer.sh <n> <prefix> <input_complex.pdb> \
#                                   [partial_T] [out_dir]
#
# WHY PARTIAL DIFFUSION RATHER THAN FOLD CONDITIONING
#
# Fold conditioning passes secondary structure and block adjacency -- an
# abstraction of the topology. It says "four helices packed like so" and
# nothing about sitting diagonally across the groove at a TCR-like angle.
# For a transfer experiment the docking POSE is the most valuable thing the
# validated binder knows, and conditioning discards it. Partial diffusion
# keeps the coordinates and perturbs them, so the pose survives and
# partial_T tunes how far the design may depart from it.
#
# WHY IT IS NOT CIRCULAR
#
# The architecture is known-good as an ARCHITECTURE, not as a solution to
# the peptide being targeted here. Conditioning on KH46 against its own
# target would measure rebuilding; conditioning on it against a different
# peptide asks whether a pMHC-binding fold transfers. It still carries
# target-CLASS information -- the fold is shaped to sit on a groove -- so
# the honest claim is architecture transfer, not unbiased discovery.
#
# THE INPUT
#
# Partial diffusion starts from coordinates, so the binder must already be
# placed on the new target: chain A = binder, B = MHC, C = peptide, built by
# superimposing the new platform onto 9MIN's and carrying the binder across.
# RFdiffusion requires a contig string exactly as long as the input, so the
# binder length is read from the file rather than assumed -- it cannot
# change, which is this method's one real limitation.
set -euo pipefail

N=${1:-100}
PREFIX=${2:-transfer}
INPUT=${3:?input complex pdb (A=binder, B=MHC, C=peptide)}
PT=${4:-20}

BASE="$HOME/pmhc"
RFD="$BASE/RFdiffusion"
PKG=/usr/local/lib/python3.9/dist-packages
MODELS="$PKG/models"
EXAMPLES="$PKG/examples"
# Explicit out_dir so a partial_T sweep can write several prefixes into one
# directory and be scored by a single run_arm.sh pass.
OUT_DIR=${5:-"$BASE/out_${PREFIX}"}
IN_DIR=$(cd "$(dirname "$INPUT")" && pwd)
IN_FILE=$(basename "$INPUT")

mkdir -p "$OUT_DIR"

# partial_T is relative to diffuser.T, which defaults to 50. RFdiffusion's
# own guidance: partial_T 20 is the equivalent of the 80-of-200 used in
# Vazquez-Torres et al., and they recommend sweeping it. Too low leaves the
# binder clashing with the new groove; too high discards the pose.
#
# For the MART-1 transfer there is no geometry to repair -- prep_transfer.py
# measured ZERO backbone clashes after superposition, and all twelve clashes
# in the placed complex are side-chain/side-chain, which diffusion never sees
# and MPNN repacks for free. So partial_T is not being set to fix sterics.
# It is being set to decide how much INTERFACE REMODELLING is allowed, and the
# answer is not guessable: the new peptide replaces KH46's largest contact
# (Trp5, 32 contacts) with glycine (4), so the binder has to rebuild over a
# cavity. Hence the chain sweeps it rather than asserting a value.
if [ "$PT" -gt 50 ]; then
  echo "partial_T=$PT exceeds diffuser.T=50" >&2
  exit 1
fi

chain_len() {
  awk -v ch="$1" '$1=="ATOM" && substr($0,22,1)==ch && substr($0,13,4)~/CA/ {n++} END {print n+0}' "$INPUT"
}
BINDER_LEN=$(chain_len A)
MHC_LEN=$(chain_len B)
PEP_LEN=$(chain_len C)

for pair in "A:$BINDER_LEN" "B:$MHC_LEN" "C:$PEP_LEN"; do
  if [ "${pair#*:}" -eq 0 ]; then
    echo "input $IN_FILE has no chain ${pair%%:*}; expected A=binder B=MHC C=peptide" >&2
    exit 1
  fi
done

echo "=== partial-diffusion transfer ==="
echo "  input      $IN_FILE  (A=$BINDER_LEN  B=$MHC_LEN  C=$PEP_LEN)"
echo "  partial_T  $PT of 50"
echo "  designs    $N"

# The contig must total exactly the input length: the binder as a fixed
# length, the target chains held by reference.
CONTIGS="[${BINDER_LEN}-${BINDER_LEN}/0 B1-${MHC_LEN}/0 C1-${PEP_LEN}]"
echo "  contigs    $CONTIGS"

docker run --rm --gpus all \
  -v "$RFD/models:${MODELS}" \
  -v "$RFD/examples:${EXAMPLES}" \
  -v "$IN_DIR:/app/in" \
  -v "$OUT_DIR:/app/out" \
  rfdiffusion:fixed \
  inference.output_prefix="/app/out/${PREFIX}" \
  inference.input_pdb="/app/in/${IN_FILE}" \
  "contigmap.contigs=${CONTIGS}" \
  diffuser.partial_T="$PT" \
  'ppi.hotspot_res=[C1,C4,C7,C8]' \
  inference.num_designs="$N" \
  inference.write_trajectory=False \
  denoiser.noise_scale_ca=0 \
  denoiser.noise_scale_frame=0

echo "  backbones: $(ls "$OUT_DIR"/${PREFIX}_*.pdb 2>/dev/null | wc -l)"
echo "TRANSFER DIFFUSION COMPLETE -> $OUT_DIR"
