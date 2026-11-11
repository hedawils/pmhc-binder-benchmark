#!/usr/bin/env bash
# MPNN + AF2 over all RFdiffusion backbones. Runs on the box, in tmux.
#
#   usage:  run_pipeline.sh [seqs_per_backbone]   (default 4)
#
# Nothing is filtered here. Every design gets every metric written to the
# table, and thresholds are applied later in analysis. Discarding failures
# during the run would make the filter-orthogonality question unanswerable
# without repeating ~22 GPU-hours.
#
# DEVIATION FROM THE PUBLISHED PROTOCOL: relax_cycles=0. Bennett et al. use
# FastRelax cycles between MPNN passes, which improves sequences but costs
# CPU time we do not have before Monday. Our pass rates are therefore a
# LOWER BOUND relative to the published pipeline. This must be stated in
# any writeup -- it is a real difference, not a detail.
set -euo pipefail

SEQS=${1:-4}
BASE="$HOME/pmhc"
DLBD="$BASE/dl_binder_design"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"

echo "=== [1/4] convert $(ls "$BASE"/out/pilot_*.pdb 2>/dev/null | wc -l) backbones to two-chain ==="
"$MM" run -n af2 python "$BASE/to_two_chain.py" "$BASE/out" "$BASE/pilot_2chain" | tail -3

echo
echo "=== [2/4] ProteinMPNN, ${SEQS} sequences per backbone ==="
cd "$DLBD/mpnn_fr"
"$MM" run -n mpnn python dl_interface_design.py \
  -pdbdir "$BASE/pilot_2chain" \
  -outpdbdir "$BASE/pilot_mpnn" \
  -checkpoint_name "$BASE/mpnn.checkpoint" \
  -relax_cycles 0 -seqs_per_struct "$SEQS" 2>&1 | grep -c 'reported success' \
  | xargs echo "  backbones designed:"

echo
echo "=== [3/4] AF2 initial guess ==="
cd "$DLBD/af2_initial_guess"
"$MM" run -n af2 python predict.py \
  -pdbdir "$BASE/pilot_mpnn" \
  -outpdbdir "$BASE/pilot_af2" \
  -checkpoint_name "$BASE/af2.checkpoint" \
  -scorefilename "$BASE/pilot.sc" 2>&1 | grep -c 'reported success' \
  | xargs echo "  predictions scored:"

echo
# Collector runs in the MPNN env, not af2: af2 is pinned to biopython 1.79
# for AlphaFold's SCOPData import, and 1.79 calls np.int, which numpy >=1.24
# removed. The mpnn env has biopython 1.88, where both work.
echo "=== [4/4] join scores with peptide_focus ==="
"$MM" run -n mpnn python "$BASE/scripts/collect_metrics.py" \
  "$BASE/pilot.sc" "$BASE/pilot_af2" "$BASE/pilot_metrics.tsv"

echo
echo "PIPELINE COMPLETE"
wc -l < "$BASE/pilot_metrics.tsv" | xargs echo "rows in metrics table:"

# Stop (not terminate) the box: compute billing ends, EBS and all results
# are preserved. Guard set at launch via instance-initiated-shutdown-behavior.
echo "stopping instance to end compute billing; results are on the EBS volume"
sudo shutdown -h +5
