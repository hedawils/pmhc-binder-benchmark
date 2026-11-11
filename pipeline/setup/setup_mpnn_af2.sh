#!/usr/bin/env bash
# Install ProteinMPNN + AF2 initial-guess (dl_binder_design) on the box.
#
# Runs alongside the RFdiffusion job. Almost all of this is download and
# dependency solving, so GPU contention with the running diffusion is
# negligible. Doing it now rather than at 23:30 means any environment
# breakage surfaces while there is still a working day to react to it --
# RFdiffusion alone threw two independent environment failures today.
#
# micromamba rather than conda: a single static binary, no bootstrap, and
# a solver that finishes in minutes instead of grinding.
set -euo pipefail

BASE="$HOME/pmhc/dl_binder_design"
cd "$BASE"

echo "=== [1/4] ProteinMPNN (submodule came up empty) ==="
if [ -z "$(ls -A mpnn_fr/ProteinMPNN 2>/dev/null)" ]; then
  rm -rf mpnn_fr/ProteinMPNN
  git clone -q https://github.com/dauparas/ProteinMPNN.git mpnn_fr/ProteinMPNN
fi
echo "  weights: $(ls mpnn_fr/ProteinMPNN/vanilla_model_weights/*.pt 2>/dev/null | wc -l) checkpoints"

echo "=== [2/4] AF2 params (~5GB) ==="
PARAMS="$BASE/af2_initial_guess/model_weights/params"
mkdir -p "$PARAMS"
if [ ! -f "$PARAMS/params_model_1_ptm.npz" ]; then
  cd "$PARAMS"
  wget -q --show-progress https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar
  tar --extract --file=alphafold_params_2022-12-06.tar
  rm -f alphafold_params_2022-12-06.tar
  cd "$BASE"
fi
echo "  params_model_1_ptm.npz: $(ls -la "$PARAMS/params_model_1_ptm.npz" 2>/dev/null | awk '{print $5}') bytes"

echo "=== [3/4] micromamba ==="
cd "$HOME"
if [ ! -x "$HOME/bin/micromamba" ]; then
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba
fi
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"
echo "  $($MM --version)"

echo "=== [4/4] environments ==="
for env_file in proteinmpnn_fastrelax af2_binder_design; do
  name=$($MM env list 2>/dev/null | awk '{print $1}' | grep -x "$env_file" || true)
  if [ -z "$name" ]; then
    echo "  creating $env_file ..."
    $MM create -y -q -f "$BASE/include/${env_file}.yml" 2>&1 | tail -3
  else
    echo "  $env_file already present"
  fi
done

echo
echo "SETUP COMPLETE"
$MM env list
