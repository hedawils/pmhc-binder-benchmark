#!/usr/bin/env bash
#
# Rebuild the de novo pMHC binder design environment from scratch.
#
# Captured from the working pilot instance on 2026-11-02. Every version below
# is the one that actually produced the results, not a "latest" that may have
# drifted since.
#
#   Instance   g5.2xlarge (NVIDIA A10G, 23 GB)
#   OS         Ubuntu 22.04.5 LTS, kernel 6.8.0-1052-aws
#   Driver     580.126.09, CUDA 12.8 on the host
#
# The host CUDA version does not need to match the container's (11.6) --
# only the driver has to be new enough, and 580.x is far beyond the minimum.
#
# Two separate dependency worlds, deliberately kept apart:
#
#   RFdiffusion      Docker. Its SE3Transformer stack (dgl + e3nn + a specific
#                    torch/CUDA pairing) is the single most reliable way to
#                    lose a day to dependency resolution. The official image
#                    pins all of it.
#   MPNN + AF2       micromamba. A static binary with a solver that finishes
#                    in minutes rather than grinding.
#
# Usage:  bash reproduce_env.sh          # full build, ~60-90 min, ~25 GB disk
#         bash reproduce_env.sh verify   # check an existing install only
#
set -euo pipefail

BASE="$HOME/pmhc"
MAMBA_ROOT_PREFIX="$HOME/micromamba"; export MAMBA_ROOT_PREFIX
MM="$HOME/bin/micromamba"

# --- pinned revisions -------------------------------------------------------
RFDIFFUSION_SHA=86507b6          # RosettaCommons/RFdiffusion, 2026-07-15
DLBINDER_SHA=cafa385             # nrbennet/dl_binder_design,  2024-01-10
PROTEINMPNN_SHA=8907e66          # dauparas/ProteinMPNN
AF2_PARAMS=alphafold_params_2022-12-06.tar
MICROMAMBA_VERSION=2.10.0
# ---------------------------------------------------------------------------

step () { echo; echo "=== $* ==="; }

verify () {
  local fail=0 n nv
  step "verification"

  if docker image inspect rfdiffusion:fixed >/dev/null 2>&1; then
    echo "  ok   docker image rfdiffusion:fixed"
    # The numpy pin is the failure that is easiest to reintroduce and hardest
    # to diagnose: torch 1.12.1 is built against the numpy 1.x C ABI, the
    # stock image ships 2.0.2, and the only symptom is "Numpy is not
    # available" raised from inside the diffusion loop.
    nv=$(docker run --rm --entrypoint python3.9 rfdiffusion:fixed \
           -c 'import numpy; print(numpy.__version__)' 2>/dev/null || echo ERROR)
    case "$nv" in
      1.*) echo "  ok   numpy $nv inside container (must be 1.x)" ;;
      *)   echo "  FAIL numpy '$nv' inside container -- must be 1.x"; fail=1 ;;
    esac
  else
    echo "  MISSING  docker image rfdiffusion:fixed"
    fail=1
  fi

  n=$( { ls "$BASE/RFdiffusion/models"/*.pt 2>/dev/null || true; } | wc -l | tr -d " ")
  if [ "$n" -ge 7 ]; then
    echo "  ok   $n RFdiffusion checkpoints"
  else
    echo "  MISSING  RFdiffusion checkpoints (found $n, want >=7)"
    fail=1
  fi

  if [ -f "$BASE/dl_binder_design/af2_initial_guess/model_weights/params/params_model_1_ptm.npz" ]; then
    echo "  ok   AF2 params"
  else
    echo "  MISSING  AF2 params"
    fail=1
  fi

  n=$( { ls "$BASE/dl_binder_design/mpnn_fr/ProteinMPNN/vanilla_model_weights"/*.pt 2>/dev/null || true; } | wc -l | tr -d " ")
  if [ "$n" -ge 1 ]; then
    echo "  ok   $n ProteinMPNN checkpoints"
  else
    echo "  MISSING  ProteinMPNN weights"
    fail=1
  fi

  for e in mpnn af2; do
    if [ -d "$MAMBA_ROOT_PREFIX/envs/$e" ]; then
      echo "  ok   micromamba env '$e'"
    else
      echo "  MISSING  micromamba env '$e'"
      fail=1
    fi
  done

  if nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null; then
    :
  else
    echo "  WARN no GPU visible on host"
  fi

  echo
  if [ "$fail" -eq 0 ]; then
    echo "ENVIRONMENT OK"
  else
    echo "ENVIRONMENT INCOMPLETE"
    return 1
  fi
}

if [ "${1:-}" = "verify" ]; then verify; exit $?; fi

# ---------------------------------------------------------------------------
step "[1/6] host packages: docker + nvidia container toolkit"
if ! command -v docker >/dev/null; then
  sudo apt-get -qq update
  sudo apt-get -qq install -y docker.io git wget curl
  sudo usermod -aG docker "$USER"
  echo "  NOTE: log out and back in for the docker group to take effect,"
  echo "        then re-run this script."
fi
if ! docker info 2>/dev/null | grep -qi nvidia; then
  curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey \
    | sudo gpg --yes --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg
  curl -fsSL https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list \
    | sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' \
    | sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list >/dev/null
  sudo apt-get -qq update && sudo apt-get -qq install -y nvidia-container-toolkit
  sudo nvidia-ctk runtime configure --runtime=docker && sudo systemctl restart docker
fi
docker --version

mkdir -p "$BASE"; cd "$BASE"

# ---------------------------------------------------------------------------
step "[2/6] RFdiffusion @ $RFDIFFUSION_SHA"
if [ ! -d RFdiffusion ]; then
  git clone -q https://github.com/RosettaCommons/RFdiffusion.git
fi
git -C RFdiffusion fetch -q --all
git -C RFdiffusion checkout -q "$RFDIFFUSION_SHA"
cd RFdiffusion
mkdir -p models
if [ ! -f models/Base_ckpt.pt ]; then
  bash scripts/download_models.sh models/ 2>&1 | tail -3
fi
echo "  $(ls models/*.pt | wc -l | tr -d " ") checkpoints"

# ---------------------------------------------------------------------------
step "[3/6] docker image + numpy patch"
# Build the stock image, then layer the numpy pin on top. Keeping the patch as
# a separate one-line Dockerfile rather than editing the upstream one means a
# future RFdiffusion bump does not silently drop it.
docker image inspect rfdiffusion:latest >/dev/null 2>&1 \
  || docker build -q -f docker/Dockerfile -t rfdiffusion . | tail -2

mkdir -p "$BASE/patch"
cat > "$BASE/patch/Dockerfile" <<'PATCH'
FROM rfdiffusion
# torch 1.12.1 is built against the numpy 1.x C ABI; the image ships 2.0.2,
# which fails at runtime with 'Numpy is not available'.
RUN python3.9 -m pip install --no-cache-dir 'numpy<2'
PATCH
docker build -q -f "$BASE/patch/Dockerfile" -t rfdiffusion:fixed "$BASE/patch" | tail -1
docker run --rm --gpus all --entrypoint nvidia-smi rfdiffusion:fixed \
  --query-gpu=name --format=csv,noheader

# ---------------------------------------------------------------------------
step "[4/6] dl_binder_design @ $DLBINDER_SHA  +  ProteinMPNN @ $PROTEINMPNN_SHA"
cd "$BASE"
if [ ! -d dl_binder_design ]; then
  git clone -q https://github.com/nrbennet/dl_binder_design.git
fi
git -C dl_binder_design fetch -q --all
git -C dl_binder_design checkout -q "$DLBINDER_SHA"

# The upstream submodule comes up empty on a plain clone, so clone it directly.
MPNN_DIR="$BASE/dl_binder_design/mpnn_fr/ProteinMPNN"
if [ -z "$( { ls -A "$MPNN_DIR" 2>/dev/null)" ]; then
  rm -rf "$MPNN_DIR"
  git clone -q https://github.com/dauparas/ProteinMPNN.git "$MPNN_DIR"
fi
git -C "$MPNN_DIR" checkout -q "$PROTEINMPNN_SHA"
echo "  $(ls "$MPNN_DIR"/vanilla_model_weights/*.pt 2>/dev/null || true; } | wc -l | tr -d " ") MPNN checkpoints"

# ---------------------------------------------------------------------------
step "[5/6] AF2 params (~5 GB)"
PARAMS="$BASE/dl_binder_design/af2_initial_guess/model_weights/params"
mkdir -p "$PARAMS"
if [ ! -f "$PARAMS/params_model_1_ptm.npz" ]; then
  ( cd "$PARAMS" \
    && wget -q --show-progress "https://storage.googleapis.com/alphafold/$AF2_PARAMS" \
    && tar --extract --file="$AF2_PARAMS" \
    && rm -f "$AF2_PARAMS" )
fi
echo "  $(ls "$PARAMS"/*.npz | wc -l | tr -d " ") parameter files"

# ---------------------------------------------------------------------------
step "[6/6] micromamba $MICROMAMBA_VERSION + environments"
if [ ! -x "$MM" ]; then
  ( cd "$HOME" && curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
      | tar -xj bin/micromamba )
fi
echo "  $($MM --version)"

# The yml files ship with dl_binder_design and name the envs 'mpnn' and 'af2'.
# For a byte-identical rebuild use the captured lock files instead -- see the
# note at the bottom of this script.
for y in proteinmpnn_fastrelax af2_binder_design; do
  $MM env export -n "$( [ "$y" = proteinmpnn_fastrelax ] && echo mpnn || echo af2 )" \
     >/dev/null 2>&1 && { echo "  env for $y already present"; continue; }
  echo "  creating env from include/${y}.yml ..."
  $MM create -y -q -f "$BASE/dl_binder_design/include/${y}.yml" 2>&1 | tail -2
done

verify

cat <<'NOTE'

--------------------------------------------------------------------------
Exact rebuild
  The resolved package sets from the working box are captured alongside this
  script as env/mpnn.lock.yml and env/af2.lock.yml. The steps above rebuild
  from the upstream ymls, which re-solve and may drift. To pin exactly:

      micromamba create -y -n mpnn -f env/mpnn.lock.yml
      micromamba create -y -n af2  -f env/af2.lock.yml

Running the pipeline
  Targets go in $HOME/pmhc/targets and are mounted at /app/targets inside the
  container; outputs land in $HOME/pmhc/out, mounted at /app/out. The run
  scripts (run_rfdiffusion*.sh, run_arm.sh, run_scaffold_arm*.sh) assume that
  layout. Model weights mount from $HOME/pmhc/RFdiffusion/models.

Known failure modes
  "Numpy is not available"        numpy >= 2 leaked into the container; the
                                  rfdiffusion:fixed layer is the fix.
  SSH dies under load             the design loop will starve sshd on a 2-vCPU
                                  box. Protect it before long runs:
                                  sudo systemctl set-property --runtime \
                                    sshd.service MemoryMin=256M IOWeight=1000 \
                                    CPUWeight=1000
  Empty ProteinMPNN directory     the upstream submodule does not populate on
                                  a plain clone; step 4 clones it directly.
--------------------------------------------------------------------------
NOTE
