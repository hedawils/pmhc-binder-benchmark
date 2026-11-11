#!/usr/bin/env bash
# Full provisioning for the pMHC binder pilot, from a bare Deep Learning AMI.
#
#   AMI: ami-012ba162b9cd2729c
#        Deep Learning OSS Nvidia Driver AMI GPU PyTorch 2.7 (Ubuntu 22.04)
#   Instance: g5.2xlarge (A10G 24GB) or g5.xlarge (same GPU, 16GB RAM)
#
# ---------------------------------------------------------------------------
# WHY THIS FILE EXISTS
#
# Standing up RFdiffusion + ProteinMPNN + AF2 in late 2026 takes ten separate
# version fixes, none of them documented in the upstream repos. Six are
# downstream of the numpy 2.x migration: these pipelines were pinned before
# it landed, and pip/conda now resolve to versions that break them at import
# or at runtime. Without this file the working box is a snowflake, and the
# knowledge dies with the instance.
#
#   1  RFdiffusion resolves checkpoints relative to its pip install location
#      (dist-packages/models), not to the repo's models/ directory.
#   2  numpy 2.0.2 vs torch 1.12.1  -> "RuntimeError: Numpy is not available"
#   3  pyrosetta is not on conda-forge; install via PyPI pyrosetta-installer
#   4  mkl >=2024 vs torch 2.0      -> "undefined symbol: iJIT_NotifyEvent"
#   5  numpy 2.x vs torch 2.0.1     -> ABI break in the numpy bridge
#   6  numpy 2.x vs jax 0.4.20      -> "numpy has no attribute 'trapz'"
#   7  cuDNN 9.x vs jaxlib cudnn89  -> JAX silently falls back to CPU
#   8  biopython >=1.80 dropped Bio.Data.SCOPData, which AlphaFold imports
#   9  dm-haiku 0.0.17 vs jax 0.4.20 -> "jax.core has no take_current_trace"
#  10  tensorflow 2.21 vs ml-dtypes 0.3.1 required by jax 0.4.20
#
# Note #7 in particular: JAX does not fail loudly, it just runs on CPU, which
# would turn a 22-hour AF2 run into a multi-week one. Always assert on
# jax.devices() rather than trusting the install.
#
# One consequence worth remembering: the metrics collector runs in the MPNN
# env, because af2 is pinned to biopython 1.79 (for SCOPData) and 1.79 calls
# np.int, which numpy >=1.24 removed. The two constraints are unsatisfiable
# in a single environment.
# ---------------------------------------------------------------------------
set -euo pipefail

BASE="$HOME/pmhc"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
MM="$HOME/bin/micromamba"
mkdir -p "$BASE"

echo "########## [1/5] RFdiffusion ##########"
cd "$BASE"
[ -d RFdiffusion ] || git clone -q https://github.com/RosettaCommons/RFdiffusion.git
cd RFdiffusion
mkdir -p models
[ -f models/Complex_base_ckpt.pt ] || bash scripts/download_models.sh models/
docker build -q -f docker/Dockerfile -t rfdiffusion . >/dev/null

# Fix 2: the official image ships numpy 2.x against torch 1.12.1.
mkdir -p "$BASE/patch" && cd "$BASE/patch"
cat > Dockerfile <<'EOF'
FROM rfdiffusion
RUN python3.9 -m pip install --no-cache-dir 'numpy<2'
EOF
docker build -q -t rfdiffusion:fixed . >/dev/null
docker run --rm --entrypoint python3.9 rfdiffusion:fixed -c \
  "import numpy,torch; assert numpy.__version__.startswith('1.'); \
   torch.from_numpy(numpy.zeros(3)); print('  rfdiffusion:fixed ok')"

echo "########## [2/5] micromamba ##########"
cd "$HOME"
[ -x "$MM" ] || curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba

echo "########## [3/5] dl_binder_design ##########"
cd "$BASE"
[ -d dl_binder_design ] || git clone -q https://github.com/nrbennet/dl_binder_design.git
cd dl_binder_design
# The ProteinMPNN submodule does not populate; clone it directly.
if [ -z "$(ls -A mpnn_fr/ProteinMPNN 2>/dev/null)" ]; then
  rm -rf mpnn_fr/ProteinMPNN
  git clone -q https://github.com/dauparas/ProteinMPNN.git mpnn_fr/ProteinMPNN
fi
PARAMS="$PWD/af2_initial_guess/model_weights/params"
mkdir -p "$PARAMS"
if [ ! -f "$PARAMS/params_model_1_ptm.npz" ]; then
  ( cd "$PARAMS"
    wget -q https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar
    tar --extract --file=alphafold_params_2022-12-06.tar
    rm -f alphafold_params_2022-12-06.tar )
fi

echo "########## [4/5] mpnn env ##########"
# Fix 3: strip pyrosetta from the yml, install it from PyPI afterwards.
grep -v pyrosetta include/proteinmpnn_fastrelax.yml > /tmp/mpnn.yml
"$MM" create -y -q -n mpnn -f /tmp/mpnn.yml
"$MM" install -y -q -n mpnn 'mkl<2024' 'numpy<2'          # fixes 4, 5
"$MM" run -n mpnn pip install -q pyrosetta-installer
"$MM" run -n mpnn python -c \
  "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta(serialization=True)"
"$MM" run -n mpnn python -c \
  "import torch,numpy,pyrosetta; assert torch.cuda.is_available(); \
   torch.from_numpy(numpy.zeros(3)); print('  mpnn env ok')"

echo "########## [5/5] af2 env ##########"
grep -v pyrosetta include/af2_binder_design.yml > /tmp/af2.yml
"$MM" create -y -q -n af2 -f /tmp/af2.yml
# Fix 7: jaxlib 0.4.20 is built against cuDNN 8.9 and CUDA 12.2. The yml's
# constraints are all >=, so pip resolves to a 12.9-era stack that fails.
"$MM" run -n af2 pip install -q \
  'nvidia-cublas-cu12==12.2.5.6' 'nvidia-cuda-cupti-cu12==12.2.142' \
  'nvidia-cuda-nvcc-cu12==12.2.140' 'nvidia-cuda-runtime-cu12==12.2.140' \
  'nvidia-cudnn-cu12==8.9.7.29' 'nvidia-cufft-cu12==11.0.8.103' \
  'nvidia-cusolver-cu12==11.5.2.141' 'nvidia-cusparse-cu12==12.1.2.141' \
  'nvidia-nccl-cu12==2.19.3'
# Fixes 8, 9, 10 plus the numpy pin they keep undoing.
"$MM" run -n af2 pip install -q --force-reinstall --no-deps 'tensorflow-cpu==2.16.1'
"$MM" run -n af2 pip install -q 'biopython==1.79' 'dm-haiku==0.0.10' \
  'chex==0.1.85' 'ml-dtypes==0.3.1' 'numpy<2' 'protobuf<5'
"$MM" run -n af2 pip install -q pyrosetta-installer
"$MM" run -n af2 python -c \
  "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta(serialization=True)"

# Fix 7 is silent on failure, so assert rather than hope.
"$MM" run -n af2 python -c "
import jax, numpy, tensorflow as tf
from Bio.Data import SCOPData
devs = jax.devices()
assert any(d.platform == 'gpu' for d in devs), f'JAX ON CPU: {devs}'
assert numpy.__version__.startswith('1.'), numpy.__version__
print('  af2 env ok -', devs)
"

echo
echo "SETUP COMPLETE"
