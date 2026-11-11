#!/usr/bin/env bash
# Provision RFdiffusion on the pilot box. Runs on the instance, in tmux.
#
# Docker rather than conda deliberately: RFdiffusion's SE3Transformer stack
# (dgl + e3nn + a specific torch/CUDA pairing) is the single most reliable
# way to lose a day to dependency resolution. The official image pins all of
# it. We have one working day, so we do not negotiate with conda.
set -euo pipefail

cd "$HOME"
mkdir -p pmhc && cd pmhc

echo "=== [1/3] clone RFdiffusion ==="
if [ ! -d RFdiffusion ]; then
  git clone -q https://github.com/RosettaCommons/RFdiffusion.git
fi
cd RFdiffusion

echo "=== [2/3] download model weights ==="
mkdir -p models
# The repo ships its own download script; using it avoids hardcoding
# checkpoint URLs that drift.
if [ -f scripts/download_models.sh ]; then
  bash scripts/download_models.sh models/ 2>&1 | tail -5
else
  echo "!! download_models.sh missing -- repo layout changed, check manually"
  exit 1
fi
echo "weights present:"
ls -la models/*.pt 2>/dev/null | awk '{print "   ", $9, $5}' || echo "   none found"

echo "=== [3/3] build docker image ==="
docker build -f docker/Dockerfile -t rfdiffusion . 2>&1 | tail -15

echo
echo "=== verifying GPU visibility inside container ==="
docker run --rm --gpus all rfdiffusion nvidia-smi --query-gpu=name --format=csv,noheader

echo
echo "SETUP COMPLETE"
