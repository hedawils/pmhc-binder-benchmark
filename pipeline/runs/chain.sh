#!/usr/bin/env bash
# Wait for the diffusion run to finish, then start MPNN+AF2 automatically.
# Avoids leaving the GPU idle between stages overnight.
while tmux has-session -t pilot 2>/dev/null; do sleep 60; done
echo "diffusion finished at $(date -u)"
echo "backbones: $(ls $HOME/pmhc/out/pilot_*.pdb 2>/dev/null | wc -l)"
bash $HOME/pmhc/run_pipeline.sh 4
