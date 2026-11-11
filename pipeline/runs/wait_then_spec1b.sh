#!/usr/bin/env bash
# Waits on KF9's PID (not a name match, which would catch this script's own
# pgrep and never terminate), then scores arm 1b for specificity.
set -uo pipefail
W=1606
echo "waiting on KF9 pid $W, started $(date -u '+%F %T UTC')"
while [ "$W" -ne 0 ] && kill -0 "$W" 2>/dev/null; do sleep 60; done
echo "KF9 done at $(date -u '+%F %T UTC'); starting arm 1b specificity"
cd $HOME/pmhc
exec ~/bin/micromamba run -n mpnn python scripts/score_specificity.py \
  kh46self_metrics.tsv kh46self_af2 kh46self_specificity.tsv \
  --samples 10 --decoys decoys_a0201.txt
