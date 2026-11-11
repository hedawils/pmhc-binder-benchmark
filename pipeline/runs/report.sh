#!/usr/bin/env bash
# The morning read: every arm's funnel, side by side, plus what went wrong.
#
#   usage:  report.sh            on the box
#
# chain.log runs to thousands of lines of RFdiffusion timestep spam, so the
# answers are not findable by scrolling. This pulls out the three things worth
# knowing: which arms produced data, how their funnels compare, and whether
# anything failed or fell back.
set -uo pipefail

BASE="$HOME/pmhc"
MM="$HOME/bin/micromamba"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"
LOG="$BASE/chain.log"

echo "######################################################################"
echo "# pMHC pilot -- run report, $(date -u)"
echo "######################################################################"

echo
echo "=== arms ==="
echo "  NOTE: 'transfer' is the only arm NOT on NY-ESO-1. Its funnel is a"
echo "        different TARGET (MART-1 AAGIGILTV on the same allele), so its"
echo "        column in the side-by-side below is not directly comparable."
printf "  %-10s %-12s %-10s %s\n" ARM TARGET BACKBONES SCORED
for pair in "kh46ctl:NY-ESO-1:out_kh46ctl/kh46ctl_*.pdb" \
            "transfer:MART-1:out_transfer/transferpt*.pdb" \
            "liu:NY-ESO-1:out_liu/liu*_*.pdb" \
            "fold:NY-ESO-1:out_scaffolded/fold_*.pdb" \
            "free:NY-ESO-1:out/pilot_*.pdb"; do
  arm=${pair%%:*}; rest=${pair#*:}; tgt=${rest%%:*}; glob=${rest#*:}
  # shellcheck disable=SC2086
  bb=$(ls $BASE/$glob 2>/dev/null | wc -l | tr -d ' ')
  if [ -f "$BASE/${arm}_metrics.tsv" ]; then
    scored=$(( $(wc -l < "$BASE/${arm}_metrics.tsv") - 1 ))
  else
    scored="not run"
  fi
  printf "  %-10s %-12s %-10s %s\n" "$arm" "$tgt" "$bb" "$scored"
done

echo
echo "=== binder length distribution (did insertion do anything?) ==="
for pair in "kh46ctl:out_kh46ctl/kh46ctl_*.pdb" "fold:out_scaffolded/fold_*.pdb"; do
  arm=${pair%%:*}; glob=${pair#*:}
  # shellcheck disable=SC2086
  files=$(ls $BASE/$glob 2>/dev/null | head -60)
  [ -n "$files" ] || { printf "  %-10s no backbones\n" "$arm"; continue; }
  printf "  %-10s " "$arm"
  for f in $files; do
    awk '$1=="ATOM" && substr($0,13,4)~/CA/ && substr($0,22,1)=="A"' "$f" | wc -l
  done | sort -n | uniq -c | awk '{printf "%s res x%s  ", $2, $1}'
  echo
done

echo
echo "=== failures, skips and fallbacks ==="
grep -nE "FAILED|unbound variable|skipping|already has|falling back|SMOKE|no scaffold met|not prepared" \
     "$LOG" 2>/dev/null | grep -v "^.*Timestep" | tail -25 | sed 's/^/  /' \
  || echo "  none found"

echo
echo "=== selection criterion outcome ==="
grep -A 14 "surveyed .* scaffolds" "$LOG" 2>/dev/null | tail -18 | sed 's/^/  /' \
  || echo "  selection did not run"

echo
echo "=== arm 3 (liu): accepted scaffolds by SOURCE allele ==="
if [ -f "$BASE/docked_liu/manifest.tsv" ]; then
  awk -F'\t' 'NR>1 && $2==1 {n[$4]++} END {for (a in n) printf "  %5d  %s\n", n[a], a}' \
      "$BASE/docked_liu/manifest.tsv" | sort -rn
  echo "  (join liu_metrics.tsv to this manifest on scaffold to ask whether"
  echo "   a scaffold's ORIGINAL allele predicts whether it transfers)"
else
  echo "  no manifest"
fi

echo
echo "=== funnels, side by side ==="
TSVS=()
for arm in liu free fold kh46ctl transfer; do
  [ -f "$BASE/${arm}_metrics.tsv" ] && TSVS+=("$BASE/${arm}_metrics.tsv")
done
if [ "${#TSVS[@]}" -eq 0 ]; then
  echo "  no metrics tables yet"
else
  "$MM" run -n mpnn python "$BASE/scripts/compare_arms.py" "${TSVS[@]}"
fi

echo
echo "######################################################################"
echo "# chain still running: $(tmux has-session -t chain 2>/dev/null && echo YES || echo no)"
echo "# instance stops itself when the chain finishes; dead-man is the backstop."
echo "######################################################################"
