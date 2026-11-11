#!/usr/bin/env bash
# Unattended overnight orchestration: four arms, in priority order.
#
# Waits for free-design diffusion to finish, then runs each arm through
# MPNN + AF2 + metrics and stops the box.
#
#   1  KH46 control   plumbing check -- can this path recover a known-good pose?
#   2  transfer        does a validated architecture transfer to a new peptide?
#   3  liu scaffolds   partial diffusion from the 152 PRIVILEGED SCAFFOLDS
#                      released with the published work -- skips the
#                      expensive backbone-discovery step entirely
#   4  fold            the published protocol's STRUCTURE -- select one fold
#                      of 50, condition it with insertion, 100 backbones --
#                      but with an automated proxy where they used an expert
#                      eye. Not a reproduction; see run_rfdiffusion_scaffolded
#   5  free            baseline, all 400 backbones
#
# ORDER IS BY WHAT SHOULD SURVIVE A TRUNCATION
#
# AF2 writes its scorefile incrementally and run_arm.sh is checkpointed, so
# an overrun truncates the LAST arm rather than losing everything. The free
# arm is therefore last, and that is a change: it used to be third, with the
# fold arm last on the grounds that conditioning was the narrowest question.
# Two things moved. The fold arm now reproduces the published protocol
# properly, which makes it the most informative arm here; and the free arm's
# second 200 backbones are the most expendable work in the run. A truncated
# free arm still yields a valid funnel at whatever n it reached -- no other
# arm degrades that gracefully, since each of the others is a single
# condition that is either measured or not.
#
# The control runs FIRST because it gates interpretation of the fold arm: if
# a known-good fold handed to the scaffold-guided path yields nothing, a null
# there is about our plumbing, not about conditioning. It already passed --
# 2/10 backbones, best design ipAE 5.85 and peptide_focus 0.357 against
# KH46's own 0.328 -- so that gate is open.
#
# It deliberately does NOT auto-skip anything on a control failure. Gating
# automatically on a threshold that has never been validated would be a
# worse error than reading the log in the morning.
#
# FAILURE ISOLATION
#
# No `set -e` around the stages: each arm must be able to fail without
# taking the others down.
set -uo pipefail

BASE="$HOME/pmhc"
LOG="$BASE/chain.log"
exec > >(tee -a "$LOG") 2>&1

MM="$HOME/bin/micromamba"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"

echo "######################################################################"
echo "# four-arm chain armed at $(date -u)"
echo "######################################################################"

step() { echo; echo "======================================================"; \
         echo ">>> $* :: $(date -u)"; echo "======================================================"; }

# Resume: an arm with a populated metrics file is finished and is not redone.
# The chain has already had to be restarted mid-run once; without this, a
# restart would re-diffuse and re-score work that was already paid for.
done_already() {
  local arm=$1
  local f="$BASE/${arm}_metrics.tsv"
  local n
  [ -f "$f" ] || return 1
  n=$(( $(wc -l < "$f") - 1 ))
  [ "$n" -gt 0 ] || return 1
  echo "arm $arm already has $n scored designs; skipping it"
  return 0
}

summarise() {   # arm name -> funnel, so the morning read is one grep
  # Each `local` gets its own line. `local arm=$1 f="${arm}_x"` looks fine and
  # is not: bash expands every word of the `local` command before running it,
  # so ${arm} resolves against the OUTER scope, where it does not exist. Under
  # `set -u` that is fatal to the whole script, not just the function -- it
  # killed this chain after arm 1 once already, with the box then idle until
  # the dead-man. The arms are isolated from each other's failures; they were
  # not isolated from a bug in this helper.
  local arm=$1
  local f="$BASE/${arm}_metrics.tsv"
  if [ -f "$f" ]; then
    echo "--- FUNNEL: $arm ---"
    "$MM" run -n mpnn python "$BASE/scripts/compare_arms.py" "$f" 2>&1 | sed 's/^/    /'
  else
    echo "--- FUNNEL: $arm -- NO METRICS ---"
  fi
}

# --- wait for free diffusion ------------------------------------------------
step "waiting for free diffusion (tmux session 'pilot')"
while tmux has-session -t pilot 2>/dev/null; do sleep 60; done
echo "free diffusion finished; backbones: $(ls "$BASE"/out/pilot_*.pdb 2>/dev/null | wc -l)"

# === ARM 1: KH46 blueprint control =========================================
step "ARM 1/5  KH46 blueprint control (10 designs, 6 seqs)"
if done_already kh46ctl; then
  summarise kh46ctl
elif [ -f "$BASE/kh46_fold/kh46_ss.pt" ]; then
  bash "$BASE/run_rfdiffusion_scaffolded.sh" 10 kh46ctl \
       "$BASE/kh46_fold" "$BASE/out_kh46ctl" || echo "KH46 CONTROL DIFFUSION FAILED"
  bash "$BASE/run_arm.sh" kh46ctl "$BASE/out_kh46ctl" 6 0 'kh46ctl_*.pdb' \
    || echo "KH46 CONTROL ARM FAILED"
  summarise kh46ctl
  echo
  echo "### READ THIS FIRST: if the control above passes nothing, the"
  echo "### scaffold-guided path is suspect and arm 4's result means little."
else
  echo "no KH46 blueprint at $BASE/kh46_fold; skipping the control"
fi

# === ARM 2: transfer, gated on a prepared target ===========================
# Gated rather than assumed: the new target needs trimming, superposition and
# a partial_T choice, and a hastily prepared target would spend 12 GPU-hours
# producing nonsense. The marker is written only after a smoke test passes.
#
# PARTIAL_T may hold several whitespace-separated values, and does: the budget
# is a fixed number of designs either way, so splitting it across a sweep costs
# nothing and answers "how much remodelling does a transfer need?" instead of
# assuming one value. Each value writes its own prefix into one directory, so a
# single run_arm.sh pass scores them all and the design names carry the value
# for free at analysis time.
step "ARM 2/5  partial-diffusion transfer (~100 designs, 6 seqs)"
TRANSFER_INPUT="$BASE/transfer/complex.pdb"
if done_already transfer; then
  summarise transfer
elif [ -f "$BASE/transfer/READY" ] && [ -f "$TRANSFER_INPUT" ]; then
  PT_LIST=$(cat "$BASE/transfer/PARTIAL_T" 2>/dev/null || echo 20)
  # shellcheck disable=SC2086
  set -- $PT_LIST
  PER=$(( 100 / $# ))

  # SMOKE GATE: one design, end to end, before committing 12 GPU-hours.
  #
  # This runs here rather than by hand because the GPU was busy with free
  # diffusion when the target was prepared, and a concurrent run risked an
  # OOM that would have killed 356 finished backbones. In the chain it costs
  # ~2 minutes, uses the real hardware in the real sequence, and a failure
  # skips the sweep instead of filling a directory with nonsense.
  echo "--- smoke test: 1 design at partial_T=$1 ---"
  SMOKE="$BASE/out_smoke"
  rm -rf "$SMOKE"
  bash "$BASE/run_partial_transfer.sh" 1 smoke "$TRANSFER_INPUT" "$1" "$SMOKE" \
    || echo "SMOKE DIFFUSION FAILED"
  SMOKE_PDB=$(find "$SMOKE" -maxdepth 1 -name 'smoke_*.pdb' -type f | head -1)
  SMOKE_OK=0
  if [ -n "$SMOKE_PDB" ]; then
    # The converter is the real test: it resolves binder/MHC/peptide by exact
    # length, so it fails loudly if diffusion emitted a layout we did not
    # expect -- which is the failure mode that would otherwise be silent.
    # Status captured directly rather than through a pipe into sed: this is a
    # gate, and `cmd | sed` would report sed's exit status without pipefail.
    CONV_LOG=$("$MM" run -n af2 python "$BASE/to_two_chain.py" \
                 "$SMOKE_PDB" "$SMOKE/2chain" 2>&1)
    CONV_RC=$?
    echo "$CONV_LOG" | sed 's/^/      /'
    if [ "$CONV_RC" -eq 0 ]; then
      SMOKE_OK=1
    else
      echo "SMOKE CONVERSION FAILED (rc=$CONV_RC) -- unexpected chain layout"
    fi
  else
    echo "SMOKE TEST produced no backbone"
  fi

  if [ "$SMOKE_OK" -ne 1 ]; then
    echo "### TRANSFER SMOKE TEST FAILED -- skipping the sweep."
    echo "### The other three arms are unaffected."
  else
    echo "smoke test passed; sweeping partial_T over: $PT_LIST ($PER each)"
    for PT in $PT_LIST; do
      bash "$BASE/run_partial_transfer.sh" "$PER" "transferpt${PT}" \
           "$TRANSFER_INPUT" "$PT" "$BASE/out_transfer" \
        || echo "TRANSFER DIFFUSION FAILED at partial_T=$PT"
    done
    bash "$BASE/run_arm.sh" transfer "$BASE/out_transfer" 6 0 \
         'transferpt*.pdb' || echo "TRANSFER ARM FAILED"
    summarise transfer
  fi
else
  echo "transfer target not prepared (no $BASE/transfer/READY); skipping."
  echo "This is the designed fallback, not a failure -- the arm is queued"
  echo "and will run once the target is prepared and smoke-tested."
fi

# === ARM 3: Liu et al.'s privileged scaffolds, by partial diffusion =======
#
# PRIORITISED over both the fold arm and the free arm. The published work's
# own stated reason for building a scaffold set is that de novo backbone
# discovery is expensive and low-yield; those 152 scaffolds ARE that step's
# output, and they are released. Starting from them skips it, which is the
# only version of this experiment that fits a single GPU, so it runs before
# the two arms that re-derive backbones from scratch.
#
# scripts/dock_scaffolds.py has already superimposed each scaffold's own pMHC
# onto our target and carried the binder across -- the same operation as the
# published align_chainB.py, plus the checks it omits: per-scaffold
# superposition RMSD (median 0.43-0.77 A, after correcting a one-residue
# N-terminal register offset that had been rejecting 80 of 152), backbone
# clashes separated from side-chain clashes, and every scaffold's source
# allele recorded in a manifest so the funnel can be split by it.
#
# Static peptide_focus of the placed complexes: median 0.354, 138 of 140
# already above the 0.236 cutoff, against KH46's own crystal at 0.328. That
# is why partial_T is 10 here but swept 10/20/30 in the KH46 transfer arm:
# these start ABOVE the filter, so the pose is the asset and the job is to
# perturb it as little as possible.
step "ARM 3/5  Liu privileged scaffolds, partial diffusion (142, 6 seqs)"
if done_already liu; then
  summarise liu
else
  LIU_DOCKED="$BASE/docked_liu"
  if [ -d "$LIU_DOCKED" ] && \
     [ "$(ls "$LIU_DOCKED"/scaf*.pdb 2>/dev/null | wc -l)" -gt 0 ]
  then
    bash "$BASE/run_scaffold_arm.sh" "$LIU_DOCKED" "$BASE/out_liu" 10 \
      || echo "LIU SCAFFOLD DIFFUSION FAILED"
    if [ "$(ls "$BASE"/out_liu/liu*_*.pdb 2>/dev/null | wc -l)" -gt 0 ]; then
      bash "$BASE/run_arm.sh" liu "$BASE/out_liu" 6 0 'liu*_*.pdb' \
        || echo "LIU ARM FAILED"
      summarise liu
    else
      echo "skipping liu arm: no backbones"
    fi
  else
    echo "no docked scaffolds at $LIU_DOCKED; skipping arm 3."
  fi
fi

# === ARM 4: fold conditioning -- the published protocol, reproduced ========
#
# SKIPPED IN THE 2026-10-31 RUN, by blocking make_scaffolds.sh on the box
# (original kept as make_scaffolds.sh.orig). The arm's own question had by
# then been answered twice: the topology survey found free design already
# produces the published 4-helix topology class unprompted (0/213 strand), so
# conditioning cannot supply a missing fold; and arm 3 showed borrowed
# scaffolds raise the per-sequence pAE rate 4.2x while leaving the fraction of
# viable BACKBONES unchanged (20% vs 23%, p=1.0), so conditioning would not
# buy the pose either. What was given up: the only arm varying binder LENGTH
# (sampled_insertion 0-30), and the only structural reproduction of the
# published protocol. The ~5-8 h went to arm 5, the unconditioned baseline
# that every other arm is measured against and that had not yet run.
# tests/dryrun_chain.sh nofold covers the skip path.
#
# RUNS BEFORE THE FREE ARM, deliberately. It was last when it was the
# narrowest question; it is now the closest thing here to a reproduction of
# the published result, and the free arm's second 200 backbones are the most
# expendable work in the run. Whatever truncates should be those, not this --
# and a truncated free arm still yields a valid funnel at whatever n it
# reached, because AF2 writes its scorefile incrementally. Nothing else here
# degrades that gracefully.
#
# WHAT THE PUBLISHED PROTOCOL ACTUALLY WAS
#   50 unconditional monomers -> SELECT ONE -> loop-masked conditioning with
#   insertion -> 100 backbones (99-129 aa from a 99-residue fold) -> 6 seqs.
#
# Two things were missing here and both are now supplied:
#
#   selection    Conditioning across all 50 folds gives any ONE fold ~2 of
#                100 backbones -- a fiftieth of the exposure -- so a null
#                would have measured "is a RANDOM fold any good", a weaker
#                question. select_scaffold.py substitutes an explicit
#                criterion for an unreproducible visual judgement.
#
#   insertion    scaffoldguided defaults leave sampled_insertion at 0, so
#                every design comes out at the scaffold's exact length;
#                measured on arm 1, all ten backbones are exactly 128. The
#                published 99 -> 99-129 aa spread IS insertion. 0-30 also
#                centres this arm near 130, the free arm's fixed length,
#                removing a length confound rather than adding one.
#
# If no scaffold meets the criterion the arm falls back to the whole library
# rather than dying, so the worst case is the behaviour we had before.
step "ARM 4/5  fold conditioning, auto-selected fold (100 designs, 6 seqs)"
if done_already fold; then
  summarise fold
else
  FOLD_OK=1
  bash "$BASE/make_scaffolds.sh" 50 100-130 \
    || { echo "SCAFFOLDS FAILED"; FOLD_OK=0; }

  if [ "$FOLD_OK" -eq 1 ]; then
    SEL="$BASE/scaffold_selected"
    SCAF_DIR="$BASE/scaffold_folds"      # fallback: the whole library
    if "$MM" run -n mpnn python "$BASE/scripts/select_scaffold.py" \
         "$BASE/scaffold_pdbs" "$SCAF_DIR" "$SEL" 1
    then
      SCAF_DIR="$SEL"
    else
      echo "### no scaffold met the selection criterion."
      echo "### falling back to library-wide conditioning for this arm."
    fi
    bash "$BASE/run_rfdiffusion_scaffolded.sh" 100 fold "$SCAF_DIR" \
         "$BASE/out_scaffolded" 30 \
      || echo "FOLD DIFFUSION FAILED"
  fi

  if [ "$(ls "$BASE"/out_scaffolded/fold_*.pdb 2>/dev/null | wc -l)" -gt 0 ]
  then
    bash "$BASE/run_arm.sh" fold "$BASE/out_scaffolded" 6 0 'fold_*.pdb' \
      || echo "FOLD ARM FAILED"
    summarise fold
  else
    echo "skipping fold arm: no backbones"
  fi
fi

# === ARM 1b: KH46 onto its own target, by partial diffusion ===============
#
# The missing cell in the 2x2. Arm 1 conditions on KH46's ss/adj; arms 2 and 3
# partially diffuse from a docked pose. So arm 1 differs from them in BOTH
# method and starting material, and isolates neither. This arm is
# method-matched to both:
#
#   vs arm 2   same source, same method, same partial_T sweep; only the
#              PEPTIDE changes -> isolates the target, and removes the
#              "arm 2 starts from a validated binder for this allele"
#              confound that currently qualifies arm 2's whole result
#   vs arm 3   same target, same method; only the SOURCE changes ->
#              isolates source quality, validated-for-this-peptide against
#              validated-for-other-peptides
#
# It is a POSITIVE CONTROL: a ceiling and a sanity check, not a discovery
# test. At partial_T=10 it largely regenerates a binder that already works, so
# a high pass rate here is the expected result and is evidence of nothing. Its
# value is as the denominator arms 2 and 3 currently lack.
#
# Mirrors arm 2 exactly -- 10/20/30, 33 each, 6 seqs -- so every comparison is
# PAIRED at each partial_T rather than pooled across them.
#
# RUNS BEFORE THE FREE ARM: it makes two results we already have
# interpretable, and the free arm is the only one that still yields a valid
# funnel when truncated.
step "ARM 1b  KH46 onto its own target, partial diffusion (99, 6 seqs)"
if done_already kh46self; then
  summarise kh46self
else
  SELF_INPUT="$BASE/kh46self/complex.pdb"
  if [ -f "$SELF_INPUT" ]; then
    for PT in 10 20 30; do
      bash "$BASE/run_partial_transfer.sh" 33 "kh46selfpt${PT}" \
           "$SELF_INPUT" "$PT" "$BASE/out_kh46self" \
        || echo "KH46SELF DIFFUSION FAILED at partial_T=$PT"
    done
    if [ "$(ls "$BASE"/out_kh46self/kh46selfpt*.pdb 2>/dev/null | wc -l)" \
         -gt 0 ]
    then
      bash "$BASE/run_arm.sh" kh46self "$BASE/out_kh46self" 6 0 \
           'kh46selfpt*.pdb' || echo "KH46SELF ARM FAILED"
      summarise kh46self
    else
      echo "skipping kh46self arm: no backbones"
    fi
  else
    echo "no self-control input at $SELF_INPUT; skipping arm 1b."
    echo "This is the designed fallback, not a failure. Build it with:"
    echo "  python3 scripts/prep_self_control.py \\"
    echo "      targets/9min_control_kh46.pdb kh46self/complex.pdb"
  fi
fi

# === ARM 5: free design, the FLOOR =========================================
#
# HALF the 400 available backbones. This arm's only job is to be the
# denominator: the unconditioned rate that arms 3 and 1b are measured against.
# It does not need to resolve its own rate precisely, only to sit clearly
# below arm 3, so the question is power, not sample size for its own sake.
#
# Arm 3 is 33/141 backbones with a hit (23.4%). Fisher, two-sided, against a
# free arm at various true rates:
#
#     free rate      n=200            n=400
#        5%          p=0.0000         p=0.0000
#       10%          p=0.0013         p=0.0002
#       13%          p=0.0138         p=0.0047
#       15%          p=0.0649         p=0.0273
#
# So 200 separates arm 3 from any free rate up to ~13%, and 400 only extends
# that to ~15%. Arm 1 -- fold-conditioned on the published winner's own fold,
# the most guided arm here -- reached 20%, so an UNCONDITIONED arm above 13%
# would be a surprising result in its own right. The second 200 backbones buy
# resolution only inside a 13-15% band where neither n settles the question
# cleanly, and cost ~4 h of AF2.
#
# The backbones already exist, so the cap saves prediction time only, no
# diffusion. run_arm.sh takes the first 200 by version sort, i.e. pilot_0 to
# pilot_199 in generation order; diffusion designs are i.i.d. so that is an
# unbiased subsample. The remaining 200 can be added later without rework if
# the margin allows.
#
# Last in the queue because it is the only arm that still yields a valid
# funnel at whatever n it reaches, since AF2 writes its scorefile
# incrementally.
step "ARM 5/5  free design, the floor (200 of 400 backbones, 4 seqs)"
if done_already free; then
  summarise free
else
  bash "$BASE/run_arm.sh" free "$BASE/out" 4 200 'pilot_*.pdb' \
    || echo "FREE ARM FAILED"
  summarise free
fi

# --- summary ----------------------------------------------------------------
echo
echo "######################################################################"
echo "# chain done at $(date -u)"
for arm in kh46ctl kh46self transfer liu fold free; do
  f="$BASE/${arm}_metrics.tsv"
  if [ -f "$f" ]; then
    echo "#   $arm: $(( $(wc -l < "$f") - 1 )) designs scored"
  else
    echo "#   $arm: not run"
  fi
done
echo "######################################################################"

# Stop, not terminate: compute billing ends, EBS and all results survive.
echo "stopping instance in 5 minutes"
sudo shutdown -h +5
