#!/usr/bin/env bash
# ARM 4 SKIPPED BY REQUEST, 2026-10-31.
#
# Not a failure. The fold-conditioning question it was built to answer has
# already been answered twice over:
#   - the topology survey found 0/213 strand residues and a modal 4-helix
#     bundle, i.e. free design already produces the published topology class,
#     so conditioning cannot supply a fold the unguided run was missing;
#   - arm 3 showed borrowed scaffolds raise the per-sequence pAE rate 4.2x
#     but leave the fraction of viable BACKBONES unchanged (20% vs 23%,
#     p=1.0), so the pose is not what conditioning would buy either.
#
# What is given up: this was the only arm varying binder LENGTH
# (sampled_insertion 0-30) and the only structural reproduction of the
# published protocol. Both are recoverable later; the GPU hours are not.
#
# Skipping frees ~5-8 h, which goes to arm 5 -- the unconditioned 400-backbone
# baseline that is the denominator every other arm is currently missing.
#
# Restore with:  mv make_scaffolds.sh.orig make_scaffolds.sh
echo "    make_scaffolds.sh: ARM 4 SKIPPED BY REQUEST -- not a failure."
echo "    Rationale and restore instructions are in this file's header."
exit 1
