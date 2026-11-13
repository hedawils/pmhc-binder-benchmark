# holo − apo vs the alanine baseline

Cut from the README to keep the specificity section proportionate. Reproduce with
`python analysis/reproduce.py`.

Both baselines share a numerator: `logp_holo`, the log probability of the native peptide residue
with the binder present. They differ in what they subtract.

    delta     = logp_holo - logp_apo    scored again with the binder chain deleted
    delta_ala = logp_holo - logp_ala    the same position mutated to alanine

`delta` cancels everything structural, leaving the binder's contribution. `delta_ala` keeps the
binder in both terms, so the groove's own preference stays in the number.

Per design, mean over contacted positions (P4/P5/P7) minus mean over anchors (P2/P9):

| arm | `delta` | `delta_ala` |
|---|---:|---:|
| repurposed | **+1.27** (p < 1e-4) | −0.01 (p = 0.83) |
| free | **+1.43** (p < 1e-4) | +0.31 (p = 0.023) |
| KH46 | **+0.90** (p < 1e-4) | **−0.25** (p < 1e-4) |

`delta` separates contacted positions from anchors in every arm. `delta_ala` does not, and on the
KH46 arm it inverts — scoring anchors higher than the positions the binder touches. Its largest
value anywhere is P3 at +2.4 ± 0.3 across all 285 designs: large, consistent, and at a position
buried in the groove that no binder reaches. That is a property of HLA-A*02:01, not of any design.

At the level of ranking whole designs the two agree much better, because the groove term is close
to constant across designs on a fixed target. The decoy panel changes design ranking more than the
choice of baseline does.

This compares methods, not results: Liu et al.'s specificity numbers are not published, so the
alanine filter here is my implementation of their one-sentence description.
