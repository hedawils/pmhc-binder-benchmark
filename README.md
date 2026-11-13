# To what extent do privileged scaffolds improve de novo design of pMHC binders?

A controlled benchmark of **partial diffusion from privileged scaffolds** against **free diffusion**,
for de novo binders to peptide–MHC targets. Five arms, two targets, ~3,400 designs, identical
filters throughout.

**Short answer: it probably depends on the target.**

---

## Results

| target | starting point | hit rate | per scaffold |
|---|---|---:|---:|
| NY-ESO-1 `SLLMWITQV` / HLA-A\*02:01 | 142 privileged scaffolds | **7.9%** (67/852) | 19.3% |
| | free diffusion | 1.9% (15/800) | 5.5% |
| | KH46, a binder already proven on this peptide | 34.2% (203/594) | 67.1% |
| SARS-CoV-2 `KAYNVTQAF` / HLA-C\*12:02 | the same 142 scaffolds | **3.6%** (21/588) | 8.2% |
| | free diffusion | 4.2% (25/600) | 11.3% |

Repurposing wins by **3.8× per design** on the first target and by **nothing** on the second
(0.86×, p = 0.65). Per-scaffold rates are placed on a common 4-sequence basis (expectation over all
4-subsets, so no RNG seed is involved).

### The four gates

**A design counts as a hit only if it clears all four.** They are applied identically to every
arm, so no comparison here depends on where the thresholds sit — but the absolute hit rates do,
and the gates are the single biggest assumption in this work.

| gate | threshold | what it asks |
|---|---|---|
| **scRMSD** | < 2.0 Å | does the design refold to the backbone it was designed as? |
| **binder pLDDT** | > 92 | is the binder confidently folded on its own? |
| **interface pAE** | < 7.0 | is the predicted *relative placement* of binder and target confident? |
| **peptide_focus** | ≥ 0.202 | is the binder gripping the peptide, or just the MHC platform? |

The first three are standard practice for AF2-filtered binder design. `peptide_focus` needs
defining:

```
peptide_focus = buried_peptide / (buried_peptide + buried_mhc)
```

the fraction of the target's buried surface area that belongs to the 9-mer peptide rather than the
MHC platform, computed as the change in solvent-accessible surface area of each target residue on
binding. It exists because a pMHC binder that engages only the groove would bind every cell
presenting that allele, regardless of which peptide is in it.

The **0.202** threshold comes from crystal structures of known pMHC binders, where 0.236 separates
genuine peptide engagement from MHC binding — KH46's own crystal scores 0.328 — shifted down
by 0.034 because AF2-predicted complexes score systematically below crystals on this measure.

Two consequences worth stating up front. This gate is a *floor* on the quantity analysed in
[Most of every interface is MHC](#most-of-every-interface-is-mhc-not-peptide) below — the two are
the same metric, not independent findings. And it is the only gate that behaves differently
between targets: inert on NY-ESO-1, where it removes nothing from either arm, but cutting both KF9
arms hard (35→21 and 38→25).

Nothing here was expressed or measured — see [Limitations](#limitations).

### Why it works when it works

Repurposed scaffolds do **not** fold better than designs built from scratch. They fold *worse*.

| all designs, medians | repurposed | free | p |
|---|---:|---:|---:|
| binder pLDDT, NY-ESO-1 | 90.9 | **92.1** | 7.6e-06 |
| binder pLDDT, KF9 | 90.8 | **93.0** | 7.1e-10 |
| designs contacting the peptide, NY-ESO-1 | **43.4%** | 30.0% | 1.6e-08 |
| buried peptide Å², NY-ESO-1 † | **297** | 46 | 9.5e-13 |
| designs contacting the peptide, KF9 | 43.7% | 43.5% | 0.94 |
| buried peptide Å², KF9 † | 107 | 103 | 0.48 |

† among designs that contact the peptide at all; unconditional medians are 0 in every arm.

The pattern suggests the advantage is **interface geometry** rather than foldability — the
borrowed pose already sits where a binder needs to be. On the second target the interface
differences are absent on every measure while the fold penalty remains, which is consistent with
the same reading: the geometry transfers only when it has somewhere to transfer to.

Treat that as the most parsimonious account of the numbers, not as a demonstrated mechanism.
Nothing here isolates geometry as the cause, and the arms differ in binder length as well as in
starting point (130 residues free, median ~85 repurposed), which could contribute to the interface
differences on its own.

### Most of every interface is MHC, not peptide

Share of buried interface contributed by the peptide, in designs that passed all four gates.
This is `peptide_focus` from [the gates above](#the-four-gates), reported as a distribution
rather than as a threshold:

| arm | peptide share | best design |
|---|---:|---:|
| NY-ESO-1 repurposed | 31.9% | 39.9% |
| NY-ESO-1 free | 28.0% | 33.1% |
| KH46 control | 28.8% | 37.9% |
| KF9 repurposed | 22.0% | 31.0% |
| KF9 free | 22.2% | 25.8% |

**~70% of every surviving interface is the MHC.** That is not by itself alarming — natural TCRs
also bury far more MHC than peptide, and are exquisitely peptide-specific — and every design here
clears the `peptide_focus` gate, so all of them engage the peptide to some degree. What the
fraction says is where the binding surface sits, not whether the binder discriminates; that
question belongs to the specificity filter below. The reason to track it is that surface on the
MHC is surface on a polymorphic, universally expressed protein, so the more of the interface it
carries, the more a design's behaviour depends on the allele rather than the antigen.

Across all 550 NY-ESO-1 designs with a real interface the median is 30%, the 95th percentile 47%,
and the maximum **69.8%** — peptide-dominant designs do get generated.

They are then slightly disfavoured by the other three gates. Among designs already clearing the
`peptide_focus` bar of 20.2%, the median peptide share is 33.4%; among those that also clear
scRMSD, pLDDT and interface pAE it falls to **30.7%** (p = 0.0002). KF9 behaves the same way,
24.0% to 22.1% (p < 1e-4). The effect is small, but it runs the wrong way: the confidence filters
mildly prefer designs that bury more MHC, which is the larger surface available.

Note also that hotspots were specified as peptide residues only (`C4,C5,C7,C8` and
`C1,C4,C7,C8`). No MHC residue was ever named, and the interfaces still came out ~70% MHC — in
this pilot, hotspot specification did not control interface composition.

### A specificity filter

Beyond the four gates, each design is scored for whether its interface *prefers* its own peptide.
The approach is **Liu et al.'s** — ProteinMPNN scoring of how predictable each peptide residue is
in complex with the binder. Added here: whole decoy peptides rather than alanine substitution
alone, and a **holo − apo baseline** — the same complex scored with the binder chain deleted, then
subtracted — so the groove's own preference cancels out.

Mean holo − apo per peptide position, NY-ESO-1:

```
position      P1     P2     P3     P4     P5     P6     P7     P8     P9
residue        S      L      L      M      W      I      T      Q      V
repurposed -0.29  +0.01  -0.11  +0.70  +1.69  -0.29  +1.27  -0.41  -0.10
free       -0.38  +0.05  -0.16  +0.92  +2.17  -0.33  +1.10  -0.13  -0.12
KH46       -0.28  +0.01  -0.18  +1.06  +0.51  -0.33  +0.92  -0.41  -0.16
```

The signal lands on **P4, P5 and P7** — the solvent-exposed positions a binder can reach — and is
~0 at **P2 and P9**, the canonical HLA-A\*02:01 anchors buried in the groove. If the metric were
reading the groove rather than the binder, the anchors would dominate. They do not. That is the
closest thing here to a validation of the filter.

Method in `pipeline/scripts/score_specificity.py`; per-design numbers from
`analysis/reproduce.py`, which also compares this baseline against the alanine difference Liu et
al. describe.

### The positive control has the worst specificity

| arm | median margin | fraction preferring on-target |
|---|---:|---:|
| NY-ESO-1 free | +1.15 | 80.0% |
| NY-ESO-1 repurposed | +0.46 | 67.2% |
| KH46 control | **−2.12** | **11.3%** |

Margin is on-target summed holo − apo minus the best decoy's. Repurposed vs free is not
significant (p = 0.15).

The KH46 arm has by far the **highest hit rate** (34.2%) and by far the **worst specificity** —
89% of its designs score some decoy above the peptide they were designed for. It is not a size
effect: KH46 designs (128 aa) and free-diffusion designs (130 aa) are nearly the same length and
sit at opposite ends of this table.

Passing the gates and being selective are different things, and here the arm that optimises the
gates hardest is the least selective.

How far that generalises is unclear: the 203 designs come from **one** source structure at one
noise level, so this is a property of re-diffusing KH46 rather than a demonstrated property of
re-diffusing proven binders. The mechanism is plausible — stay close to a pose selected for
affinity and you inherit its interface, selective or not — but a second source structure is what
would test it.

---

## Limitations

Stated here rather than buried, because they bound what the numbers above can mean.

**n = 2 targets.** Every claim about when repurposing helps rests on one target where it did
and one where it did not. Two points cannot distinguish "depends on the target" from "worked
once", and nothing here says which of NY-ESO-1 or KF9 is the unusual one. The per-design
numbers within each target are well powered; the generalisation across targets is not.

**n = 1 positive control, at one noise level.** The KH46 arm re-diffuses a single published
binder at a single `partial_T`. It is that design's ceiling, not the target's, and its
specificity result — the lowest of any arm — is a property of one source structure. Whether
re-diffusing a proven binder generally produces poorly selective designs cannot be answered
from one of them.

**Nothing was expressed or measured.** Every result is a computational proxy that has never been
checked against binding. The gates are standard practice, not validated ground truth.

**Binder length is confounded with method.** Free-diffusion arms are fixed at 130 residues;
repurposed arms inherit their donor's length, median ~85. Interface-size comparisons between those
arms are therefore confounded. A length-matched free arm at ~85 residues is the single cleanest
follow-up and has not been run.

**The donor cohort is confounded with HLA locus.** Hit rate varies enormously by the target a
scaffold was originally built for — PHOX2B 52.9% (18/34) down to PAP 0% (0/22) on NY-ESO-1, with
free diffusion at 5.5% in between, so some donors are worse than no donor. But all 34 HLA-C donors
come from a **single source target**, so donor identity and locus cannot be separated in this
design. No claim about HLA locus is supportable in either direction.

**The specificity filter has never been checked against binding.** Liu et al.'s equivalent
sits upstream of yeast display and characterisation; this one sits upstream of nothing.

**The decoys are easy.** None of the four shares a residue with the target at the positions the
binders actually read. The filter has not been tested against near-neighbour peptides.

**MPNN relax cycles were disabled** (`-relax_cycles 0`). Applied identically across all arms so it
cannot bias the comparison, but it is not the standard protocol and probably costs absolute hit
rate everywhere.

**The second target differs in both allele and peptide**, so the NY-ESO-1 vs KF9 contrast alone
cannot attribute the failure to either.

**One more asymmetry worth noting:** the fourth gate removes nothing in either NY-ESO-1 arm, but
cuts both KF9 arms hard (35→21 and 38→25). On that target a real share of interface-confident
designs grip the MHC rather than the peptide.

---

## Layout

```
data/           every design, every gate — 15 TSVs, one row per design
analysis/       scripts that turn data/ into the results above
pipeline/       what ran on the GPU box to produce data/
structures/     NOT IN THIS REPO -- see below
targets/        the three target structures
reproduce_env.sh   rebuild the environment from pinned commits
docs/           internal working notes
```

**The 337 surviving structures are not in this repository.** They are ~131 MB, which would
make cloning unpleasant for everyone who only wants the data and the code. They are not
re-creatable either — diffusion is stochastic, so re-running the pipeline produces different
designs, not these. If you want them, open an issue and I will attach them to a release.

`data/manifest_{liu,kf9}.tsv` is the join key for donor provenance — source target, allele,
peptide and transplant-stage geometry for each scaffold.

## Reproducing

```bash
pip install -r requirements.txt
python analysis/reproduce.py        # regenerates every number in this README from data/
```

Rebuilding the design environment itself — RFdiffusion, ProteinMPNN, AF2 initial guess — needs a
GPU box:

```bash
bash reproduce_env.sh          # full build
bash reproduce_env.sh verify   # check an existing install, non-destructive
```

Pinned: RFdiffusion `86507b6`, dl_binder_design `cafa385`, ProteinMPNN `8907e66`,
AF2 params `2022-12-06`. Run on a `g5.2xlarge` (A10G), Ubuntu 22.04.

One non-obvious trap is handled in there: torch 1.12.1 is built against the numpy 1.x C ABI while
the official RFdiffusion image ships 2.0.2, and the only symptom is `Numpy is not available` raised
from inside the diffusion loop.

## Prior work

The 152 privileged scaffolds come from Liu B et al., *Science*,
[doi:10.1126/science.adv0185](https://doi.org/10.1126/science.adv0185) — 10 of them could not be
transplanted onto the NY-ESO-1 groove without backbone clashes, leaving 142 (147 for KF9).

Their stated rationale for partial diffusion is reduced computational cost; no per-design
comparison against free diffusion is reported. Their readout is pooled yeast display over
200–12,000 designs per target, which has no per-design denominator and so cannot measure hit-rate
efficiency in either direction. That gap is what this benchmark begins to explore.

Their Data S1 is consistent with the donor-heterogeneity result here: 8 of their 18 validated
binders come from partial diffusion of a single scaffold (`mage-513`), which is the sole source of
binders for 3 of their 10 non-MAGE targets.

## License

Code under MIT. The scaffold set and target structures carry their original terms — see
`docs/internal_notes.md` for provenance.
