# De novo pMHC binder design pilot

Does starting RFdiffusion from a *privileged scaffold* (partial diffusion) beat
designing from scratch (free diffusion)? Two targets, five arms, ~3,400 designs.

Consolidated 2026-11-02 from the `g5.2xlarge` pilot instance and the local
analysis session.

## Headline results

| target | arm | per design | per scaffold |
|---|---|---|---|
| NY-ESO-1 `SLLMWITQV` / HLA-A\*02:01 | partial, 142 privileged | **7.9%** 67/852 | 18.3% 26/142 |
| | free diffusion | 1.9% 15/800 | 5.5% 11/200 |
| | partial from KH46 (positive control) | 34.2% 203/594 | 67.7% 67/99 |
| SARS-CoV-2 KF9 `KAYNVTQAF` / HLA-C\*12:02 | partial, same 142 | **3.6%** 21/588 | 8.2% 12/147 |
| | free diffusion | 4.2% 25/600 | 11.3% 17/150 |

Borrowing wins 3.8× per design on NY-ESO-1 and **does not replicate** on KF9
(0.86×, p = 0.65). Per-scaffold rates are on a common 4-sequence basis; the
NY-ESO-1 arms sampled 6 sequences per scaffold, so raw per-scaffold numbers
there read higher (23.2%, 67/142) and are not comparable across targets.

Four gates: scRMSD < 2.0 Å, binder pLDDT > 92, interface pAE < 7.0,
peptide_focus ≥ 0.202.

## Layout

```
reproduce_env.sh      rebuild the whole environment from pinned commits
                      `bash reproduce_env.sh verify` checks an install, non-destructively
env/                  resolved micromamba lock files + the raw capture from the box
targets/              the three target structures (NY-ESO-1, KF9, KH46 control)
data/                 every metrics TSV, both scaffold manifests, specificity tables
pipeline/             code that ran on the GPU box and produced data/
analysis/             code that ran locally and turned data/ into the findings
deck/                 source HTML for the slide deck
```

## data/

One row per design. Columns: `backbone_id, seq_id, description, scrmsd,
plddt_binder, pae_interaction, peptide_focus, buried_peptide, buried_mhc,
plddt_target, target_aligned_rmsd`.

| file | arm |
|---|---|
| `liu_metrics.tsv` | NY-ESO-1, partial diffusion from 142 privileged scaffolds |
| `free_metrics.tsv` | NY-ESO-1, free diffusion |
| `kh46self_metrics.tsv` | NY-ESO-1, partial from KH46 (positive control) |
| `kf9_metrics.tsv` | KF9, partial diffusion from the same 142 |
| `kf9free_metrics.tsv` | KF9, free diffusion |
| `ny1b04_metrics.tsv` | dropped arm — the source turned out to be an idealised scaffold, not a second NY-ESO-1 binder |
| `kh46ctl_metrics.tsv` | fold-conditioning arm (shape description, not docked pose) — 3.3%, indistinguishable from free |
| `transfer_metrics.tsv`, `test_metrics.tsv`, `kf9smoke_metrics.tsv` | smoke tests and intermediates |
| `manifest_liu.tsv`, `manifest_kf9.tsv` | per-scaffold transplant record: `source_target`, `source_allele`, `source_peptide`, clash counts, transplant-stage contacts. **This is the join key for donor provenance.** |
| `liu_specificity.tsv`, `free_specificity.tsv`, `kh46self_specificity.tsv` | ProteinMPNN holo−apo log-prob margins, on-target vs 4 decoy peptides |

152 → 142 scaffolds: 10 could not be transplanted onto the NY-ESO-1 groove
without backbone clashes (5 more fail on KF9, giving 147).

## pipeline/

Ran on the box. `scripts/` holds the logic, `runs/` the per-arm drivers.

- `scripts/score_specificity.py` — the specificity filter. Relabels peptide
  residue *names* on a fixed backbone and runs one batched ProteinMPNN pass per
  sequence condition (`--conditional_probs_only`). No structural remodelling.
- `scripts/dock_scaffolds.py` — superposes Baker's scaffolds onto our groove
- `scripts/prep_transfer.py`, `prep_self_control*.py` — arm setup
- `scripts/to_two_chain.py` — collapses output to the binder/target layout the
  metrics collector expects (`PMHC_PEPTIDE` env var selects the peptide)
- `scripts/compare_arms.py`, `collect_metrics.py` — gate evaluation, TSV output
- `runs/run_rfdiffusion*.sh`, `run_arm.sh`, `run_scaffold_arm*.sh` — the arms
- `setup/` — original provisioning, superseded by `reproduce_env.sh`

## analysis/

Ran locally against `data/`. These produced the numbers in the deck, and are
**not** reproducible from `pipeline/` alone.

Deck figures and tables:

| script | produces |
|---|---|
| `src.py` | hit rate by scaffold **source locus**, per target |
| `conf.py` | donor breakdown by source target; transplant-stage features of winners vs losers |
| `donor.py` | the donor-bias table (PHOX2B 52.9% … PAP 0/22) |
| `mw.py` | fold vs interface metrics, Mann-Whitney, all designs and survivors |
| `mw2.py` | the same conditioned on the design contacting the peptide at all |
| `gen_donors.py` | the screen-shallow-then-deepen SVG |
| `gen_metrics.py` | inline-styled metric tables for the deck |
| `analyse_design.py` | contact and interaction analysis for one complex (the scaffold-128 figure) |
| `clash_detail.py` | which atoms clash — the 152 → 142 trim |
| `score_specificity.py` | local copy of the box scorer |

Figure pipeline (PyMOL → PowerPoint):
`turntable.py` (360° ray-traced frames) → `export_parts.py` (per-component OBJ)
→ `build_glb.py` (single .glb for PowerPoint's 3D model insert).

Provenance and composition work:
`provenance.py`, `prov_full.py`, `final_prov.py`, `rank_hlac.py`,
`cognate_check.py`, `ctarget.py`, `peps.py`, `inspect6ulr.py`, `pdftxt.py`.
This is the line of work that established the HLA-C / PHOX2B source identity —
later **confirmed outright** by the paper's own Data S1, which lists
`phox2b-5` and `phox2b-11` as `C*07:02` / `QYNPIRTTF`.

Interface composition: `bbsc.py`, `bbsc2.py` (backbone vs side-chain contact
split, and whether PHOX2B's higher side-chain fraction is a Trp5 packing
artefact).

Remaining scripts — `arm_stats.py`, `depth.py`, `diversity.py`, `exptest.py`,
`final.py`, `plddt.py`, `survey.py`, `tails.py`, `tailpos.py` — are earlier
exploratory analyses with no docstring. Read the source before trusting any
output; some were one-off checks superseded by the scripts above.

## structures/

The 337 structures that survive all four gates, plus the scaffold-128 model
behind the interface figure. **These are not regenerable.** Diffusion is
stochastic, so re-running the pipeline produces different designs, not these.

```
structures/survivors/liu/        67   NY-ESO-1, partial from privileged scaffolds
structures/survivors/free/       15   NY-ESO-1, free diffusion
structures/survivors/kh46self/  203   NY-ESO-1, partial from KH46
structures/survivors/kf9/        21   KF9, partial
structures/survivors/kf9free/    25   KF9, free diffusion
structures/survivors/ny1b04/      5   dropped arm
structures/scaf128.pdb                the design in the structural figure
```

Filenames match the `description` column in `data/*_metrics.tsv`, so the
metrics join straight onto the coordinates.

## Not included

Left on the pilot instance, which was **stopped, not terminated** (2026-11-02):
the ~4,100 AF2 predictions for designs that failed a gate, all RFdiffusion
backbones, the ProteinMPNN sequence pools, and the specificity caches —
roughly 16 GB, in `~/pmhc/` on the 194 GB root volume.

They are still recoverable: restart the instance and copy them off. They are
not re-creatable, since diffusion is stochastic. If the instance is ever
terminated with `DeleteOnTermination` set, they go with the volume.

## Environment

`g5.2xlarge` (NVIDIA A10G), Ubuntu 22.04.5, driver 580.126.09.
RFdiffusion `86507b6` in Docker; dl_binder_design `cafa385` and ProteinMPNN
`8907e66` under micromamba; AF2 params `2022-12-06`.

The one non-obvious thing: the working image is `rfdiffusion:fixed`, not the
stock build. torch 1.12.1 is compiled against the numpy 1.x C ABI while the
official image ships 2.0.2, and the only symptom is `Numpy is not available`
raised from inside the diffusion loop. `reproduce_env.sh` rebuilds that layer
and its verify mode executes python inside the container to confirm it.

## Reference

Liu B et al., *Science* — doi:10.1126/science.adv0185. Data S1 (18 validated
binders, 11 targets) flags 8 of 18 as partial diffusion from one scaffold,
`mage-513`, which is the sole source of binders for 3 of their 10 non-MAGE
targets. Their pooled yeast-display readout has no per-design denominator, so
it cannot measure hit-rate efficiency in either direction.
