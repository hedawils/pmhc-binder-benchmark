#!/usr/bin/env python
"""Side-by-side funnel for the free-design and fold-conditioned arms.

Both arms share target, hotspots, MPNN settings and AF2 settings, so
conditioning is the single moving variable and the pass rates are directly
comparable. Scaffold length is matched by construction (see
env/make_scaffolds.sh), which is what makes that claim defensible -- the
shipped PPI library would have moved length at the same time.

Also runs the pose-limited vs sequence-limited diagnostic, which asks a
question the pass rates alone cannot: when designs fail, is the backbone
wrong or is the sequence wrong?

    usage:  compare_arms.py <free_metrics.tsv> <fold_metrics.tsv>
            compare_arms.py <metrics.tsv>                 # one arm only
"""

from __future__ import annotations

import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pmhc_pilot.funnel import (  # noqa: E402
    AF2_THRESHOLDS,
    STAGE_ORDER,
    format_funnel,
    funnel,
    passes,
    unique_rejections,
)

NUMERIC = ("scrmsd", "plddt_binder", "pae_interaction", "peptide_focus")


def load(path: Path) -> list[dict]:
    """Read a metrics table, coercing the numeric columns to float."""
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        rec = dict(zip(header, line.split("\t")))
        for key in NUMERIC:
            try:
                rec[key] = float(rec[key])
            except (ValueError, KeyError):
                rec[key] = math.nan
        rows.append(rec)
    return rows


def pose_vs_sequence(rows: list[dict]) -> None:
    """Within-backbone spread of pae_interaction.

    ProteinMPNN cannot move the binder: the backbone fixes the docking pose,
    and every sequence from one backbone inherits it. So the spread of
    pae_interaction *within* a backbone separates two failure modes that an
    aggregate pass rate cannot.

      tight, at the ceiling   -> pose-limited. The sequences are irrelevant;
                                 the pose is wrong. Deeper MPNN sampling
                                 cannot help. Needs different backbones.
      wide, with near-misses  -> sequence-limited. The pose can work and
                                 MPNN has not found the sequence for it yet.
                                 Deeper sampling is the right spend.

    This is the measurement that decides where the next GPU-hours go, so it
    is reported whether or not anything passed.
    """
    by_bb: dict[str, list[float]] = {}
    for r in rows:
        v = r.get("pae_interaction", math.nan)
        if not math.isnan(v):
            by_bb.setdefault(r["backbone_id"], []).append(v)

    multi = {k: v for k, v in by_bb.items() if len(v) > 1}
    if not multi:
        print("  no backbone has >1 scored sequence; diagnostic unavailable")
        return

    spreads = [max(v) - min(v) for v in multi.values()]
    bests = [min(v) for v in multi.values()]
    cutoff = AF2_THRESHOLDS["pae_interaction"][1]

    print(f"  backbones with >1 sequence: {len(multi)}")
    print(f"  within-backbone pae spread: median {statistics.median(spreads):.2f}, "
          f"max {max(spreads):.2f}")
    print(f"  best-per-backbone pae:      median {statistics.median(bests):.2f}, "
          f"min {min(bests):.2f}")

    spread = statistics.median(spreads)
    passing = sum(1 for b in bests if b < cutoff)

    # How far a backbone could plausibly travel on sequence choice alone.
    # Scaled by the OBSERVED within-backbone spread rather than a fixed
    # window: where sequence barely moves pae, "close to the cutoff" has to
    # mean close on that scale. Two typical spreads of headroom.
    reach = max(spread, 0.5) * 2
    rescuable = sum(1 for b in bests if cutoff <= b < cutoff + reach)
    hopeless = len(bests) - passing - rescuable

    print(f"  backbones with a passing sequence (pae < {cutoff}): {passing}")
    print(f"  within reach of the cutoff (+{reach:.1f}):"
          f"{'':<16}{rescuable}")
    print(f"  out of reach:{'':<38}{hopeless}")

    # Proportions, not a binary label. The earlier rule was
    #   POSE-LIMITED iff median spread < 1.0 AND nothing passed
    # which forced any arm with a single hit into "sequence-limited" and
    # recommended deeper MPNN sampling -- the wrong spend when most
    # backbones sit far above the cutoff and sequence moves pae by ~1. The
    # realistic case is bimodal: a few backbones work, the rest cannot be
    # rescued at any sampling depth, and both facts matter.
    frac_hopeless = hopeless / len(bests)
    if frac_hopeless >= 0.6:
        verdict = (f"MOSTLY POSE-LIMITED -- {hopeless}/{len(bests)} backbones "
                   f"out of reach (sequence moves pae by ~{spread:.1f}; their "
                   f"best sits at {statistics.median(bests):.1f} vs cutoff "
                   f"{cutoff}). MORE BACKBONES beat deeper MPNN sampling.")
    elif rescuable >= passing and rescuable > 0:
        verdict = (f"SEQUENCE-LIMITED -- {rescuable} backbones sit within "
                   f"{reach:.1f} of the cutoff; deeper MPNN sampling is the "
                   f"right spend.")
    else:
        verdict = (f"MIXED -- {passing} passing, {rescuable} in reach, "
                   f"{hopeless} out of reach; no single bottleneck.")
    print(f"  verdict: {verdict}")
    if passing and frac_hopeless >= 0.6:
        print(f"  note: the {passing} passing backbone(s) are a separate "
              f"population, not the tail of a continuum.")


def report(name: str, rows: list[dict]) -> None:
    print(f"\n{'=' * 72}\n{name}  (n={len(rows)})\n{'=' * 72}")
    if not rows:
        print("  no rows")
        return

    print(format_funnel(funnel(rows, AF2_THRESHOLDS)))

    print("\n-- per-filter marginal pass rates (each filter alone) --")
    for metric in STAGE_ORDER:
        ok = sum(1 for r in rows if passes(r, metric, AF2_THRESHOLDS))
        print(f"  {metric:<18} {ok:>5} / {len(rows):<5} {ok / len(rows):>7.1%}")

    print("\n-- unique rejections (designs each filter rejects alone) --")
    for metric, n in unique_rejections(rows, AF2_THRESHOLDS).items():
        print(f"  {metric:<18} {n:>5}")

    print("\n-- pose-limited vs sequence-limited --")
    pose_vs_sequence(rows)


def main() -> None:
    paths = [Path(a) for a in sys.argv[1:]]
    if not paths:
        sys.exit(__doc__)

    loaded = []
    for p in paths:
        if not p.exists():
            print(f"skipping missing {p}")
            continue
        loaded.append((p.stem, load(p)))

    if not loaded:
        sys.exit("no metrics tables found")

    for name, rows in loaded:
        report(name, rows)

    # Any number of arms, not just two -- the run has four, and four
    # individual funnels with no side-by-side is the one view that cannot
    # answer "which route is more efficient".
    if len(loaded) >= 2:
        print(f"\n{'=' * 72}\nARM COMPARISON\n{'=' * 72}")
        short = [n.replace("_metrics", "")[:17] for n, _ in loaded]
        header = f"{'stage':<16}" + "".join(f"{n:>18}" for n in short)
        print(header)
        print("-" * len(header))
        funnels = [dict((s.metric, s) for s in funnel(r, AF2_THRESHOLDS))
                   for _, r in loaded]
        for metric in STAGE_ORDER:
            cells = ""
            for f in funnels:
                s = f[metric]
                cells += (f"{s.seqs_out:>5}/{s.seqs_in:<5}"
                          f"{s.seq_pass_rate:>7.1%}")
            print(f"{metric:<16}{cells}")
        print("\nAll-filter survivors:")
        for name, rows in loaded:
            n = sum(1 for r in rows
                    if all(passes(r, m, AF2_THRESHOLDS) for m in STAGE_ORDER))
            print(f"  {name:<22} {n:>5} / {len(rows)}")


if __name__ == "__main__":
    main()
