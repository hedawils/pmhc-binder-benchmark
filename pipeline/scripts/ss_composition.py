#!/usr/bin/env python
"""Secondary-structure composition of a set of RFdiffusion backbones.

    usage:  ss_composition.py <label>=<dir_of_ss_pt> [<label>=<dir> ...]

Answers a question the pass rates cannot: what topologies is RFdiffusion
actually producing? Free binder design is unconstrained in principle, but
RFdiffusion's sampling prior is not uniform over fold space, and if the
free arm turns out to be as helical as a fold-conditioned one, then
"conditioning" is a much narrower manipulation than its name suggests.

Reads the _ss.pt tensors emitted by RFdiffusion's make_secstruc_adj.py, so
whatever is compared here was labelled by the same assignment -- P-SEA,
which uses CA positions only and is therefore valid on side-chain-free
diffusion backbones.

Class indices come from ss_to_tensor in that helper: 0=Helix, 1=Strand,
2=Loop, 3=Mask.
"""

from __future__ import annotations

import glob
import statistics as st
import sys
from pathlib import Path

import torch

LABEL = {0: "H", 1: "E", 2: "L", 3: "-"}


def runs(seq: str, want: str) -> int:
    """Number of maximal runs of `want` -- i.e. how many helices or strands."""
    count, prev = 0, None
    for char in seq:
        if char == want and prev != want:
            count += 1
        prev = char
    return count


def load(directory: str) -> list[dict]:
    rows = []
    for path in sorted(glob.glob(f"{directory}/*_ss.pt")):
        ss = torch.load(path, map_location="cpu")
        seq = "".join(LABEL.get(int(v), "?") for v in ss.tolist())
        length = len(seq)
        if not length:
            continue
        rows.append({
            "name": Path(path).name,
            "len": length,
            "helix": seq.count("H") / length,
            "strand": seq.count("E") / length,
            "loop": seq.count("L") / length,
            "n_helix": runs(seq, "H"),
            "n_strand": runs(seq, "E"),
        })
    return rows


def summarise(label: str, rows: list[dict]) -> None:
    if not rows:
        print(f"\n=== {label} ===\n  no _ss.pt files found")
        return

    n = len(rows)
    print(f"\n=== {label}  (n={n}) ===")

    lens = [r["len"] for r in rows]
    print(f"  length           min {min(lens)}  median {int(st.median(lens))}"
          f"  max {max(lens)}")

    for key in ("helix", "strand", "loop"):
        vals = [r[key] for r in rows]
        print(f"  {key:<8} frac    mean {st.mean(vals):.3f}"
              f"  median {st.median(vals):.3f}  max {max(vals):.3f}")

    # Strand content is the discriminating question: an all-alpha prior means
    # free design is not exploring beta topologies at all, whatever its
    # nominal freedom.
    for cutoff in (0.0, 0.05, 0.10, 0.20):
        hit = sum(1 for r in rows if r["strand"] > cutoff)
        tag = "any strand" if cutoff == 0 else f"strand > {cutoff:.0%}"
        print(f"  {tag:<16} {hit:>5} / {n}  ({hit / n:>6.1%})")

    hist: dict[int, int] = {}
    for r in rows:
        hist[r["n_helix"]] = hist.get(r["n_helix"], 0) + 1
    counts = sorted(r["n_helix"] for r in rows)
    print(f"  helices/design   median {int(st.median(counts))}"
          f"  range {min(counts)}-{max(counts)}")
    print("  histogram        " +
          "  ".join(f"{k}H:{v}" for k, v in sorted(hist.items())))


def main() -> None:
    args = sys.argv[1:]
    if not args or any("=" not in a for a in args):
        sys.exit(__doc__)
    for arg in args:
        label, _, directory = arg.partition("=")
        summarise(label, load(directory))


if __name__ == "__main__":
    main()
