#!/usr/bin/env python
"""Pick the scaffold fold(s) to condition on, standing in for an expert eye.

    usage:  select_scaffold.py <pdb_dir> <ss_dir> <out_dir> [how_many]

WHY THIS EXISTS

The published procedure on this target generated 50 unconditional monomers,
SELECTED ONE BY EYE, and conditioned 100 backbones on that single fold. Our
fold arm generated 50 and conditioned across all 50, which is not the same
experiment: with random scaffold sampling, 100 backbones over 50 folds gives
any individual fold about 2 backbones instead of 100 -- a fiftieth of the
exposure. A null result from that would mostly measure "is a RANDOM fold any
good for this target", a much weaker question whose answer is probably no.

The selection step is the active ingredient and it is unreproducible as
published, because "it looked right" is not a procedure. This script
substitutes an EXPLICIT criterion so the step can at least be run and
reported. It is a stand-in, not a reconstruction of anyone's judgement.

THE CRITERION, AND WHY IT IS NOT CIRCULAR

Two rules, both derived from the TARGET's geometry and the published class,
never from KH46's solution -- conditioning on KH46's own blueprint is a
separate arm, and reusing it here would collapse the two:

  1. helix count in [3, 5].  Every binder validated on a pMHC groove is a
     small all-alpha bundle: the published TCR mimic has 4 helices, KH46 4,
     NY1-B04 5. This is a class constraint, not a shape match.

  2. maximise span along the longest principal axis, subject to a floor.
     The groove is long and narrow -- a 9mer peptide is ~25 A end to end and
     the alpha1/alpha2 helices flank it -- so a binder has to reach across
     it. A compact globule cannot, however good its topology looks.

What this cannot claim: that this criterion is what the published authors
were doing, or that it beats random. The arm is split so that is MEASURED
rather than assumed -- half the backbones on the selected fold, half across
the whole library, same budget.
"""

import shutil
import sys
from pathlib import Path

import numpy as np
import torch

# ss_to_tensor in RFdiffusion's make_secstruc_adj.py: 0=Helix 1=Strand 2=Loop
LABEL = {0: "H", 1: "E", 2: "L", 3: "-"}

# Helix-count window: the all-alpha bundle class every validated pMHC binder
# belongs to (published TCR mimic 4, KH46 4, NY1-B04 5).
HELIX_MIN, HELIX_MAX = 3, 5

# A 9mer in the groove spans ~25 A; a binder must reach across it and engage
# the flanking helices. Angstroms, along the longest principal axis.
MIN_SPAN = 35.0


def runs(seq: str, want: str) -> int:
    count, prev = 0, None
    for ch in seq:
        if ch == want and prev != want:
            count += 1
        prev = ch
    return count


def ca_coords(pdb: Path) -> np.ndarray:
    out = []
    for ln in pdb.read_text().splitlines():
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
            out.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return np.array(out)


def span(xyz: np.ndarray) -> float:
    """Extent along the longest principal axis, in Angstroms."""
    if len(xyz) < 3:
        return 0.0
    centred = xyz - xyz.mean(axis=0)
    # Right singular vectors are the principal axes, largest first.
    axis = np.linalg.svd(centred, full_matrices=False)[2][0]
    proj = centred @ axis
    return float(proj.max() - proj.min())


def survey(pdb_dir: Path, ss_dir: Path) -> list[dict]:
    rows = []
    for ss_path in sorted(ss_dir.glob("*_ss.pt")):
        stem = ss_path.name[:-6]              # drop "_ss.pt"
        pdb = pdb_dir / f"{stem}.pdb"
        adj = ss_dir / f"{stem}_adj.pt"
        if not pdb.exists() or not adj.exists():
            continue
        try:
            ss = torch.load(ss_path, map_location="cpu")
            seq = "".join(LABEL.get(int(v), "?") for v in ss.tolist())
            xyz = ca_coords(pdb)
        except Exception as exc:                # a corrupt file is not fatal
            print(f"  skipping {stem}: {exc}")
            continue
        if not seq or len(xyz) < 3:
            continue
        rows.append({
            "stem": stem,
            "len": len(seq),
            "n_helix": runs(seq, "H"),
            "n_strand": runs(seq, "E"),
            "helix_frac": seq.count("H") / len(seq),
            "span": span(xyz),
            "ss": ss_path,
            "adj": adj,
        })
    return rows


def main() -> None:
    if len(sys.argv) < 4:
        sys.exit("usage: select_scaffold.py <pdb_dir> <ss_dir> <out_dir> [n]")
    pdb_dir, ss_dir, out_dir = (Path(sys.argv[1]), Path(sys.argv[2]),
                                Path(sys.argv[3]))
    want = int(sys.argv[4]) if len(sys.argv) > 4 else 1

    rows = survey(pdb_dir, ss_dir)
    if not rows:
        sys.exit(f"no usable scaffolds found in {ss_dir} / {pdb_dir}")

    print(f"surveyed {len(rows)} scaffolds")
    print(f"criterion: {HELIX_MIN}-{HELIX_MAX} helices, span >= "
          f"{MIN_SPAN:.0f} A, then longest span wins")

    ok = [r for r in rows
          if HELIX_MIN <= r["n_helix"] <= HELIX_MAX and r["span"] >= MIN_SPAN]
    print(f"  {len(ok)} of {len(rows)} pass the class + span filter")
    if not ok:
        # Report WHY before failing, so the morning read explains the fallback.
        helices = sorted(r["n_helix"] for r in rows)
        spans = sorted(r["span"] for r in rows)
        print(f"  helix counts seen: min {helices[0]} median "
              f"{helices[len(helices) // 2]} max {helices[-1]}")
        print(f"  spans seen: min {spans[0]:.1f} median "
              f"{spans[len(spans) // 2]:.1f} max {spans[-1]:.1f} A")
        sys.exit("no scaffold satisfies the criterion")

    ok.sort(key=lambda r: -r["span"])

    print(f"\n{'rank':<6}{'scaffold':<16}{'len':>5}{'helices':>9}"
          f"{'strands':>9}{'helix%':>8}{'span A':>9}")
    for i, r in enumerate(ok[:12], 1):
        mark = "  <-- selected" if i <= want else ""
        print(f"{i:<6}{r['stem']:<16}{r['len']:>5}{r['n_helix']:>9}"
              f"{r['n_strand']:>9}{100 * r['helix_frac']:>7.1f}%"
              f"{r['span']:>9.1f}{mark}")

    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob("*.pt"):
        stale.unlink()
    for r in ok[:want]:
        shutil.copy2(r["ss"], out_dir / r["ss"].name)
        shutil.copy2(r["adj"], out_dir / r["adj"].name)

    names = ", ".join(r["stem"] for r in ok[:want])
    print(f"\nselected {want}: {names}")
    print(f"wrote {len(list(out_dir.glob('*.pt')))} tensors -> {out_dir}")


if __name__ == "__main__":
    main()
