#!/usr/bin/env python
"""Build the ny1b04-onto-our-target input: the SECOND self-transfer source.

WHY A SECOND ONE

Arm 1b used a single structure, KH46. One source cannot tell you whether
34.2% is "what this target allows" or "what KH46 allows", and the deck was
leaning on the first reading. This arm supplies n=2 so the ceiling claim is
either corroborated or killed.

THE PEPTIDE WRINKLE, AND WHY SUPERPOSITION RATHER THAN RELABELLING

ny1b04 was solved against native NY-ESO-1 157-165, SLLMWITQ*C*. Every other
arm targets the 9V anchor analog, SLLMWITQ*V*. Rather than carry a one-residue
target difference through the whole arm, this script superposes ny1b04's
platform onto OUR standard target and keeps our peptide, so the target is
atom-identical to arms 1b/3/5 and the comparison is exact.

That is defensible on our own contact data: P9 is the C-terminal anchor,
buried in the F pocket pointing down into the groove, and it makes ZERO
side-chain contacts with the binder in any design we measured. C9 vs V9 is
invisible from the binder's side. The cost is that the binder is now sitting
on a groove 0.x A from the one it was solved on -- reported below as the
platform RMSD, which is the number to check before trusting this arm.

Unlike prep_self_control.py, this script DOES move atoms: the binder is
rigid-body transformed. The platform and peptide are copied from our target
untouched, so only the binder's frame changes, and the internal geometry of
the binder is preserved exactly (a rotation+translation cannot alter it).

    usage:  prep_self_control2.py <ny1b04_control.pdb> <target.pdb> <out.pdb>
"""

import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

PLATFORM_LAST = 180          # chain B residues <= this are the alpha-1/alpha-2 platform
BACKBONE = (" N  ", " CA ", " C  ", " O  ")
CLASH = 2.5                  # A, backbone-backbone


def read_chains(path):
    chains = OrderedDict()
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[16] not in (" ", "A"):
                continue
            chains.setdefault(line[21], OrderedDict()) \
                  .setdefault(int(line[22:26]), []).append(line)
    return chains


def seq(residues):
    return "".join(THREE_TO_ONE.get(v[0][17:20].strip(), "X")
                   for v in residues.values())


def ca(residues):
    """resnum -> CA coordinate."""
    out = {}
    for num, lines in residues.items():
        for ln in lines:
            if ln[12:16] == " CA ":
                out[num] = np.array([float(ln[30:38]), float(ln[38:46]),
                                     float(ln[46:54])])
    return out


def kabsch(P, Q):
    """Rotation+translation taking P onto Q. Returns R, t, rmsd."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = qc - R @ pc
    rmsd = float(np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean()))
    return R, t, rmsd


def main():
    if len(sys.argv) != 4:
        sys.exit("usage: prep_self_control2.py <ny1b04.pdb> <target.pdb> <out.pdb>")
    src, tgt, dst = map(Path, sys.argv[1:4])
    for p in (src, tgt):
        if not p.is_file():
            sys.exit(f"no such input: {p}")

    s, t = read_chains(src), read_chains(tgt)
    if "A" not in s or "B" not in s:
        sys.exit(f"{src.name}: expected chains A (binder) and B (MHC+peptide)")
    if "A" not in t or "C" not in t:
        sys.exit(f"{tgt.name}: expected chains A (platform) and C (peptide)")

    binder = s["A"]
    s_plat = OrderedDict((n, v) for n, v in sorted(s["B"].items())
                         if n <= PLATFORM_LAST)
    s_pep = OrderedDict((n, v) for n, v in sorted(s["B"].items())
                        if n > PLATFORM_LAST)
    t_plat, t_pep = t["A"], t["C"]

    print(f"{src.name}:  binder {len(binder)}  platform {len(s_plat)}  "
          f"peptide {len(s_pep)} = {seq(s_pep)}")
    print(f"{tgt.name}:  platform {len(t_plat)}  peptide {len(t_pep)} = {seq(t_pep)}")

    if seq(s_plat) != seq(t_plat):
        n = sum(a != b for a, b in zip(seq(s_plat), seq(t_plat)))
        sys.exit(f"platform sequences differ at {n} positions; this script "
                 "assumes the same construct and matches by residue number")
    print(f"  platforms: {len(t_plat)} residues, sequence identical")

    s_ca, t_ca = ca(s_plat), ca(t_plat)
    common = sorted(set(s_ca) & set(t_ca))
    if len(common) < 170:
        sys.exit(f"only {len(common)} platform CA pairs; need >= 170")
    P = np.array([s_ca[n] for n in common])
    Q = np.array([t_ca[n] for n in common])
    R, tr, rmsd = kabsch(P, Q)
    print(f"  superposed on {len(common)} platform CA: RMSD {rmsd:.3f} A")
    if rmsd > 1.5:
        sys.exit(f"platform RMSD {rmsd:.3f} A is too high to trust; the two "
                 "structures are not the same groove in the same state")

    # Rigid-body the binder across. Internal geometry is preserved by
    # construction; only its frame changes.
    moved = OrderedDict()
    for num, lines in binder.items():
        new = []
        for ln in lines:
            xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            x, y, z = R @ xyz + tr
            new.append(ln[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + ln[54:])
        moved[num] = new

    # Did the transfer land the binder on top of the groove?
    tgt_bb = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                       for res in (t_plat, t_pep) for l in
                       [x for v in res.values() for x in v]
                       if l[12:16] in BACKBONE])
    bnd_bb = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                       for v in moved.values() for l in v
                       if l[12:16] in BACKBONE])
    d = np.linalg.norm(bnd_bb[:, None, :] - tgt_bb[None, :, :], axis=2)
    clashes = int((d < CLASH).sum())
    print(f"  backbone-backbone clashes < {CLASH} A: {clashes}")
    print(f"  closest binder-target backbone approach: {d.min():.2f} A")
    if clashes:
        print("  NOTE: partial diffusion repairs backbone clashes; "
              "side-chain overlap is repacked by MPNN and is not counted here.")

    # Emit A=binder B=platform C=peptide, each renumbered from 1.
    out, serial = [], 1
    groups = [("A", list(moved.values())),
              ("B", list(t_plat.values())),
              ("C", list(t_pep.values()))]
    for dest, residues in groups:
        for index, lines in enumerate(residues, start=1):
            for ln in lines:
                out.append(f"{ln[:6]}{serial:5d}{ln[11:21]}{dest}{index:4d}"
                           f"{ln[26]}{ln[27:].rstrip()}\n")
                serial += 1
        out.append(f"TER   {serial:5d}\n")
        serial += 1
        print(f"  chain {dest}: {len(residues)} residues")
    out.append("END\n")

    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("".join(out))
    print(f"wrote {dst}  (A=binder B=MHC C=peptide {seq(t_pep)})")


if __name__ == "__main__":
    main()
