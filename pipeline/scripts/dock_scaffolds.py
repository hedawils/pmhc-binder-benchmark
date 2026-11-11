#!/usr/bin/env python
"""Dock Liu et al.'s privileged scaffolds onto our target, ready for partial diffusion.

    usage:  dock_scaffolds.py <scaffold_dir> <target_pdb> <out_dir> [--no-sasa]

WHAT THIS IS

Liu et al. (Science 2025, doi 10.1126/science.adv0185) released 152 designed
pMHC-I binders as docked complexes: chain A is the binder, chain B is the
pMHC it was designed against, with the MHC and peptide FUSED under one
continuous numbering. Their own pipeline repurposes these by superimposing
chain B onto a new target, keeping the binder, then partial-diffusing
(their `align_chainB.py`, which does this in PyMOL).

This is the same operation in Biopython, with the checks their script does
not do: the superposition RMSD is recorded, backbone clashes are measured
separately from side-chain clashes, and every scaffold's source allele is
carried through to a manifest so the arm can be split by it afterwards.

WHY THE BINDER MOVES, NOT THE TARGET

prep_transfer.py moves the target under a fixed binder, because there it was
one binder and one new target. Here one target is shared by 152 scaffolds, so
the target is held fixed and each binder is brought into its frame. Same
rigid-body operation, fewer copies of the target.

THE SUPERPOSITION IS BY RELATIVE INDEX, AND THAT IS CHECKED

Chain B's numbering starts after the binder (scaffold-dependent), and its
peptide is fused to its MHC with no chain break, so there is no reliable way
to split them by numbering alone. Instead the first SUPERPOSE_N residues of
chain B are used -- far from the C-terminal peptide, so the peptide is never
included -- matched by relative index against the target's MHC.

Relative index assumes the two MHCs align residue-for-residue without
insertions. That holds for the same allele and approximately across HLA-A/B/C,
but it is an assumption, so the RMSD is reported and anything above
MAX_RMSD is rejected rather than silently docked out of register. 42% of the
library carries an MHC that is none of the paper's four alleles, so this
guard is doing real work.
"""

import sys
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, Superimposer

# Residues of chain B used for superposition. Chain B is 182-191 residues of
# MHC+peptide; 170 stays clear of the peptide for every scaffold in the set.
SUPERPOSE_N = 170

# Two copies of a class I platform should superimpose well below this. Above
# it, the residue-for-residue assumption has broken down.
MAX_RMSD = 2.5

CLASH = 2.5      # heavy-atom steric clash, Angstroms
CONTACT = 4.5
BACKBONE = {"N", "CA", "C", "O"}

# The paper's targets, keyed by the peptide fused to the end of chain B.
# Lets the arm be split by the allele each scaffold was DESIGNED against,
# which is the whole question of whether a privileged scaffold is portable.
KNOWN = {
    "RMFPNAPYL": ("WT1", "A*02:01"),
    "KLTPLCVTL": ("HIV", "A*02:01"),
    "GLMWLSYFV": ("SARS-CoV1", "A*02:01"),
    "LLWNGPIAV": ("YFV", "A*02:01"),
    "TLMSAMTNL": ("PAP", "A*02:01"),
    "YLEPGPVTA": ("gp100", "A*02:01"),
    "ELAGIGILTV": ("MART-1", "A*02:01"),
    "ALYVDSLFFL": ("PRAME", "A*02:01"),
    "EVDPIGHLY": ("MAGE-A3", "A*01:01"),
    "ESDPIVAQY": ("Titin", "A*01:01"),
    "QYNPIRTTF": ("PHOX2B", "C*07:02"),
    "TTAPFLSGK": ("CTNNb1-S5F", "A*03:01"),
    "ALHGGWTTK": ("seed scaffold", "A*03:01"),
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def residues(chain):
    return [r for r in chain if r.id[0] == " "]


def seq(res):
    return "".join(THREE_TO_ONE.get(r.get_resname(), "X") for r in res)


def heavy(res, backbone_only=False):
    out = []
    for r in res:
        for a in r:
            if a.element == "H" or a.get_altloc() not in (" ", "A"):
                continue
            if backbone_only and a.get_name() not in BACKBONE:
                continue
            out.append(a)
    return out


def contacts(a_atoms, b_atoms):
    if not a_atoms or not b_atoms:
        return float("nan"), 0, 0
    a = np.array([x.coord for x in a_atoms])
    b = np.array([x.coord for x in b_atoms])
    d = np.linalg.norm(a[:, None] - b[None], axis=-1)
    return float(d.min()), int((d < CLASH).sum()), int((d < CONTACT).sum())


def identify(chain_b_res):
    """Which pMHC was this scaffold designed against?"""
    tail = seq(chain_b_res)
    for n in (10, 9, 8):
        if tail[-n:] in KNOWN:
            name, allele = KNOWN[tail[-n:]]
            return tail[-n:], name, allele
    return "", "unknown", "unknown"


def emit(fh, res, chain_id, start, serial):
    for offset, r in enumerate(res):
        for a in r:
            if a.element == "H" or a.get_altloc() not in (" ", "A"):
                continue
            x, y, z = a.coord
            nm = a.get_name()
            nm = f" {nm:<3}" if len(nm) < 4 else nm
            fh.write(f"ATOM  {serial:>5} {nm} {r.get_resname():>3} "
                     f"{chain_id}{start + offset:>4}    "
                     f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00\n")
            serial += 1
    fh.write("TER\n")
    return serial


def main() -> None:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    do_sasa = "--no-sasa" not in sys.argv
    if len(args) < 3:
        sys.exit("usage: dock_scaffolds.py <scaffold_dir> <target_pdb> "
                 "<out_dir> [--no-sasa]")
    scaf_dir, target_pdb, out_dir = Path(args[0]), Path(args[1]), Path(args[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    focus_fn = None
    if do_sasa:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
        try:
            from pmhc_pilot.interface import analyse_interface
            focus_fn = analyse_interface
        except Exception as exc:
            print(f"peptide_focus unavailable ({exc}); continuing without it")

    parser = PDBParser(QUIET=True)
    tgt = parser.get_structure("t", target_pdb)[0]
    # prep_target.py writes the trimmed target as A = MHC(180), C = peptide(9).
    t_mhc = residues(tgt["A"])
    t_pep = residues(tgt["C"])
    t_all = heavy(t_mhc + t_pep)
    print(f"target {target_pdb.name}: MHC {len(t_mhc)}  peptide "
          f"{len(t_pep)} ({seq(t_pep)})")

    rows = []
    for path in sorted(scaf_dir.glob("scaf*.pdb"),
                       key=lambda p: int("".join(c for c in p.stem
                                                 if c.isdigit()) or 0)):
        rec = {"scaffold": path.stem, "accepted": 0, "reason": ""}
        try:
            m = parser.get_structure("s", path)[0]
            binder = residues(m["A"])
            chain_b = residues(m["B"])
        except Exception as exc:
            rec["reason"] = f"parse failed: {exc}"
            rows.append(rec)
            continue

        pep, name, allele = identify(chain_b)
        rec.update(source_peptide=pep or "-", source_target=name,
                   source_allele=allele, binder_len=len(binder),
                   chainB_len=len(chain_b))

        # Search a small N-terminal offset rather than assuming chain B
        # starts at MHC residue 1.
        #
        # Without this, 80 of 152 scaffolds were rejected at RMSDs clustered
        # between 3.78 and 4.14 A -- far too tight a spread to be genuine
        # allele divergence, which would scatter. The cause is that the
        # library's MHCs are trimmed differently at the N-terminus: the
        # 189-residue chain B starts GSHSMRY..., the 183-residue one starts
        # PHSMRY..., one residue later. A constant register error of a few
        # residues produces exactly this signature, and it is a plumbing
        # artefact, not a property of the scaffold.
        best = None
        for off in range(-4, 5):
            n = min(SUPERPOSE_N, len(chain_b) - max(off, 0),
                    len(t_mhc) + min(off, 0))
            if n < 100:
                continue
            ti = t_mhc[max(-off, 0):max(-off, 0) + n]
            si = chain_b[max(off, 0):max(off, 0) + n]
            f = [r["CA"] for r in ti if "CA" in r]
            mv = [r["CA"] for r in si if "CA" in r]
            if len(f) != len(mv) or len(f) < 100:
                continue
            s = Superimposer()
            s.set_atoms(f, mv)
            if best is None or s.rms < best[0]:
                best = (float(s.rms), off, s)
        if best is None:
            rec["reason"] = "too few CA pairs at any offset"
            rows.append(rec)
            continue
        rec["rmsd"], rec["offset"] = round(best[0], 3), best[1]
        sup = best[2]
        if sup.rms > MAX_RMSD:
            rec["reason"] = f"platform RMSD {sup.rms:.2f} > {MAX_RMSD}"
            rows.append(rec)
            continue

        # Move the binder into the target's frame; the target never moves.
        sup.apply(heavy(binder))
        b_bb = heavy(binder, backbone_only=True)
        _, bb_clash, _ = contacts(b_bb, heavy(t_mhc + t_pep,
                                              backbone_only=True))
        pmin, pclash, pcontact = contacts(heavy(binder), heavy(t_pep))
        mmin, mclash, mcontact = contacts(heavy(binder), heavy(t_mhc))
        rec.update(bb_clashes=bb_clash, pep_clashes=pclash,
                   mhc_clashes=mclash, pep_contacts=pcontact,
                   mhc_contacts=mcontact, pep_min=round(pmin, 2))

        if bb_clash:
            rec["reason"] = f"{bb_clash} BACKBONE clashes -- diffusion sees these"
            rows.append(rec)
            continue
        if pcontact == 0:
            rec["reason"] = "no contact with the peptide at all"
            rows.append(rec)
            continue

        dst = out_dir / f"{path.stem}.pdb"
        with dst.open("w") as fh:
            serial = emit(fh, binder, "A", 1, 1)
            serial = emit(fh, t_mhc, "B", 1, serial)
            emit(fh, t_pep, "C", 1, serial)
            fh.write("END\n")

        if focus_fn is not None:
            try:
                r = focus_fn(parser.get_structure("d", dst)[0],
                             ["A"], ["C"], ["B"])
                rec["peptide_focus"] = round(r.peptide_focus, 3)
                rec["buried_peptide"] = round(r.buried_peptide, 1)
            except Exception as exc:
                rec["peptide_focus"] = ""
                rec["focus_error"] = str(exc)[:40]

        rec["accepted"] = 1
        rows.append(rec)

    cols = ["scaffold", "accepted", "source_target", "source_allele",
            "source_peptide", "binder_len", "chainB_len", "rmsd",
            "offset", "bb_clashes", "pep_clashes", "mhc_clashes", "pep_contacts",
            "mhc_contacts", "pep_min", "peptide_focus", "buried_peptide",
            "reason"]
    man = out_dir / "manifest.tsv"
    with man.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    ok = [r for r in rows if r["accepted"]]
    print(f"\ndocked {len(ok)} of {len(rows)} scaffolds -> {out_dir}")
    print(f"manifest: {man}")
    if len(ok) < len(rows):
        from collections import Counter
        print("\nrejected:")
        for reason, c in Counter(r["reason"].split("--")[0].strip()
                                 for r in rows
                                 if not r["accepted"]).most_common():
            print(f"  {c:>4}  {reason}")
    if ok:
        from collections import Counter
        print("\naccepted by source allele:")
        for a, c in Counter(r["source_allele"] for r in ok).most_common():
            print(f"  {c:>4}  {a}")
        lens = sorted(r["binder_len"] for r in ok)
        print(f"\nbinder length: min {lens[0]}  median "
              f"{lens[len(lens) // 2]}  max {lens[-1]}")
        foc = sorted(r["peptide_focus"] for r in ok
                     if isinstance(r.get("peptide_focus"), float))
        if foc:
            print(f"static peptide_focus: min {foc[0]:.3f}  median "
                  f"{foc[len(foc) // 2]:.3f}  max {foc[-1]:.3f}  "
                  f"(crystal cutoff 0.236)")
            print(f"  already above 0.236: {sum(1 for f in foc if f >= 0.236)}"
                  f" / {len(foc)}")


if __name__ == "__main__":
    main()
