import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _repo import ROOT
"""NOTE: this script needs an input file that is not shipped in this repo
(targets/transfer_2guo_kh46.pdb). It is kept for provenance; the path below shows where it
would live in a full working tree.
"""
"""Which atoms clash, and do they survive a backbone-only view?

RFdiffusion represents residues as backbone frames (N, CA, C) plus O; it does
not see side chains at all. So a side-chain clash in the placed complex is not
something diffusion has to resolve -- it is something ProteinMPNN repacks away
for free. Only a BACKBONE clash is a real geometric problem.
"""
import numpy as np
from Bio.PDB import PDBParser

BB = {"N", "CA", "C", "O"}
m = PDBParser(QUIET=True).get_structure("x", str(ROOT / "targets") + "/transfer_2guo_kh46.pdb")[0]

def atoms(chain, bb_only):
    return [(r.id[1], r.get_resname(), a.get_name(), a.coord)
            for r in m[chain] if r.id[0] == " " for a in r
            if a.element != "H" and (not bb_only or a.get_name() in BB)]

for bb_only in (False, True):
    tag = "backbone only" if bb_only else "all heavy atoms"
    print(f"\n=== {tag} ===")
    A = atoms("A", bb_only)
    for ch, label in (("B", "MHC"), ("C", "peptide")):
        T = atoms(ch, bb_only)
        a = np.array([x[3] for x in A]); t = np.array([x[3] for x in T])
        d = np.linalg.norm(a[:, None] - t[None], axis=-1)
        nclash = int((d < 2.5).sum())
        print(f"  KH46/{label:<8} min {d.min():>5.2f}  clashes {nclash:>3}  "
              f"contacts {int((d<4.5).sum()):>4}")
        if nclash and not bb_only:
            ii, jj = np.where(d < 2.5)
            for i, j in list(zip(ii, jj))[:14]:
                print(f"      {A[i][1]}{A[i][0]}:{A[i][2]:<3} -- "
                      f"{T[j][1]}{T[j][0]}:{T[j][2]:<3}  {d[i,j]:.2f}")

# per-peptide-residue contact census, old vs new
print("\n=== peptide contacts per position ===")
ref = PDBParser(QUIET=True).get_structure("r", str(ROOT / "targets") + "/9min_control_kh46.pdb")[0]
def census(model, bchain, pchain):
    b = np.array([a.coord for r in model[bchain] if r.id[0]==" " for a in r if a.element!="H"])
    out = []
    for r in model[pchain]:
        if r.id[0] != " ": continue
        p = np.array([a.coord for a in r if a.element != "H"])
        d = np.linalg.norm(b[:, None] - p[None], axis=-1)
        out.append((r.get_resname(), int((d < 4.5).sum()), float(d.min())))
    return out
old = census(ref, "D", "C")
new = census(m, "A", "C")
print(f"  {'pos':<4} {'old':<16} {'contacts':<9}   {'new':<16} {'contacts'}")
for i, (o, n) in enumerate(zip(old, new), 1):
    print(f"  {i:<4} {o[0]:<6} {o[2]:>5.2f}A {o[1]:>4}        "
          f"{n[0]:<6} {n[2]:>5.2f}A {n[1]:>4}")
