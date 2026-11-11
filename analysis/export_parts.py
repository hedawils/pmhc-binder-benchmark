"""Export each coloured component of the scaffold-128 complex as its own OBJ.

Run:  pymol -cq export_parts.py

PyMOL's OBJ writer emits geometry only -- no .mtl, no vertex colours -- so one
combined OBJ would arrive in PowerPoint as a single flat grey blob. Exporting
one file per component lets build_glb.py attach a known material to each,
which is what makes alpha1 and alpha2 distinguishable without baked text
labels (text geometry faces backwards through half of any turntable).

Side chains go out as SURFACES, not sticks: PyMOL ray-traces sticks and
spheres as primitives and its OBJ writer emits nothing for them -- every stick
export came back 0 bytes. Surfaces are real triangle meshes.

Chain map (scaffold 128; differs from scaffold 26's numbering):
    A        1-79    binder
    B       80-259   platform, with groove helices 128-170 and 215-259
    B      260-268   SLLMWITQV
"""

import os
from pymol import cmd

from _repo import ROOT
B = str(ROOT / "structures")
SRC = os.path.join(B, "structures", "scaf128.pdb")
OUT = os.path.join(B, "parts")
os.makedirs(OUT, exist_ok=True)

READ = "263+264+266+267"           # P4 P5 P7 P8
ANCHOR = "260+261+262+265+268"     # P1 P2 P3 P6 P9

cmd.load(SRC, "cx")
cmd.remove("hydrogens")
cmd.set("cartoon_sampling", 6)
cmd.set("surface_quality", 0)
cmd.set("solvent_radius", 1.2)
cmd.set("cartoon_fancy_helices", 1)

PARTS = [
    ("binder",     "chain A",                                              "cartoon"),
    ("platform",   "chain B and resi 80-259 and not resi 128-170+215-259", "cartoon"),
    ("helix_a1",   "chain B and resi 128-170",                             "cartoon"),
    ("helix_a2",   "chain B and resi 215-259",                             "cartoon"),
    ("peptide",    "chain B and resi 260-268",                             "cartoon"),
    ("pep_read",   f"chain B and resi {READ} and sidechain",               "surface"),
    ("pep_anchor", f"chain B and resi {ANCHOR} and sidechain",             "surface"),
    ("partners",   "chain A and resi 39 and sidechain",                    "surface"),
]

for name, sel, rep in PARTS:
    cmd.hide("everything")
    n = cmd.count_atoms(sel)
    if n == 0:
        print(f"  {name}: EMPTY selection, skipped")
        continue
    cmd.show(rep, sel)
    path = os.path.join(OUT, f"{name}.obj")
    cmd.save(path)
    size = os.path.getsize(path) if os.path.exists(path) else 0
    print(f"  {name:11s} {rep:8s} atoms={n:4d}  {size/1e6:.2f} MB")

print("done")
