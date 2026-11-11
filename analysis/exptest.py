"""NOTE: this script needs an input file that is not shipped in this repo
(structures/scaf26.pdb). It is kept for provenance; the path below shows where it
would live in a full working tree.
"""
from pymol import cmd
import os
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _repo import ROOT
B = str(ROOT / "structures")
cmd.load(os.path.join(B,"scaf26.pdb"),"cx")
cmd.hide("everything"); cmd.show("cartoon","chain A")
for ext in ("obj","dae","stl","wrl","gltf","glb"):
    p=os.path.join(B,"exptest."+ext)
    try:
        cmd.save(p)
        print(f"  {ext}: OK  {os.path.getsize(p) if os.path.exists(p) else 0} bytes")
    except Exception as e:
        print(f"  {ext}: FAILED  {type(e).__name__}: {str(e)[:60]}")
