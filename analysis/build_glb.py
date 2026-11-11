"""Assemble the per-component OBJs into one PowerPoint-ready .glb.

Run:  glbenv/bin/python build_glb.py

Why a GLB and not the OBJ PyMOL writes: PyMOL's OBJ export carries geometry
only -- no .mtl, no vertex colours -- so PowerPoint would render the whole
complex as one flat grey mass. GLB embeds materials in a single file, which is
also the format PowerPoint handles most reliably.

ORIENTATION MATTERS FOR THE TURNTABLE

PowerPoint's Turntable animation spins the model about its vertical (Y) axis.
So the model is rotated here such that:

    +Y   runs from the peptide up through the binder (the stacking direction)
    +X   runs along the peptide, i.e. down the groove

That way a turntable sweep keeps the binder sitting on the groove throughout
and shows the interface from every side, instead of tumbling end over end.

NO BAKED TEXT LABELS, DELIBERATELY

alpha1 and alpha2 are distinguished by colour, not by 3D text. Text in a model
file is geometry: it would be mirrored and unreadable through half of every
rotation. The legend belongs in PowerPoint, as static text beside the model.
"""

import os
import numpy as np
import trimesh
from trimesh.visual.material import PBRMaterial

B = os.path.dirname(os.path.abspath(__file__))
PARTS_DIR = os.path.join(B, "parts")
PDB = os.path.join(B, "structures", "scaf128.pdb")
OUT = os.path.join(B, "scaffold128_3d.glb")

# name -> (hex colour, what it is)
COLOURS = {
    "platform":   ("#D6D5CE", "HLA platform (beta sheet floor)"),
    "helix_a1":   ("#5F6B73", "alpha-1 groove helix"),
    "helix_a2":   ("#B3A994", "alpha-2 groove helix"),
    "binder":     ("#2D6A9F", "de novo binder, 79 aa"),
    "peptide":    ("#16202B", "SLLMWITQV backbone"),
    "pep_read":   ("#C45026", "peptide positions the binder contacts"),
    "pep_anchor": ("#E0B01B", "buried anchors, zero contacts"),
    "partners":   ("#33A07F", "binder residues making polar contacts"),
}


def hex_rgba(h, a=1.0):
    h = h.lstrip("#")
    return [int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4)] + [a]


def frame_from_structure():
    """Build a rotation putting the groove along X and the stack along Y."""
    ca = {}
    for line in open(PDB):
        if line.startswith("ATOM") and line[12:16] == " CA ":
            ca.setdefault(line[21], {})[int(line[22:26])] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    pep = np.array([v for k, v in sorted(ca["B"].items()) if k >= 260])
    mhc = np.array([v for k, v in sorted(ca["B"].items()) if k < 260])
    binder = np.array(list(ca["A"].values()))

    # X: the peptide's long axis, from its first principal component
    centred = pep - pep.mean(0)
    x = np.linalg.svd(centred, full_matrices=False)[2][0]
    if np.dot(pep[-1] - pep[0], x) < 0:
        x = -x
    # Y: peptide/MHC -> binder, i.e. the stacking direction, orthogonalised
    y = binder.mean(0) - mhc.mean(0)
    y = y - np.dot(y, x) * x
    y /= np.linalg.norm(y)
    z = np.cross(x, y)
    R = np.eye(4)
    R[:3, :3] = np.vstack([x, y, z])        # world -> model frame
    return R, pep.mean(0)


def main():
    R, pivot = frame_from_structure()
    T = np.eye(4)
    T[:3, 3] = -pivot                        # centre on the peptide
    M = R @ T

    scene = trimesh.Scene()
    total_f = 0
    for name, (hexcol, what) in COLOURS.items():
        path = os.path.join(PARTS_DIR, f"{name}.obj")
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            print(f"  {name:11s} MISSING, skipped")
            continue
        mesh = trimesh.load(path, force="mesh", process=False)
        if mesh.faces.shape[0] == 0:
            print(f"  {name:11s} no faces, skipped")
            continue
        mesh.apply_transform(M)
        mesh.visual = trimesh.visual.TextureVisuals(
            material=PBRMaterial(
                name=name,
                baseColorFactor=hex_rgba(hexcol),
                metallicFactor=0.0,
                roughnessFactor=0.62,
                doubleSided=True,
            ))
        scene.add_geometry(mesh, geom_name=name)
        total_f += mesh.faces.shape[0]
        print(f"  {name:11s} {mesh.faces.shape[0]:6d} faces  {hexcol}  {what}")

    scene.export(OUT)
    print(f"\nwrote {OUT}")
    print(f"  {total_f:,} triangles, {os.path.getsize(OUT)/1e6:.1f} MB")
    ext = scene.bounds[1] - scene.bounds[0]
    print(f"  extent X {ext[0]:.0f} A (groove)  Y {ext[1]:.0f} A (stack)  "
          f"Z {ext[2]:.0f} A")


if __name__ == "__main__":
    main()
