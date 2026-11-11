"""Render a 360 degree turntable of the scaffold-128 interface for PowerPoint.

Run headless:  pymol -cq turntable.py

Light background and the deck palette, so it drops onto the light slides with
no dark box around it. Frames are ray-traced; the GIF is assembled separately.

WHY SCAFFOLD 128 AND NOT 26

Scaffold 26 had the best interface pAE of the 67 survivors, which is why it
was the first choice. It also FAILS its own specificity check: it scores the
MART-1 decoy above NY-ESO-1 and ranks 57/67 on specificity margin. Scaffold
128 is a hair worse on pAE (3.54 vs 3.31), better on peptide_focus (0.268 vs
0.226), and ranks FIRST of 67 on margin -- its on-target score is positive
while all four decoys are negative. It is the design that survives questions.

Chain map (to_two_chain.py layout, note this differs from scaffold 26's):
    A        1-79    de novo binder
    B       80-259   HLA-A*02:01 alpha1/alpha2 platform
      a1   128-170   groove helix
      a2   215-259   groove helix
    B      260-268   SLLMWITQV  (P1..P9)
"""

import os
from pymol import cmd

from _repo import ROOT
BASE = str(ROOT / "structures")
SRC = os.path.join(BASE, "structures", "scaf128.pdb")
OUT = os.path.join(BASE, "frames")

# 96 frames at 3.75 deg. Halving the speed by doubling the frame delay would
# expose the angular steps instead of hiding them, so the step is halved and
# the cadence kept at 14 fps.
FRAMES = 96
W, H = 640, 480

PLATFORM = "80-259"
HELICES = "128-170+215-259"
PEPTIDE = "260-268"
READ = "263+264+266+267"          # P4 P5 P7 P8 -- all 137 contacts land here
ANCHOR = "260+261+262+265+268"    # P1 P2 P3 P6 P9 -- zero contacts
PARTNERS = "39"                   # binder Asn39, the one polar contact

os.makedirs(OUT, exist_ok=True)
cmd.load(SRC, "cx")
cmd.remove("hydrogens")
cmd.bg_color("white")
cmd.hide("everything")
cmd.show("cartoon", "cx")

cmd.set_color("deepblue", [0.176, 0.416, 0.624])   # binder
cmd.set_color("palegrey", [0.867, 0.863, 0.839])   # platform body
cmd.set_color("midgrey",  [0.573, 0.569, 0.533])   # groove helices, neutral
cmd.set_color("rust",     [0.769, 0.310, 0.149])   # anchors: zero contacts
cmd.set_color("gold",     [0.878, 0.690, 0.106])   # positions actually read
cmd.set_color("ink",      [0.086, 0.125, 0.169])   # peptide backbone
cmd.set_color("teal",     [0.200, 0.620, 0.510])   # polar contact

cmd.color("deepblue", "chain A")
cmd.color("palegrey", f"chain B and resi {PLATFORM}")
cmd.color("midgrey",  f"chain B and resi {HELICES}")
cmd.color("ink",      f"chain B and resi {PEPTIDE}")

cmd.show("sticks", f"chain B and resi {PEPTIDE} and not (name C+N+O)")
cmd.set("cartoon_side_chain_helper", 1)
# Inverted on request: the CONTACTED positions take the strong warm colour so
# the eye lands on what the slide is arguing about; the untouched anchors take
# the softer gold.
cmd.color("rust", f"chain B and resi {READ} and sidechain")
cmd.color("gold", f"chain B and resi {ANCHOR} and sidechain")

cmd.show("sticks", f"chain A and resi {PARTNERS} and not (name C+N+O)")
cmd.color("teal", f"chain A and resi {PARTNERS} and sidechain")

# the binder is the bulk of the picture; translucent so the peptide it grips
# stays visible through it as the model turns
cmd.set("cartoon_transparency", 0.22, "chain A")
cmd.set("stick_radius", 0.23, f"chain B and resi {PEPTIDE}")

# Thr7 OG1 ... Asn39 ND2, 3.44 A -- a real H-bond distance, unlike scaffold
# 26's 2.22 A contact which was an unrelaxed-geometry artifact
cmd.distance("hb0", "chain B and resi 266 and name OG1",
             "chain A and resi 39 and name ND2")
cmd.color("teal", "hb0")
cmd.set("dash_width", 2.6)
cmd.set("dash_gap", 0.33)
cmd.hide("labels")

cmd.set("ray_opaque_background", 1)
cmd.set("antialias", 2)
cmd.set("ray_shadows", 0)
cmd.set("specular", 0.2)
cmd.set("ambient", 0.22)
cmd.set("direct", 0.55)
cmd.set("cartoon_fancy_helices", 1)
cmd.set("cartoon_highlight_color", "grey70")

cmd.orient("cx")
cmd.zoom("cx", 3.0)
cmd.turn("x", -12)

step = 360.0 / FRAMES
for i in range(FRAMES):
    cmd.png(os.path.join(OUT, f"f{i:03d}.png"), width=W, height=H,
            dpi=150, ray=1)
    cmd.turn("y", step)

print(f"wrote {FRAMES} frames to {OUT}")
