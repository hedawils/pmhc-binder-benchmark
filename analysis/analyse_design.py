#!/usr/bin/env python3
"""Contact and interaction analysis for one designed pMHC complex.

    usage:  analyse_design.py <design.pdb> [--peptide SLLMWITQV] [--json out.json]

Written to run LOCALLY on an AF2 output, so figure work never depends on the
GPU box staying reachable -- which it has not, three times in twelve hours.

Input layout is to_two_chain.py's: chain A = binder, chain B = MHC platform
followed by the peptide as the trailing residues.

WHAT IT REPORTS, AND THE LIMITS OF EACH

  contacts      binder heavy atoms within 4.5 A of a peptide SIDE-CHAIN atom.
                Backbone-only proximity is excluded, because a binder lying
                across the groove without reading anything would otherwise
                score contacts.

  H-bond        N/O ... N/O within 3.5 A. There are no hydrogens in an AF2
                model, so donor-H...acceptor angles cannot be computed: these
                are distance criteria, i.e. "close enough and chemically
                capable", not validated bonds.

  salt bridge   charged carboxylate/ammonium pairs within 4.0 A. For
                SLLMWITQV this must return zero -- the peptide carries no
                charged side chain -- which makes it a free correctness check
                on the classifier.

  pi-stacking   aromatic ring centroids within 5.5 A.

  CAVEAT THAT MATTERS: AF2 side-chain rotamers are the least reliable part of
  the model, and run_arm.sh uses relax_cycles=0, so nothing is energy
  minimised. A sub-2.5 A polar contact is a clash, not a strong bond, and is
  flagged as such. Treat the PATTERN as the result and individual distances
  as unreliable.
"""

import argparse
import collections
import json
import sys

import numpy as np

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
BACKBONE = {"N", "CA", "C", "O"}
ANION = {("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")}
CATION = {("LYS", "NZ"), ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
          ("HIS", "ND1"), ("HIS", "NE2")}
RINGS = {"TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
         "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
         "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
         "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"]}

CONTACT_CUT, HBOND_CUT, SALT_CUT, PI_CUT, CLASH = 4.5, 3.5, 4.0, 5.5, 2.5


def read(path):
    atoms = []
    for l in open(path):
        if not l.startswith("ATOM") or l[16] not in (" ", "A"):
            continue
        atoms.append(dict(
            ch=l[21], num=int(l[22:26]), name=l[12:16].strip(),
            res=l[17:20].strip(), el=(l[76:78].strip() or l[12:16].strip()[0]),
            xyz=np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return atoms


def centroid(atoms, ch, num):
    g = [a for a in atoms if a["ch"] == ch and a["num"] == num]
    if not g:
        return None, None
    rn = g[0]["res"]
    if rn not in RINGS:
        return None, rn
    pts = [a["xyz"] for a in g if a["name"] in RINGS[rn]]
    return (np.mean(pts, axis=0) if len(pts) >= 5 else None), rn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb")
    ap.add_argument("--peptide", default="SLLMWITQV")
    ap.add_argument("--json")
    a = ap.parse_args()
    pep_len = len(a.peptide)

    atoms = read(a.pdb)
    if not atoms:
        sys.exit(f"no ATOM records in {a.pdb}")
    A = [x for x in atoms if x["ch"] == "A"]
    Ball = [x for x in atoms if x["ch"] == "B"]
    if not A or not Ball:
        sys.exit("expected chain A (binder) and chain B (MHC + peptide)")

    nums = sorted({x["num"] for x in Ball})
    pep_nums = nums[-pep_len:]
    seq = "".join(THREE_TO_ONE.get(
        next(x["res"] for x in Ball if x["num"] == n), "X") for n in pep_nums)
    if seq != a.peptide:
        sys.exit(f"trailing {pep_len} residues of chain B read {seq}, "
                 f"expected {a.peptide} -- wrong layout or wrong peptide")

    plat = [n for n in nums if n not in pep_nums]
    print(f"{a.pdb}")
    print(f"  binder   chain A {min(x['num'] for x in A)}-"
          f"{max(x['num'] for x in A)}  ({len({x['num'] for x in A})} res)")
    print(f"  platform chain B {plat[0]}-{plat[-1]}  ({len(plat)} res)")
    print(f"  peptide  chain B {pep_nums[0]}-{pep_nums[-1]}  {seq}\n")

    out = {"design": a.pdb, "peptide": seq, "positions": {}}
    total, polar_all = 0, []
    for i, pn in enumerate(pep_nums, start=1):
        side = [x for x in Ball if x["num"] == pn and x["name"] not in BACKBONE]
        if not side:
            out["positions"][i] = {"residue": seq[i-1], "contacts": 0,
                                   "partners": [], "polar": []}
            print(f"  P{i} {seq[i-1]} ({'GLY' if seq[i-1]=='G' else '?'})"
                  f"   contacts=  0   (no side chain)")
            continue
        rn = side[0]["res"]
        partners, n_contact, polar = set(), 0, []
        for b in A:
            hit = False
            for p in side:
                d = float(np.linalg.norm(b["xyz"] - p["xyz"]))
                if d > CONTACT_CUT:
                    continue
                hit = True
                if b["el"] in "NO" and p["el"] in "NO" and d <= HBOND_CUT:
                    kind = ("salt bridge"
                            if (((p["res"], p["name"]) in ANION
                                 and (b["res"], b["name"]) in CATION)
                                or ((p["res"], p["name"]) in CATION
                                    and (b["res"], b["name"]) in ANION))
                            else "H-bond")
                    polar.append(dict(kind=kind, pep_atom=p["name"],
                                      binder=f"{THREE_TO_ONE.get(b['res'],'X')}"
                                             f"{b['num']}:{b['name']}",
                                      dist=round(d, 2),
                                      clash=d < CLASH))
            if hit:
                partners.add((b["num"], b["res"]))
                n_contact += 1
        total += n_contact
        polar_all += polar
        out["positions"][i] = {
            "residue": seq[i-1], "resname": rn, "contacts": n_contact,
            "partners": sorted(p[0] for p in partners),
            "polar": polar}
        bar = "#" * min(40, n_contact // 2)
        print(f"  P{i} {seq[i-1]} ({rn})  contacts={n_contact:3d}  {bar}")
        for q in sorted(polar, key=lambda z: z["dist"]):
            flag = "   <-- TOO SHORT, unrelaxed geometry" if q["clash"] else ""
            print(f"       {q['kind']:12s} {q['pep_atom']:4s} ... "
                  f"{q['binder']:12s} {q['dist']:.2f} A{flag}")

    # pi-stacking between peptide aromatics and binder aromatics
    pi = []
    for i, pn in enumerate(pep_nums, start=1):
        cp, rp = centroid(atoms, "B", pn)
        if cp is None:
            continue
        for bn in sorted({x["num"] for x in A}):
            cb, rb = centroid(atoms, "A", bn)
            if cb is None:
                continue
            d = float(np.linalg.norm(cp - cb))
            if d <= PI_CUT:
                pi.append(dict(position=i, pep=rp, binder=f"{rb}{bn}",
                               dist=round(d, 2)))
    out["pi_stacking"] = pi

    read_pos = [i for i, v in out["positions"].items() if v["contacts"] > 0]
    zero_pos = [i for i, v in out["positions"].items() if v["contacts"] == 0]
    on_read = sum(out["positions"][i]["contacts"] for i in read_pos)
    print(f"\n  total side-chain contacts : {total}")
    print(f"  positions contacted       : {read_pos}")
    print(f"  positions with ZERO       : {zero_pos}")
    if total:
        print(f"  concentration             : {on_read}/{total} "
              f"({100*on_read/total:.0f}%) on {len(read_pos)} positions")
    print(f"  salt bridges              : "
          f"{sum(1 for q in polar_all if q['kind']=='salt bridge')} "
          f"(expect 0 for a peptide with no charged side chain)")
    for q in pi:
        print(f"  pi-stacking               : P{q['position']} {q['pep']} "
              f"... {q['binder']}  {q['dist']} A")
    n_clash = sum(1 for q in polar_all if q["clash"])
    if n_clash:
        print(f"  WARNING: {n_clash} polar contact(s) below {CLASH} A -- "
              f"unrelaxed AF2 geometry, distances not quantitative")

    out["total_contacts"] = total
    if a.json:
        json.dump(out, open(a.json, "w"), indent=2)
        print(f"\nwrote {a.json}")


if __name__ == "__main__":
    main()
