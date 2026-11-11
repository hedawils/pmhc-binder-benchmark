import csv, math, statistics as st

import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _repo import B, ROOT, sibling

FOCUS = 0.202

def load(key):
    rows = list(csv.DictReader(open(B + key + "_metrics.tsv"), delimiter="\t"))
    out = []
    for r in rows:
        d = {}
        for k in ("scrmsd","plddt_binder","pae_interaction","peptide_focus",
                  "buried_peptide","buried_mhc","plddt_target"):
            try: d[k] = float(r[k])
            except (TypeError, ValueError): d[k] = None
        out.append(d)
    return out

def survivors(rows):
    return [r for r in rows
            if r["scrmsd"] is not None and r["scrmsd"] < 2.0
            and r["plddt_binder"] is not None and r["plddt_binder"] > 92
            and r["pae_interaction"] is not None and r["pae_interaction"] < 7.0
            and r["peptide_focus"] is not None and r["peptide_focus"] >= FOCUS]

def ranks(v):
    idx = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0]*len(v); i = 0; ties = []
    while i < len(idx):
        j = i
        while j+1 < len(idx) and v[idx[j+1]] == v[idx[i]]: j += 1
        avg = (i + j)/2.0 + 1.0
        for k in range(i, j+1): r[idx[k]] = avg
        ties.append(j-i+1); i = j+1
    return r, ties

def mannwhitney(a, b):
    """two-sided, tie-corrected normal approximation; returns (U, p)"""
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0: return None, None
    r, ties = ranks(a + b)
    R1 = sum(r[:n1])
    U1 = R1 - n1*(n1+1)/2.0
    mu = n1*n2/2.0
    N = n1 + n2
    tie = sum(t**3 - t for t in ties)
    var = n1*n2/12.0 * ((N + 1) - tie/float(N*(N-1))) if N > 1 else 0.0
    if var <= 0: return U1, None
    z = (U1 - mu) / math.sqrt(var)
    # correct for continuity
    z = (abs(U1 - mu) - 0.5) / math.sqrt(var)
    p = math.erfc(z / math.sqrt(2.0))
    return U1, p

def med(rows, k):
    v = [r[k] for r in rows if r[k] is not None]
    return (st.median(v), v) if v else (None, [])

def stars(p):
    if p is None: return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"

METRICS = [
    ("FOLD",      "binder pLDDT",     "plddt_binder",    "higher better"),
    ("FOLD",      "scRMSD (Å)",       "scrmsd",          "lower better"),
    ("INTERFACE", "interface pAE",    "pae_interaction", "lower better"),
    ("INTERFACE", "peptide focus",    "peptide_focus",   "higher better"),
    ("INTERFACE", "buried peptide Å²","buried_peptide",  "higher better"),
    ("INTERFACE", "buried MHC Å²",    "buried_mhc",      "higher better"),
]

PAIRS = [("NY-ESO-1 / A*02:01", "liu", "free"),
         ("SARS-CoV-2 KF9 / C*12:02", "kf9", "kf9free")]

for pop_name, sel in (("ALL DESIGNS", lambda x: x), ("FOUR-GATE SURVIVORS", survivors)):
    print("\n" + "="*104)
    print(f"  {pop_name}")
    print("="*104)
    for tname, pkey, fkey in PAIRS:
        P, F = sel(load(pkey)), sel(load(fkey))
        print(f"\n{tname}   partial n={len(P)}   free n={len(F)}")
        print(f"  {'metric':20} {'partial':>12} {'free':>12} {'delta':>10} {'p':>12}  {'':4} favours")
        for grp, label, key, dirn in METRICS:
            mp, vp = med(P, key); mf, vf = med(F, key)
            if mp is None or mf is None: continue
            _, p = mannwhitney(vp, vf)
            d = mp - mf
            better = "partial" if ((d > 0) == ("higher" in dirn)) else "free"
            if p is not None and p >= 0.05: better = "—"
            ps = "n/a" if p is None else (f"{p:.2e}" if p < 1e-4 else f"{p:.4f}")
            print(f"  {grp[0]} {label:18} {mp:12.3f} {mf:12.3f} {d:+10.3f} {ps:>12}  {stars(p):4} {better}")
