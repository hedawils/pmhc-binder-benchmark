#!/usr/bin/env python3
"""Regenerate every number in the README from data/, and nothing else.

    python analysis/reproduce.py              # all sections
    python analysis/reproduce.py --section 3  # just one

Standard library only, on purpose: the point of this script is that a reader
can check the claims without installing anything or trusting a notebook.

Paths resolve from the repository root, found relative to this file, so it
runs from any working directory.
"""

import argparse
import csv
import math
import re
import statistics as st
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

# The four gates. A design is a hit only if it clears all of them.
GATES = dict(scrmsd=2.0, plddt_binder=92.0, pae_interaction=7.0, peptide_focus=0.202)

ARMS = [
    # key,        label,                  target
    ("liu",      "NY-ESO-1 repurposed",   "NY-ESO-1"),
    ("free",     "NY-ESO-1 free",         "NY-ESO-1"),
    ("kh46self", "NY-ESO-1 KH46 control", "NY-ESO-1"),
    ("kf9",      "KF9 repurposed",        "KF9"),
    ("kf9free",  "KF9 free",              "KF9"),
]


# --------------------------------------------------------------------------
# loading

def rows(key):
    with open(DATA / f"{key}_metrics.tsv") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return None


def is_hit(r):
    v = {k: num(r, k) for k in GATES}
    if any(x is None for x in v.values()):
        return False
    return (v["scrmsd"] < GATES["scrmsd"]
            and v["plddt_binder"] > GATES["plddt_binder"]
            and v["pae_interaction"] < GATES["pae_interaction"]
            and v["peptide_focus"] >= GATES["peptide_focus"])


def scaffold_of(backbone_id):
    """liuscaf100_0 / kf9scaf100_0 -> scaf100; free arms have no scaffold prefix."""
    m = re.search(r"(scaf\d+)", backbone_id)
    return m.group(1) if m else backbone_id


# --------------------------------------------------------------------------
# statistics, implemented here so the script stays dependency-free

def _ranks(v):
    idx = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    ties = []
    i = 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and v[idx[j + 1]] == v[idx[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[idx[k]] = avg
        ties.append(j - i + 1)
        i = j + 1
    return r, ties


def mannwhitney(a, b):
    """Two-sided, tie-corrected normal approximation with continuity correction."""
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return None
    r, ties = _ranks(list(a) + list(b))
    u = sum(r[:n1]) - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    n = n1 + n2
    if n < 2:
        return None
    var = n1 * n2 / 12.0 * ((n + 1) - sum(t ** 3 - t for t in ties) / (n * (n - 1)))
    if var <= 0:
        return None
    z = (abs(u - mu) - 0.5) / math.sqrt(var)
    return math.erfc(z / math.sqrt(2.0))


def fisher(a, b, c, d):
    """Two-sided Fisher exact on a 2x2 table, by summing tables no more likely."""
    def lc(n):
        return math.lgamma(n + 1)

    def p(a_, b_, c_, d_):
        return math.exp(lc(a_ + b_) + lc(c_ + d_) + lc(a_ + c_) + lc(b_ + d_)
                        - lc(a_) - lc(b_) - lc(c_) - lc(d_) - lc(a_ + b_ + c_ + d_))

    p0 = p(a, b, c, d)
    tot = 0.0
    for i in range(0, a + b + 1):
        j, k, l = a + b - i, a + c - i, d - a + i
        if k < 0 or l < 0:
            continue
        pi = p(i, j, k, l)
        if pi <= p0 * (1 + 1e-9):
            tot += pi
    return min(tot, 1.0)


def fmt_p(p):
    if p is None:
        return "n/a"
    return "ns" if p >= 0.05 else (f"{p:.1e}" if p < 1e-4 else f"{p:.4f}")


# --------------------------------------------------------------------------

def capped_per_design(key, cap):
    """Per-design hits using at most `cap` sequences per scaffold.

    Equalises sampling effort between arms. Sequences are taken in seq_id
    order rather than sampled, so the result does not depend on a seed.
    """
    by = defaultdict(list)
    for r in rows(key):
        by[r["backbone_id"]].append(r)
    n = h = 0
    for rs in by.values():
        rs = sorted(rs, key=lambda r: r["seq_id"])[:cap]
        n += len(rs)
        h += sum(1 for r in rs if is_hit(r))
    return h, n


def section_1_hit_rates():
    """Headline table: per-design and per-scaffold hit rates."""
    print("\n" + "=" * 78)
    print("1. HIT RATES")
    print("=" * 78)
    print("\n   A design is a hit when it clears all four gates:")
    print("   scRMSD < 2.0 | binder pLDDT > 92 | interface pAE < 7.0 | peptide_focus >= 0.202\n")
    print(f"   {'arm':24} {'per design':>18} {'per scaffold':>16} {'4-seq basis':>14}")

    store = {}
    for key, label, _ in ARMS:
        rs = rows(key)
        hits = [r for r in rs if is_hit(r)]
        store[key] = (len(hits), len(rs))

        tot, win = defaultdict(int), defaultdict(int)
        for r in rs:
            s = scaffold_of(r["backbone_id"])
            tot[s] += 1
            if is_hit(r):
                win[s] += 1
        raw = sum(1 for s in tot if win[s] > 0)

        # Arms sampled 6 sequences per scaffold are not comparable to arms that
        # sampled 4: more draws mechanically raise the chance a scaffold yields
        # at least one hit. Rather than subsample once (which depends on a seed),
        # take the expectation over all 4-subsets:
        #     P(scaffold yields >=1 hit | 4 drawn) = 1 - C(n-k, 4) / C(n, 4)
        exp4 = 0.0
        for s, n in tot.items():
            k = win[s]
            if n >= 4:
                exp4 += 1 - math.comb(n - k, 4) / math.comb(n, 4)
            else:
                exp4 += 1.0 if k else 0.0

        print(f"   {label:24} {len(hits):5}/{len(rs):<5} {100*len(hits)/len(rs):5.2f}%"
              f"   {raw:4}/{len(tot):<4} {100*raw/len(tot):5.1f}%"
              f"   {100*exp4/len(tot):8.1f}%")

    hl, nl = store["liu"]
    hf, nf = store["free"]
    hk, nk = store["kf9"]
    hkf, nkf = store["kf9free"]

    # The NY-ESO-1 repurposed arm sampled 6 sequences per scaffold; every other
    # arm sampled 4. Comparing all-sequence rates therefore compares unequal
    # sampling effort, so the headline ratio caps every arm at 4 sequences per
    # scaffold. The uncapped ratio is reported too, since it uses all the data.
    cl, dl = capped_per_design("liu", 4)
    cf, df = capped_per_design("free", 4)
    print(f"\n   NY-ESO-1, capped at 4 seqs/scaffold: "
          f"{cl}/{dl} = {100*cl/dl:.2f}%  vs  {cf}/{df} = {100*cf/df:.2f}%")
    print(f"   NY-ESO-1  repurposed / free = {(cl/dl)/(cf/df):.2f}x  (capped, headline)   "
          f"Fisher p = {fmt_p(fisher(cl, dl-cl, cf, df-cf))}")
    print(f"                               = {(hl/nl)/(hf/nf):.2f}x  (all sequences)      "
          f"Fisher p = {fmt_p(fisher(hl, nl-hl, hf, nf-hf))}")
    print(f"   KF9       repurposed / free = {(hk/nk)/(hkf/nkf):.2f}x  (both arms are "
          f"already 4/scaffold)   Fisher p = {fmt_p(fisher(hk, nk-hk, hkf, nkf-hkf))}")


def section_2_fold_vs_interface():
    """Why it works when it works: fold penalty, interface gain."""
    print("\n" + "=" * 78)
    print("2. FOLD vs INTERFACE  (all designs, medians, Mann-Whitney)")
    print("=" * 78)

    for target, rep, free in (("NY-ESO-1", "liu", "free"), ("KF9", "kf9", "kf9free")):
        P, F = rows(rep), rows(free)
        print(f"\n   {target}   repurposed n={len(P)}  free n={len(F)}")

        for label, col in (("binder pLDDT", "plddt_binder"), ("scRMSD", "scrmsd")):
            a = [num(r, col) for r in P if num(r, col) is not None]
            b = [num(r, col) for r in F if num(r, col) is not None]
            print(f"     {label:28} {st.median(a):8.2f} {st.median(b):8.2f}   "
                  f"p = {fmt_p(mannwhitney(a, b))}")

        # Unconditional peptide-burial medians are 0 in every arm, because most
        # designs never reach the peptide. Report the contact rate, then the
        # geometry among designs that do make contact.
        def touching(rs):
            return [r for r in rs if (num(r, "buried_peptide") or 0) > 0]

        pc, fc = touching(P), touching(F)
        print(f"     {'designs touching peptide':28} {100*len(pc)/len(P):7.1f}% "
              f"{100*len(fc)/len(F):7.1f}%   "
              f"p = {fmt_p(fisher(len(pc), len(P)-len(pc), len(fc), len(F)-len(fc)))}")

        for label, col in (("buried peptide A^2 (contact)", "buried_peptide"),
                           ("interface pAE (contact)", "pae_interaction")):
            a = [num(r, col) for r in pc if num(r, col) is not None]
            b = [num(r, col) for r in fc if num(r, col) is not None]
            print(f"     {label:28} {st.median(a):8.2f} {st.median(b):8.2f}   "
                  f"p = {fmt_p(mannwhitney(a, b))}")


def section_3_interface_composition():
    """How much of the interface is peptide rather than MHC platform."""
    print("\n" + "=" * 78)
    print("3. INTERFACE COMPOSITION  (peptide share of buried surface)")
    print("=" * 78)
    print(f"\n   {'arm':24} {'n':>5} {'median':>9} {'IQR':>14} {'max':>8}")

    for key, label, _ in ARMS:
        frac = []
        for r in rows(key):
            if not is_hit(r):
                continue
            p, m = num(r, "buried_peptide"), num(r, "buried_mhc")
            if p is not None and m is not None and p + m > 0:
                frac.append(100 * p / (p + m))
        if not frac:
            continue
        s = sorted(frac)
        print(f"   {label:24} {len(s):5} {st.median(s):8.1f}% "
              f"{s[len(s)//4]:6.1f}-{s[3*len(s)//4]:<6.1f} {max(s):7.1f}%")

    # Peptide-dominant designs are generated; they just do not survive the gates.
    allfrac = []
    for key in ("liu", "free"):
        for r in rows(key):
            p, m = num(r, "buried_peptide"), num(r, "buried_mhc")
            if p is not None and m is not None and p + m > 100:
                allfrac.append(100 * p / (p + m))
    allfrac.sort()
    print(f"\n   all NY-ESO-1 designs with a real interface (n={len(allfrac)}): "
          f"median {st.median(allfrac):.1f}%, "
          f"95th pct {allfrac[int(.95*len(allfrac))]:.1f}%, max {max(allfrac):.1f}%")


def _specificity(key):
    """design -> condition -> position -> holo-minus-apo delta."""
    per = defaultdict(lambda: defaultdict(dict))
    with open(DATA / f"{key}_specificity.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                per[r["design"]][r["condition"]][int(r["position"])] = float(r["delta"])
            except (TypeError, ValueError):
                pass
    return per


def section_4_positions():
    """Which peptide positions the binder actually reads."""
    print("\n" + "=" * 78)
    print("4. WHICH POSITIONS DOES THE BINDER READ?")
    print("=" * 78)
    print("\n   Mean holo-minus-apo per position. Large values mean the binder changes")
    print("   the groove's own preference there; ~0 means it does not reach it.\n")

    peptide = "SLLMWITQV"
    print("   " + " " * 13 + " ".join(f"{'P%d' % p:>6}" for p in range(1, 10)))
    print("   " + " " * 13 + " ".join(f"{c:>6}" for c in peptide))
    for key, label in (("liu", "repurposed"), ("free", "free"), ("kh46self", "KH46")):
        pos = defaultdict(list)
        for _, cond in _specificity(key).items():
            for p, v in cond.get("holo", {}).items():
                pos[p].append(v)
        print(f"   {label:13}" + " ".join(f"{st.mean(pos[p]):+6.2f}" for p in sorted(pos)))

    print("\n   P2 and P9 are the canonical HLA-A*02:01 anchors, buried in the groove.")
    print("   Near-zero there is the holo-apo subtraction working as intended.")


def section_5_specificity_margin():
    """Does each design prefer its own peptide over the decoys?"""
    print("\n" + "=" * 78)
    print("5. SPECIFICITY MARGIN BY ARM")
    print("=" * 78)
    print("\n   margin = summed delta for the on-target peptide")
    print("            minus the best decoy's summed delta\n")
    print(f"   {'arm':24} {'n':>5} {'median':>10} {'prefers on-target':>20}")

    store = {}
    for key, label in (("free", "NY-ESO-1 free"), ("liu", "NY-ESO-1 repurposed"),
                       ("kh46self", "NY-ESO-1 KH46 control")):
        margins = []
        for _, cond in _specificity(key).items():
            if "holo" not in cond:
                continue
            on = sum(cond["holo"].values())
            dec = [sum(v.values()) for k, v in cond.items() if k.startswith("decoy")]
            if dec:
                margins.append(on - max(dec))
        store[key] = margins
        print(f"   {label:24} {len(margins):5} {st.median(margins):+9.2f} "
              f"{100*sum(1 for m in margins if m > 0)/len(margins):18.1f}%")

    print(f"\n   repurposed vs free: p = {fmt_p(mannwhitney(store['liu'], store['free']))}")
    print("   The control arm has the highest hit rate and the worst specificity.")
    print("   Not a size effect: KH46 designs are 128 aa, free-diffusion 130 aa.")


def section_6_donors():
    """Hit rate by the target each scaffold was originally designed against."""
    print("\n" + "=" * 78)
    print("6. DONOR HETEROGENEITY")
    print("=" * 78)
    print("\n   NOTE: every HLA-C donor comes from one source target (PHOX2B), so")
    print("   donor identity and HLA locus cannot be separated in this design.\n")

    res = defaultdict(dict)
    meta = {}
    for key in ("liu", "kf9"):
        rs = rows(key)
        tot, win = defaultdict(int), defaultdict(int)
        for r in rs:
            s = scaffold_of(r["backbone_id"])
            tot[s] += 1
            if is_hit(r):
                win[s] += 1
        with open(DATA / f"manifest_{key}.tsv") as fh:
            for m in csv.DictReader(fh, delimiter="\t"):
                if m["accepted"] != "1" or m["scaffold"] not in tot:
                    continue
                src = m["source_target"]
                meta[src] = m["source_allele"] or "-"
                a, b = res[src].get(key, (0, 0))
                res[src][key] = (a + 1, b + (1 if win[m["scaffold"]] else 0))

    print(f"   {'donor target':16} {'allele':10} {'NY-ESO-1':>16} {'KF9':>16}")
    def rate(t):
        n, w = t
        return f"{w}/{n} = {100*w/n:5.1f}%" if n else "-"
    for src in sorted(res, key=lambda s: -(res[s].get("liu", (1, 0))[1]
                                           / max(res[s].get("liu", (1, 0))[0], 1))):
        print(f"   {src:16} {meta[src]:10} {rate(res[src].get('liu', (0,0))):>16} "
              f"{rate(res[src].get('kf9', (0,0))):>16}")


SECTIONS = [section_1_hit_rates, section_2_fold_vs_interface,
            section_3_interface_composition, section_4_positions,
            section_5_specificity_margin, section_6_donors]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--section", type=int, choices=range(1, len(SECTIONS) + 1),
                    help="run one section instead of all")
    args = ap.parse_args()

    if not DATA.is_dir():
        raise SystemExit(f"data/ not found at {DATA} -- run from a full checkout")

    for i, fn in enumerate(SECTIONS, 1):
        if args.section in (None, i):
            fn()
    print()


if __name__ == "__main__":
    main()
