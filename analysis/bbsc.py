"""Backbone vs side-chain split of binder-peptide contacts, per design.

Every arm-3 design faces the SAME peptide (SLLMWITQV), so the number of
available backbone and side-chain atoms is identical across designs and the
side-chain contact fraction is directly comparable -- no normalisation needed.

HYPOTHESIS UNDER TEST. The published PHOX2B design is described as making
"hydrogen bonds to the peptide backbone ... rigidifying the peptide". If that
cohort reads BACKBONE rather than side chains, promiscuous transfer is what you
would predict, because backbone is conserved across peptides and side chains
are not. That would mean our best-transferring cohort is the one least likely
to be peptide-SPECIFIC -- and peptide_focus cannot see the difference, since it
measures buried area, not which atoms are buried.
"""
import sys, csv, collections, random, statistics
import numpy as np

BB = {"N", "CA", "C", "O", "OXT"}
CONTACT = 4.5
PEP_LEN = 9
T = {"scrmsd": ("<", 2.0), "plddt_binder": (">", 92.0),
     "pae_interaction": ("<", 7.0), "peptide_focus": (">=", 0.202)}

def passes(r):
    try:
        for k, (op, c) in T.items():
            v = float(r[k])
            if op == "<" and not v < c: return False
            if op == ">" and not v > c: return False
            if op == ">=" and not v >= c: return False
        return True
    except (ValueError, KeyError):
        return False

def load(path):
    """chain A heavy atoms; chain B's last PEP_LEN residues' heavy atoms."""
    a, b = [], collections.OrderedDict()
    for l in open(path):
        if not l.startswith("ATOM"): continue
        if l[76:78].strip() == "H" or l[12:16].strip().startswith("H"): continue
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        if l[21] == "A":
            a.append(xyz)
        elif l[21] == "B":
            b.setdefault((l[22:26], l[26]), []).append((l[12:16].strip(), xyz))
    keys = list(b)[-PEP_LEN:]
    pep = [(i + 1, nm, xyz) for i, k in enumerate(keys) for nm, xyz in b[k]]
    return np.array(a), pep

man = {r["scaffold"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
rows = [r for r in csv.DictReader(open(sys.argv[2]), delimiter="\t") if passes(r)]
print("survivors: %d\n" % len(rows))

recs, per_pos = [], collections.defaultdict(lambda: [0, 0])
for r in rows:
    try:
        A, pep = load("%s/%s.pdb" % (sys.argv[3], r["description"]))
    except OSError:
        continue
    P = np.array([x for _, _, x in pep])
    d = np.linalg.norm(A[:, None] - P[None], axis=-1)
    hit = (d < CONTACT).sum(axis=0)          # contacts per peptide atom
    bb = sc = 0
    for (pos, nm, _), n in zip(pep, hit):
        if nm in BB:
            bb += n; per_pos[pos][0] += n
        else:
            sc += n; per_pos[pos][1] += n
    if bb + sc == 0: continue
    scaf = r["backbone_id"].replace("liu", "", 1).rsplit("_", 1)[0]
    m = man.get(scaf, {})
    recs.append(dict(d=r["description"], bb=bb, sc=sc, f=sc / (bb + sc),
                     tgt=m.get("source_target", "?"),
                     allele=m.get("source_allele", "?")))

print("=== side-chain fraction of peptide contacts, by source cohort ===")
print("  %-26s %4s %9s %9s %10s" % ("cohort", "n", "med s.c.", "mean s.c.", "mean contacts"))
g = collections.defaultdict(list)
for x in recs: g[x["tgt"]].append(x)
for k, v in sorted(g.items(), key=lambda kv: -len(kv[1])):
    print("  %-26s %4d %8.3f %9.3f %10.1f"
          % (k, len(v), statistics.median(x["f"] for x in v),
             statistics.mean(x["f"] for x in v),
             statistics.mean(x["bb"] + x["sc"] for x in v)))

ph = [x["f"] for x in recs if x["tgt"] == "PHOX2B"]
ot = [x["f"] for x in recs if x["tgt"] != "PHOX2B"]
print()
if ph and ot:
    obs = statistics.mean(ph) - statistics.mean(ot)
    pool = ph + ot
    random.seed(0)
    hits = 0
    for _ in range(20000):
        random.shuffle(pool)
        if abs(statistics.mean(pool[:len(ph)]) - statistics.mean(pool[len(ph):])) >= abs(obs):
            hits += 1
    print("  PHOX2B %.3f (n=%d) vs rest %.3f (n=%d)   diff %+.3f   permutation p=%.4f"
          % (statistics.mean(ph), len(ph), statistics.mean(ot), len(ot), obs, hits / 20000))
print()
print("=== where on the peptide the contacts land (all survivors pooled) ===")
seq = "SLLMWITQV"
print("  %-6s %10s %10s %8s" % ("pos", "backbone", "side chain", "s.c. %"))
for p in sorted(per_pos):
    bb, sc = per_pos[p]
    tot = bb + sc
    print("  %s%-5d %10d %10d %7.0f%%" % (seq[p-1], p, bb, sc, 100*sc/tot if tot else 0))
