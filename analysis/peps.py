import sys, os, collections
A3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
      "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
      "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
d = sys.argv[1]
rows = []
for fn in sorted(os.listdir(d)):
    if not fn.endswith(".pdb"): continue
    ch = collections.OrderedDict()
    with open(os.path.join(d, fn), errors="replace") as fh:
        for l in fh:
            if l.startswith("ATOM") and l[12:16].strip() == "CA":
                ch.setdefault(l[21], []).append(A3.get(l[17:20].strip(), "X"))
    # a class-I peptide is 8-11 residues; the MHC platform is ~180, binder 50-130
    peps = {c: "".join(s) for c, s in ch.items() if 8 <= len(s) <= 11}
    rows.append((fn, {c: len(s) for c, s in ch.items()}, peps))
print("files: %d" % len(rows))
seqs = collections.Counter()
nopep = []
for fn, lens, peps in rows:
    if peps:
        for c, s in peps.items(): seqs[s] += 1
    else:
        nopep.append((fn, lens))
print()
print("=== every 8-11mer sequence found, with count ===")
for s, n in seqs.most_common():
    print("  %-12s %3d" % (s, n))
print()
print("=== NY-ESO-1 check ===")
TARGETS = {"SLLMWITQV": "NY-ESO-1 157-165 (our target)",
           "SLLMWITQC": "NY-ESO-1 157-165 C165 native",
           "SLLMWITQCFL": "NY-ESO-1 157-167"}
hit = False
for t, label in TARGETS.items():
    n = seqs.get(t, 0)
    print("  %-14s %-34s %s" % (t, label, "FOUND x%d" % n if n else "absent"))
    if n: hit = True
# substring sweep, in case of a longer or truncated variant
for s in seqs:
    if "LLMWITQ" in s or "SLLMWIT" in s:
        print("  PARTIAL MATCH: %s" % s); hit = True
print()
print("VERDICT:", "NY-ESO-1 PRESENT" if hit else "no NY-ESO-1 peptide in any scaffold file")
print()
print("=== files with no 8-11mer chain (%d) ===" % len(nopep))
for fn, lens in nopep[:20]: print("  %-44s chains=%s" % (fn, lens))
