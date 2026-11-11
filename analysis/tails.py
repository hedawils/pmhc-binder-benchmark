import sys, os, collections
A3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
      "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
      "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
KNOWN = {"RMFPNAPYL":"WT1","KLTPLCVTL":"HIV","GLMWLSYFV":"SARS-CoV1",
         "LLWNGPIAV":"YFV","TLMSAMTNL":"PAP","YLEPGPVTA":"gp100",
         "ELAGIGILTV":"MART-1","ALYVDSLFFL":"PRAME","EVDPIGHLY":"MAGE-A3",
         "ESDPIVAQY":"Titin","QYNPIRTTF":"PHOX2B","TTAPFLSGK":"CTNNb1-S5F",
         "ALHGGWTTK":"seed scaffold"}
d = sys.argv[1]
tails = collections.Counter()
per = {}
for fn in sorted(os.listdir(d)):
    if not fn.endswith(".pdb"): continue
    ch = collections.OrderedDict()
    with open(os.path.join(d, fn), errors="replace") as fh:
        for l in fh:
            if l.startswith("ATOM") and l[12:16].strip() == "CA":
                ch.setdefault(l[21], []).append(A3.get(l[17:20].strip(), "X"))
    b = "".join(ch.get("B", []))
    # the longest tail that matches a known target, else the last 11
    label, pep = None, b[-11:]
    for n in (10, 9, 8):
        if b[-n:] in KNOWN:
            label, pep = KNOWN[b[-n:]], b[-n:]; break
    tails[(pep, label)] += 1
    per[fn] = (pep, label, len(b))
print("scaffolds read: %d\n" % len(per))
print("=== C-terminal tail of chain B, grouped ===")
print("  %-13s %-16s %s" % ("tail", "typed as", "n"))
for (pep, label), n in tails.most_common():
    print("  %-13s %-16s %d" % (pep, label or "UNTYPED", n))
print("\n=== NY-ESO-1 sweep over every tail (any length) ===")
bad = []
for fn, (pep, label, bl) in per.items():
    full = pep
    for motif in ("SLLMWITQ", "LLMWITQ", "LMWITQ", "MWITQ"):
        if motif in full:
            bad.append((fn, full, motif)); break
if bad:
    for fn, full, m in bad: print("  MATCH %-14s tail=%s  motif=%s" % (fn, full, m))
else:
    print("  no scaffold's chain-B tail contains any NY-ESO-1 157-165 motif")
print("\nVERDICT:", "CONTAMINATED" if bad else
      "none of the %d scaffolds derives from an NY-ESO-1 157-165 binder" % len(per))
