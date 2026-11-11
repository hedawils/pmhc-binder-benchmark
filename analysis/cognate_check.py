"""Does any candidate target peptide appear ANYWHERE in any scaffold's chain B?

The tail-based typing left 63 scaffolds untyped: their chain B is 183-184
residues and the tail does not sit after the same platform boundary as the
typed ones, so I could not say where platform ends and peptide begins. A
cognate target could therefore hide in that group. Searching the whole chain B
as a substring is boundary-independent and settles it.
"""
import sys, os, collections
A3={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
CAND = [("KF9  SARS-CoV-2 N266-274", "KAYNVTQAF"),
        ("QYDDAVYKL  C*04:01",       "QYDDAVYKL"),
        ("GILGFVFTL  C*08:01",       "GILGFVFTL"),
        ("GADGVGKSA  KRAS G12D",     "GADGVGKSA"),
        ("GADGVGKSAL KRAS 10-mer",   "GADGVGKSAL"),
        ("SLLMWITQV  NY-ESO-1",      "SLLMWITQV"),
        ("QYNPIRTTF  PHOX2B",        "QYNPIRTTF")]
d = sys.argv[1]
chains = {}
for fn in sorted(os.listdir(d)):
    if not fn.endswith(".pdb"): continue
    s = []
    for l in open(os.path.join(d, fn), errors="replace"):
        if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[21] == "B":
            s.append(A3.get(l[17:20].strip(), "X"))
    chains[fn[:-4]] = "".join(s)
print("scaffolds searched: %d   chain-B lengths: %s\n"
      % (len(chains), dict(collections.Counter(len(v) for v in chains.values()))))
print("  %-28s %-11s %s" % ("candidate", "peptide", "scaffolds containing it"))
for label, pep in CAND:
    hits = [k for k, v in chains.items() if pep in v]
    print("  %-28s %-11s %s" % (label, pep, ("%d  %s" % (len(hits), " ".join(sorted(hits)[:6]))) if hits else "0  — not cognate"))
print()
# also: the last 12 of every distinct untyped-length chain, for the record
untyped = sorted({v[-12:] for v in chains.values() if len(v) < 186})
print("  tails of the short (untyped) chain Bs, for reference:")
for t in untyped: print("    ...%s" % t)
