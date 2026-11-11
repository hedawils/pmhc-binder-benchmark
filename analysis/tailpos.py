import sys, os, collections
A3={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
d=sys.argv[1]; groups=collections.defaultdict(list)
for fn in sorted(os.listdir(d)):
    if not fn.endswith(".pdb"): continue
    ch=collections.OrderedDict()
    for l in open(os.path.join(d,fn),errors="replace"):
        if l.startswith("ATOM") and l[12:16].strip()=="CA":
            ch.setdefault(l[21],[]).append(A3.get(l[17:20].strip(),"X"))
    b="".join(ch.get("B",[]))
    groups[b[-11:]].append((fn,len(b),b))
print("=== each tail group: chain-B length, and where the platform ends ===")
for tail,items in sorted(groups.items(),key=lambda kv:-len(kv[1])):
    lens=collections.Counter(l for _,l,_ in items)
    ex=items[0][2]
    print("  %-12s n=%-3d chainB len=%s" % (tail,len(items),dict(lens)))
    print("      residues 176-end: %s" % ex[175:])
