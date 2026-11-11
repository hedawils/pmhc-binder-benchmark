"""Can we build an HLA-C target straight out of the PHOX2B scaffold files,
and derive its hotspots from the native binder's own contacts?"""
import sys, os, csv, collections
import numpy as np
A3={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
BB={"N","CA","C","O","OXT"}; CONTACT=4.5
man={r["scaffold"]:r for r in csv.DictReader(open(sys.argv[1]),delimiter="\t")}
scafdir=sys.argv[2]
ph=[s for s,r in man.items() if r["source_target"]=="PHOX2B" and r["accepted"]=="1"]
print("PHOX2B scaffolds available as C*07:02 target donors: %d\n" % len(ph))

lens=collections.Counter(); hot=collections.Counter(); tot=collections.Counter(); seqs=collections.Counter()
for s in ph:
    p=os.path.join(scafdir,s+".pdb")
    if not os.path.exists(p): continue
    a=[]; b=collections.OrderedDict()
    for l in open(p):
        if not l.startswith("ATOM"): continue
        if l[76:78].strip()=="H": continue
        xyz=(float(l[30:38]),float(l[38:46]),float(l[46:54]))
        if l[21]=="A": a.append(xyz)
        elif l[21]=="B": b.setdefault((l[22:26],l[26]),[]).append((l[12:16].strip(),l[17:20].strip(),xyz))
    lens[len(b)]+=1
    keys=list(b)[-9:]
    seqs["".join(A3.get(b[k][0][1],"X") for k in keys)]+=1
    P=[]; idx=[]
    for i,k in enumerate(keys):
        for nm,_,xyz in b[k]:
            P.append(xyz); idx.append((i+1,nm))
    A=np.array(a); P=np.array(P)
    d=np.linalg.norm(A[:,None]-P[None],axis=-1); hits=(d<CONTACT).sum(axis=0)
    for (pos,nm),n in zip(idx,hits):
        tot[pos]+=n
        if nm not in BB: hot[pos]+=n

print("chain B lengths:", dict(lens))
print("peptide read from chain B tail:", dict(seqs))
pep=max(seqs,key=seqs.get)
print()
print("=== native-binder contacts on QYNPIRTTF (the hotspot derivation) ===")
print("  %-6s %10s %12s" % ("pos","side chain","all contacts"))
for p in sorted(tot):
    print("  %s%-5d %10d %12d" % (pep[p-1],p,hot[p],tot[p]))
top=sorted(hot,key=lambda k:-hot[k])[:4]
print()
print("  -> hotspots by side-chain contact: %s" % ", ".join("C%d (%s%d)"%(p,pep[p-1],p) for p in sorted(top)))
print("  -> ours for NY-ESO-1 were C4 C5 C7 C8 (M4 W5 T7 Q8)")
