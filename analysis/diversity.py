import sys, os, csv, re, collections, statistics
A3={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
PAE,PL,RM = 7.0, 92.0, 2.0
scafdir, man_path, sc_path = sys.argv[1], sys.argv[2], sys.argv[3]

binder, chainb = {}, {}
for fn in os.listdir(scafdir):
    if not fn.endswith(".pdb"): continue
    ch=collections.OrderedDict()
    for l in open(os.path.join(scafdir,fn),errors="replace"):
        if l.startswith("ATOM") and l[12:16].strip()=="CA":
            ch.setdefault(l[21],[]).append(A3.get(l[17:20].strip(),"X"))
    binder[fn[:-4]]="".join(ch.get("A",[])); chainb[fn[:-4]]=len(ch.get("B",[]))

man={r["scaffold"]:r for r in csv.DictReader(open(man_path),delimiter="\t")}
rows=[]
with open(sc_path) as fh:
    hdr=fh.readline().split()[1:]
    for line in fh:
        f=line.split()
        if f and f[0]=="SCORE:": rows.append(dict(zip(hdr,f[1:])))
best=collections.defaultdict(list)
for r in rows:
    m=re.match(r"(scaf\d+)", r["description"].replace("liu",""))
    if m: best[m.group(1)].append(float(r["pae_interaction"])<PAE
          and float(r["plddt_binder"])>PL and float(r["binder_aligned_rmsd"])<RM)

recs=[]
for s,oks in best.items():
    m=man.get(s)
    if not m or m["accepted"]!="1": continue
    recs.append(dict(s=s, hit=any(oks), tgt=m["source_target"],
                     allele=m["source_allele"], bl=int(m["binder_len"]),
                     cb=chainb.get(s,0), seq=binder.get(s,"")))

print("=== CONFOUND 1: chain-B construction (platform truncation) ===")
g=collections.defaultdict(list)
for r in recs: g["full (>=189)" if r["cb"]>=189 else "truncated (<=184)"].append(r)
for k,v in sorted(g.items()):
    h=sum(1 for r in v if r["hit"])
    print("  %-20s n=%-4d hits=%-3d %5.1f%%" % (k,len(v),h,100*h/len(v)))
print()
print("  same split, within each source cohort:")
byt=collections.defaultdict(lambda: collections.defaultdict(list))
for r in recs: byt[r["tgt"]]["full" if r["cb"]>=189 else "trunc"].append(r)
for t,sub in sorted(byt.items(), key=lambda kv:-sum(len(x) for x in kv[1].values())):
    parts=[]
    for k in ("full","trunc"):
        v=sub.get(k,[])
        if v: parts.append("%s n=%d hits=%d (%.0f%%)" % (k,len(v),sum(1 for r in v if r["hit"]),100*sum(1 for r in v if r["hit"])/len(v)))
    print("    %-12s %s" % (t, "   ".join(parts)))
print()

def ident(a,b):
    n=min(len(a),len(b))
    return sum(1 for i in range(n) if a[i]==b[i])/n if n else 0.0

print("=== CONFOUND 2: are the winning scaffolds redundant? ===")
print("  pairwise binder-sequence identity within each cohort (n>=5)")
print("  %-14s %4s %-22s %s" % ("cohort","n","median / max identity","distinct at <90%"))
for t,v in sorted(byt.items(), key=lambda kv:-sum(len(x) for x in kv[1].values())):
    allv=[r for sub in v.values() for r in sub]
    if len(allv)<5: continue
    ids=[ident(allv[i]["seq"],allv[j]["seq"]) for i in range(len(allv)) for j in range(i+1,len(allv))]
    # greedy clustering at 90% identity
    reps=[]
    for r in allv:
        if not any(ident(r["seq"],q["seq"])>=0.90 for q in reps): reps.append(r)
    print("  %-14s %4d %6.2f / %.2f            %d" %
          (t,len(allv),statistics.median(ids),max(ids),len(reps)))
print()
print("=== the PHOX2B cohort, unpacked ===")
ph=[r for r in recs if r["tgt"]=="PHOX2B"]
hits=[r for r in ph if r["hit"]]
reps=[]
for r in hits:
    if not any(ident(r["seq"],q["seq"])>=0.90 for q in reps): reps.append(r)
print("  %d scaffolds, %d hits" % (len(ph), len(hits)))
print("  hit binder lengths: %s" % sorted(r["bl"] for r in hits))
print("  distinct hit families at <90%% identity: %d" % len(reps))
ids=[ident(hits[i]["seq"],hits[j]["seq"]) for i in range(len(hits)) for j in range(i+1,len(hits))]
if ids: print("  pairwise identity among hits: median %.2f  max %.2f" % (statistics.median(ids),max(ids)))
