"""Is PHOX2B's higher side-chain fraction just a Trp5 packing effect?

The side-chain fraction is weighted by side-chain SIZE (Trp has 10 heavy atoms,
Ala 1, Gly 0), so a cohort that merely sits over W5 more would score higher
without reading more chemistry. Recompute excluding W5, and report each
cohort's W5 share, to separate the two.
"""
import sys, csv, collections, random, statistics
import numpy as np
BB = {"N","CA","C","O","OXT"}; CONTACT = 4.5; PEP_LEN = 9
T = {"scrmsd":("<",2.0),"plddt_binder":(">",92.0),
     "pae_interaction":("<",7.0),"peptide_focus":(">=",0.202)}
def passes(r):
    try:
        for k,(op,c) in T.items():
            v=float(r[k])
            if op=="<" and not v<c: return False
            if op==">" and not v>c: return False
            if op==">=" and not v>=c: return False
        return True
    except (ValueError,KeyError): return False
def load(path):
    a,b=[],collections.OrderedDict()
    for l in open(path):
        if not l.startswith("ATOM"): continue
        if l[76:78].strip()=="H" or l[12:16].strip().startswith("H"): continue
        xyz=(float(l[30:38]),float(l[38:46]),float(l[46:54]))
        if l[21]=="A": a.append(xyz)
        elif l[21]=="B": b.setdefault((l[22:26],l[26]),[]).append((l[12:16].strip(),xyz))
    keys=list(b)[-PEP_LEN:]
    return np.array(a), [(i+1,nm,xyz) for i,k in enumerate(keys) for nm,xyz in b[k]]
man={r["scaffold"]:r for r in csv.DictReader(open(sys.argv[1]),delimiter="\t")}
rows=[r for r in csv.DictReader(open(sys.argv[2]),delimiter="\t") if passes(r)]
recs=[]
for r in rows:
    try: A,pep=load("%s/%s.pdb"%(sys.argv[3],r["description"]))
    except OSError: continue
    P=np.array([x for _,_,x in pep]); d=np.linalg.norm(A[:,None]-P[None],axis=-1)
    hit=(d<CONTACT).sum(axis=0)
    c=collections.Counter()
    for (pos,nm,_),n in zip(pep,hit):
        c["bb" if nm in BB else "sc"]+=n
        if pos==5: c["w5_sc" if nm not in BB else "w5_bb"]+=n
    if c["bb"]+c["sc"]==0: continue
    scaf=r["backbone_id"].replace("liu","",1).rsplit("_",1)[0]
    bb_x=c["bb"]-c["w5_bb"]; sc_x=c["sc"]-c["w5_sc"]
    recs.append(dict(tgt=man.get(scaf,{}).get("source_target","?"),
                     f=c["sc"]/(c["bb"]+c["sc"]),
                     f_noW5=(sc_x/(bb_x+sc_x)) if (bb_x+sc_x) else float("nan"),
                     w5share=(c["w5_sc"]+c["w5_bb"])/(c["bb"]+c["sc"]),
                     tot=c["bb"]+c["sc"]))
def perm(a,b,n=20000):
    obs=statistics.mean(a)-statistics.mean(b); pool=a+b; random.seed(0); h=0
    for _ in range(n):
        random.shuffle(pool)
        if abs(statistics.mean(pool[:len(a)])-statistics.mean(pool[len(a):]))>=abs(obs): h+=1
    return obs,h/n
print("  %-14s %4s %9s %12s %11s" % ("cohort","n","s.c. frac","s.c. ex-W5","W5 share"))
g=collections.defaultdict(list)
for x in recs: g[x["tgt"]].append(x)
for k,v in sorted(g.items(),key=lambda kv:-len(kv[1])):
    ok=[x for x in v if x["f_noW5"]==x["f_noW5"]]
    print("  %-14s %4d %9.3f %12.3f %10.3f" % (k,len(v),
        statistics.mean(x["f"] for x in v),
        statistics.mean(x["f_noW5"] for x in ok) if ok else float("nan"),
        statistics.mean(x["w5share"] for x in v)))
ph=[x for x in recs if x["tgt"]=="PHOX2B"]; ot=[x for x in recs if x["tgt"]!="PHOX2B"]
print()
for lbl,key in (("side-chain fraction","f"),("side-chain fraction EXCLUDING W5","f_noW5"),
                ("W5 share of all contacts","w5share")):
    a=[x[key] for x in ph if x[key]==x[key]]; b=[x[key] for x in ot if x[key]==x[key]]
    o,p=perm(a,b)
    print("  %-34s PHOX2B %.3f vs rest %.3f   diff %+.3f   p=%.4f" % (lbl,statistics.mean(a),statistics.mean(b),o,p))
