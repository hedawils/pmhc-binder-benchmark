import sys, os, csv, re, collections, statistics
from math import comb

A3={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
TARGET = "SLLMWITQV"          # NY-ESO-1 157-165
HOTSPOT = (3, 4, 6, 7)        # 0-based: C4 C5 C7 C8 = M W T Q
PAE_CUT, PLDDT_CUT, RMSD_CUT = 7.0, 92.0, 2.0

# BLOSUM62, enough of it to score any pair
B62_RAW = """A R N D C Q E G H I L K M F P S T W Y V
A 4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0
R -1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3
N -2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3
D -2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3
C 0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2
E -1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2
G 0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3
H -2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3
L -1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1
K -1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2
M -1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1
F -2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2
S 1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2
T 0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1
V 0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4"""
lines = B62_RAW.strip().split("\n"); cols = lines[0].split()
B62 = {}
for row in lines[1:]:
    f = row.split()
    for c, v in zip(cols, f[1:]): B62[(f[0], c)] = int(v)

def blosum(a, b):
    n = min(len(a), len(b))
    return sum(B62.get((a[i], b[i]), 0) for i in range(n)) / n

def fisher(a,b,c,d):
    n=a+b+c+d; r1,r2,c1=a+b,c+d,a+c
    p=lambda x: comb(r1,x)*comb(r2,c1-x)/comb(n,c1)
    obs=p(a); lo=max(0,c1-r2)
    return sum(p(x) for x in range(lo,min(r1,c1)+1) if p(x)<=obs*(1+1e-9))

def wilson(k,n,z=1.96):
    if not n: return (0.0,0.0)
    ph=k/n; d=1+z*z/n; c=(ph+z*z/(2*n))/d
    h=z*((ph*(1-ph)/n+z*z/(4*n*n))**0.5)/d
    return (max(0,c-h)*100, min(1,c+h)*100)

scafdir, man_path, sc_path = sys.argv[1], sys.argv[2], sys.argv[3]

# --- source peptide for EVERY scaffold, including the untyped ones ---------
tail = {}
for fn in os.listdir(scafdir):
    if not fn.endswith(".pdb"): continue
    ch=collections.OrderedDict()
    for l in open(os.path.join(scafdir,fn),errors="replace"):
        if l.startswith("ATOM") and l[12:16].strip()=="CA":
            ch.setdefault(l[21],[]).append(A3.get(l[17:20].strip(),"X"))
    b="".join(ch.get("B",[]))
    tail[fn[:-4]] = b

man = {r["scaffold"]: r for r in csv.DictReader(open(man_path), delimiter="\t")}

# --- AF2 scores -------------------------------------------------------------
rows=[]
with open(sc_path) as fh:
    hdr=fh.readline().split()[1:]
    for line in fh:
        f=line.split()
        if f and f[0]=="SCORE:": rows.append(dict(zip(hdr,f[1:])))
best=collections.defaultdict(list)
for r in rows:
    m=re.match(r"(scaf\d+)", r["description"].replace("liu",""))
    if not m: continue
    ok = (float(r["pae_interaction"])<PAE_CUT and float(r["plddt_binder"])>PLDDT_CUT
          and float(r["binder_aligned_rmsd"])<RMSD_CUT)
    best[m.group(1)].append((float(r["pae_interaction"]), ok))

print("scaffolds with AF2 scores: %d of %d accepted\n" % (len(best), sum(1 for v in man.values() if v["accepted"]=="1")))

# --- peptide assignment, resolving the untyped cohorts ----------------------
recs=[]
for s, scores in best.items():
    m = man.get(s)
    if not m or m["accepted"]!="1": continue
    pep = m["source_peptide"]
    allele = m["source_allele"]
    if pep in ("-", ""):
        t = tail.get(s, "")
        pep = t[-8:] if t else "?"      # untyped: use the chain-B tail as-is
    hit = any(ok for _, ok in scores)
    recs.append(dict(scaf=s, pep=pep, allele=allele,
                     target=m["source_target"], hit=hit,
                     bestpae=min(p for p,_ in scores),
                     nseq=len(scores),
                     focus=float(m["peptide_focus"] or "nan"),
                     blen=int(m["binder_len"]), rmsd=float(m["rmsd"])))

def table(title, keyfn, order=None):
    print("=== %s ===" % title)
    g=collections.defaultdict(list)
    for r in recs: g[keyfn(r)].append(r)
    keys = order or sorted(g, key=lambda k: -len(g[k]))
    print("  %-26s %5s %5s %8s %-16s %s" % ("group","n","hits","rate","95% CI","median best pAE"))
    for k in keys:
        v=g.get(k,[])
        if not v: continue
        h=sum(1 for r in v if r["hit"]); lo,hi=wilson(h,len(v))
        print("  %-26s %5d %5d %7.1f%% [%4.1f-%5.1f]   %.1f"
              % (k,len(v),h,100*h/len(v),lo,hi,statistics.median(r["bestpae"] for r in v)))
    print()
    return g

g_all = table("hit rate by SOURCE ALLELE", lambda r: r["allele"])
table("hit rate by HLA FAMILY (class-I locus)",
      lambda r: ("HLA-A" if r["allele"].startswith("A*") else
                 "HLA-C" if r["allele"].startswith("C*") else "untyped"))
table("hit rate by SOURCE TARGET / peptide cohort",
      lambda r: "%s %s" % (r["target"][:14], r["pep"][:11]))

# --- allele match test ------------------------------------------------------
print("=== does matching our allele (A*02:01) help? ===")
m_ = [r for r in recs if r["allele"]=="A*02:01"]
o_ = [r for r in recs if r["allele"] not in ("A*02:01","unknown")]
u_ = [r for r in recs if r["allele"]=="unknown"]
for lbl, grp in (("A*02:01 (cognate allele)",m_),("other typed allele",o_),("untyped",u_)):
    h=sum(1 for r in grp if r["hit"]); lo,hi=wilson(h,len(grp)) if grp else (0,0)
    print("  %-26s %3d scaffolds  %3d hits  %5.1f%%  [%.1f-%.1f]"
          % (lbl,len(grp),h,100*h/len(grp) if grp else 0,lo,hi))
if m_ and o_:
    a=sum(1 for r in m_ if r["hit"]); c=sum(1 for r in o_ if r["hit"])
    print("  Fisher A*02:01 vs other typed: p = %.3f" % fisher(a,len(m_)-a,c,len(o_)-c))
print()

# --- peptide similarity -----------------------------------------------------
print("=== does cognate-peptide similarity to %s predict a hit? ===" % TARGET)
for name, fn in (("whole peptide, BLOSUM62/residue", lambda p: blosum(p, TARGET)),
                 ("hotspot C4/C5/C7/C8 only",
                  lambda p: blosum("".join(p[i] for i in HOTSPOT if i < len(p)),
                                   "".join(TARGET[i] for i in HOTSPOT))),
                 ("identities, whole peptide",
                  lambda p: sum(1 for i in range(min(len(p),9)) if p[i]==TARGET[i]))):
    hits=[fn(r["pep"]) for r in recs if r["hit"] and "?" not in r["pep"]]
    miss=[fn(r["pep"]) for r in recs if not r["hit"] and "?" not in r["pep"]]
    if not hits or not miss: continue
    print("  %-34s hits %5.2f (n=%d)   misses %5.2f (n=%d)   diff %+.2f"
          % (name, statistics.mean(hits), len(hits),
             statistics.mean(miss), len(miss),
             statistics.mean(hits)-statistics.mean(miss)))
print()

# --- covariates, so similarity is not credited for something else -----------
print("=== covariates: what else differs between hits and misses? ===")
for name, key in (("static peptide_focus","focus"),("binder length","blen"),
                  ("superposition RMSD","rmsd"),("sequences scored","nseq")):
    h=[r[key] for r in recs if r["hit"] and r[key]==r[key]]
    m2=[r[key] for r in recs if not r["hit"] and r[key]==r[key]]
    if h and m2:
        print("  %-22s hits %7.3f   misses %7.3f   diff %+.3f"
              % (name, statistics.mean(h), statistics.mean(m2),
                 statistics.mean(h)-statistics.mean(m2)))
