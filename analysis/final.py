import sys, csv, re, collections
from math import comb
sys.path.insert(0,"/dev/null")
def fisher(a,b,c,d):
    n=a+b+c+d; r1,r2,c1=a+b,c+d,a+c
    p=lambda x: comb(r1,x)*comb(r2,c1-x)/comb(n,c1)
    obs=p(a); lo=max(0,c1-r2)
    return sum(p(x) for x in range(lo,min(r1,c1)+1) if p(x)<=obs*(1+1e-9))
def wilson(k,n,z=1.96):
    if not n: return (0,0)
    ph=k/n; d=1+z*z/n; c=(ph+z*z/(2*n))/d
    h=z*((ph*(1-ph)/n+z*z/(4*n*n))**0.5)/d
    return (max(0,c-h)*100,min(1,c+h)*100)
A3={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
import os
scafdir,man_path,sc_path=sys.argv[1],sys.argv[2],sys.argv[3]
cb={}
for fn in os.listdir(scafdir):
    if fn.endswith(".pdb"):
        n=0
        for l in open(os.path.join(scafdir,fn),errors="replace"):
            if l.startswith("ATOM") and l[12:16].strip()=="CA" and l[21]=="B": n+=1
        cb[fn[:-4]]=n
man={r["scaffold"]:r for r in csv.DictReader(open(man_path),delimiter="\t")}
rows=[]
with open(sc_path) as fh:
    hdr=fh.readline().split()[1:]
    for line in fh:
        f=line.split()
        if f and f[0]=="SCORE:": rows.append(dict(zip(hdr,f[1:])))
best=collections.defaultdict(list)
for r in rows:
    m=re.match(r"(scaf\d+)",r["description"].replace("liu",""))
    if m: best[m.group(1)].append(float(r["pae_interaction"])<7.0
         and float(r["plddt_binder"])>92.0 and float(r["binder_aligned_rmsd"])<2.0)
recs=[dict(s=s,hit=any(v),tgt=man[s]["source_target"],cb=cb.get(s,0))
      for s,v in best.items() if s in man and man[s]["accepted"]=="1"]

full=[r for r in recs if r["cb"]>=189]
ph=[r for r in full if r["tgt"]=="PHOX2B"]
ot=[r for r in full if r["tgt"]!="PHOX2B"]
a=sum(1 for r in ph if r["hit"]); c=sum(1 for r in ot if r["hit"])
print("=== within FULL-PLATFORM scaffolds only, is PHOX2B still special? ===")
for lbl,grp,k in (("PHOX2B",ph,a),("every other full-platform source",ot,c)):
    lo,hi=wilson(k,len(grp))
    print("  %-34s %2d/%2d  %5.1f%%  [%.1f-%.1f]" % (lbl,k,len(grp),100*k/len(grp),lo,hi))
print("  Fisher p = %.4f" % fisher(a,len(ph)-a,c,len(ot)-c))
print()
print("=== overall ===")
h=sum(1 for r in recs if r["hit"])
print("  %d of %d scaffolds yield at least one 3-filter hit (%.1f%%)" % (h,len(recs),100*h/len(recs)))
print("  PHOX2B contributes %d of those %d hits (%.0f%%) from %d%% of the scaffolds"
      % (a,h,100*a/h,round(100*len(ph)/len(recs))))
