import sys, csv, collections, statistics
from math import comb
def fisher(a,b,c,d):
    n=a+b+c+d; r1,r2,c1=a+b,c+d,a+c
    p=lambda x: comb(r1,x)*comb(r2,c1-x)/comb(n,c1)
    obs=p(a); lo=max(0,c1-r2)
    return sum(p(x) for x in range(lo,min(r1,c1)+1) if p(x)<=obs*(1+1e-9))
man={r["scaffold"]:r for r in csv.DictReader(open(sys.argv[1]),delimiter="\t")}
rows=list(csv.DictReader(open(sys.argv[2]),delimiter="\t"))
T={"scrmsd":("<",2.0),"plddt_binder":(">",92.0),
   "pae_interaction":("<",7.0),"peptide_focus":(">=",0.202)}
def ok(r):
    try:
        for k,(op,c) in T.items():
            v=float(r[k])
            if op=="<" and not v<c: return False
            if op==">" and not v>c: return False
            if op==">=" and not v>=c: return False
        return True
    except (ValueError,KeyError): return False
bb=collections.defaultdict(list)
for r in rows:
    bb[r["backbone_id"]].append((int(r["seq_id"]), ok(r)))
def scaf(b): return b.replace("liu","",1).rsplit("_",1)[0]
recs=[]
for b,v in bb.items():
    s=scaf(b)
    if s not in man or man[s]["accepted"]!="1": continue
    recs.append(dict(s=s, hit6=any(o for _,o in v), hit4=any(o for i,o in v if i<4),
                     allele=man[s]["source_allele"], tgt=man[s]["source_target"]))
print("FINAL four-filter table: %d designs, %d backbones\n" % (len(rows), len(recs)))
def fam(a): return "HLA-A (ours)" if a.startswith("A*") else "HLA-C" if a.startswith("C*") else "not recorded"
for name,key in (("depth 6","hit6"),("depth 4 (matched to arm 5)","hit4")):
    g=collections.defaultdict(list)
    for r in recs: g[fam(r["allele"])].append(r[key])
    tot=sum(len(v) for v in g.values()); th=sum(sum(v) for v in g.values())
    print("=== by locus family, %s ===" % name)
    for k in ("HLA-C","not recorded","HLA-A (ours)"):
        v=g[k]; h=sum(v)
        print("  %-16s %3d shapes  %3d worked  %5.1f%%   1 in %.1f" % (k,len(v),h,100*h/len(v),len(v)/h if h else float("inf")))
    print("  %-16s %3d shapes  %3d worked  %5.1f%%\n" % ("TOTAL",tot,th,100*th/tot))
a=[r for r in recs if r["allele"]=="A*02:01"]; o=[r for r in recs if r["allele"] not in ("A*02:01","unknown")]
ah=sum(r["hit6"] for r in a); oh=sum(r["hit6"] for r in o)
print("A*02:01 %d/%d (%.1f%%) vs other typed %d/%d (%.1f%%)  Fisher p=%.4f"
      % (ah,len(a),100*ah/len(a),oh,len(o),100*oh/len(o),fisher(ah,len(a)-ah,oh,len(o)-oh)))
pap=[r for r in recs if r["tgt"]=="PAP"]
print("PAP (A*02:01): %d shapes, %d worked" % (len(pap), sum(r["hit6"] for r in pap)))
