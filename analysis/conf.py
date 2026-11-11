import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _repo import B, ROOT, sibling
import csv, math, re, statistics as st
from collections import defaultdict
exec(open(sibling("src.py")).read().split("for tgt,key")[0])   # norm, surv_scaffolds, manifest, fisher, locus

FEATS=["binder_len","pep_contacts","mhc_contacts","peptide_focus","buried_peptide","rmsd","pep_min"]

def table(key):
    allb,win=surv_scaffolds(key); m={}
    for r in csv.DictReader(open(B+"manifest_"+key+".tsv"),delimiter="\t"):
        if r["accepted"]!="1" or r["scaffold"] not in allb: continue
        try: d={f: float(r[f]) for f in FEATS}
        except ValueError: continue
        d["locus"]=locus(r["source_allele"]); d["src"]=r["source_target"]
        d["win"]= r["scaffold"] in win
        m[r["scaffold"]]=d
    return m

def mw(a,b):
    r,ties=ranks(a+b); n1,n2=len(a),len(b)
    if n1==0 or n2==0: return None
    U=sum(r[:n1])-n1*(n1+1)/2.0; mu=n1*n2/2.0; N=n1+n2
    tie=sum(t**3-t for t in ties)
    var=n1*n2/12.0*((N+1)-tie/float(N*(N-1)))
    if var<=0: return None
    z=(abs(U-mu)-0.5)/math.sqrt(var)
    return math.erfc(z/math.sqrt(2.0))
def ranks(v):
    idx=sorted(range(len(v)),key=lambda i:v[i]); r=[0.0]*len(v); i=0; ties=[]
    while i<len(idx):
        j=i
        while j+1<len(idx) and v[idx[j+1]]==v[idx[i]]: j+=1
        a=(i+j)/2.0+1.0
        for k in range(i,j+1): r[idx[k]]=a
        ties.append(j-i+1); i=j+1
    return r,ties

for tgt,key in (("NY-ESO-1","liu"),("KF9","kf9")):
    m=table(key)
    W=[d for d in m.values() if d["win"]]; L=[d for d in m.values() if not d["win"]]
    print(f"\n######## {tgt}  winners n={len(W)}  losers n={len(L)} ########")
    print(f"  {'transplant-stage feature':24} {'winners':>10} {'losers':>10} {'p':>10}")
    for f in FEATS:
        a=[d[f] for d in W]; b=[d[f] for d in L]
        p=mw(a,b)
        print(f"  {f:24} {st.median(a):10.2f} {st.median(b):10.2f} {p:10.4g}")
    # does the feature explain the A/C gap?
    print("  -- same features, C-sourced vs A-sourced scaffolds (regardless of outcome) --")
    C=[d for d in m.values() if d["locus"]=="C"]; A=[d for d in m.values() if d["locus"]=="A"]
    for f in FEATS:
        p=mw([d[f] for d in C],[d[f] for d in A])
        print(f"  {f:24} {st.median([d[f] for d in C]):10.2f} {st.median([d[f] for d in A]):10.2f} {p:10.4g}")
    # which source targets are the C cohort?
    print("  -- source targets by locus --")
    g=defaultdict(lambda:[0,0])
    for d in m.values():
        g[(d["locus"],d["src"])][0]+=1
        if d["win"]: g[(d["locus"],d["src"])][1]+=1
    for (L_,s),(n,w) in sorted(g.items()):
        print(f"     {L_:11} {s:14} {w:3}/{n:3} = {100*w/n:5.1f}%")

print("\n\n######## within-cohort check: drop PHOX2B entirely ########")
for tgt,key in (("NY-ESO-1","liu"),("KF9","kf9")):
    m={k:v for k,v in table(key).items() if v["src"]!="PHOX2B"}
    W=[d for d in m.values() if d["win"]]; L=[d for d in m.values() if not d["win"]]
    print(f"\n{tgt}: non-PHOX2B scaffolds, winners n={len(W)}, losers n={len(L)}")
    for f in ("pep_contacts","mhc_contacts","peptide_focus","rmsd"):
        a=[d[f] for d in W]; b=[d[f] for d in L]
        if not a or not b: continue
        print(f"  {f:16} win {st.median(a):8.2f}  lose {st.median(b):8.2f}  p={mw(a,b):.4g}")
    # spread across donor source targets
    g=defaultdict(lambda:[0,0])
    for d in table(key).values():
        g[d["src"]][0]+=1
        if d["win"]: g[d["src"]][1]+=1
    rates=[(100*w/n,s,w,n) for s,(n,w) in g.items() if n>=6]
    rates.sort(reverse=True)
    print(f"  donor spread (n>=6): "+"  ".join(f"{s} {r:.0f}% ({w}/{n})" for r,s,w,n in rates))
