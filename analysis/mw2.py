import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _repo import B, ROOT, sibling
import csv, math, statistics as st
exec(open(sibling("mw.py")).read().split("METRICS = [")[0])   # reuse load/ranks/mannwhitney/med/stars

def frac(rows, key, thr=0.0):
    v=[r[key] for r in rows if r[key] is not None]
    return sum(1 for x in v if x>thr)/len(v)*100, len(v)

def twoprop(a,b):
    """chi-square 1df on two proportions, returns p"""
    (p1,n1),(p2,n2)=a,b
    x1,x2=round(p1*n1/100),round(p2*n2/100)
    n=n1+n2; x=x1+x2
    if x==0 or x==n: return None
    e=[n1*x/n, n1*(n-x)/n, n2*x/n, n2*(n-x)/n]
    o=[x1,n1-x1,x2,n2-x2]
    chi=sum((o[i]-e[i])**2/e[i] for i in range(4))
    return math.erfc(math.sqrt(chi/2.0))

for tname,pkey,fkey in (("NY-ESO-1 / A*02:01","liu","free"),
                        ("SARS-CoV-2 KF9 / C*12:02","kf9","kf9free")):
    P,F=load(pkey),load(fkey)
    print(f"\n=== {tname}  (all designs: partial n={len(P)}, free n={len(F)}) ===")
    for key,label in (("buried_peptide","any peptide burial"),):
        a,b=frac(P,key),frac(F,key)
        print(f"  {label:24} partial {a[0]:5.1f}%   free {b[0]:5.1f}%   p={twoprop(a,b):.2e}")
    # conditional on having peptide contact
    Pc=[r for r in P if (r['buried_peptide'] or 0)>0]
    Fc=[r for r in F if (r['buried_peptide'] or 0)>0]
    print(f"  -- among designs that DO touch the peptide (partial n={len(Pc)}, free n={len(Fc)}) --")
    for key,label,dirn in (("buried_peptide","buried peptide A^2","hi"),
                           ("peptide_focus","peptide focus","hi"),
                           ("pae_interaction","interface pAE","lo")):
        mp,vp=med(Pc,key); mf,vf=med(Fc,key)
        _,p=mannwhitney(vp,vf)
        print(f"  {label:24} partial {mp:8.3f}   free {mf:8.3f}   d={mp-mf:+8.3f}  p={p:.4g} {stars(p)}")
