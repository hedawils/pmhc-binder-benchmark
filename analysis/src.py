import csv, math, re
from collections import defaultdict
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _repo import B, ROOT, sibling
FOCUS=0.202

def norm(bid):
    """liuscaf100_0 / kf9scaf100_0 -> scaf100"""
    m=re.search(r"(scaf\d+)", bid)
    return m.group(1) if m else bid

def surv_scaffolds(key):
    """scaffolds with >=1 four-gate survivor, and all scaffolds attempted"""
    rows=list(csv.DictReader(open(B+key+"_metrics.tsv"),delimiter="\t"))
    allb=set(); win=set()
    for r in rows:
        allb.add(norm(r["backbone_id"]))
        try:
            if (float(r["scrmsd"])<2.0 and float(r["plddt_binder"])>92
                and float(r["pae_interaction"])<7.0 and float(r["peptide_focus"])>=FOCUS):
                win.add(norm(r["backbone_id"]))
        except (TypeError,ValueError): pass
    return allb, win

def manifest(key):
    m={}
    for r in csv.DictReader(open(B+"manifest_"+key+".tsv"),delimiter="\t"):
        if r["accepted"]=="1": m[r["scaffold"]]=(r["source_allele"], r["source_target"])
    return m

def fisher(a,b,c,d):
    def lc(n): return math.lgamma(n+1)
    def p(a,b,c,d):
        return math.exp(lc(a+b)+lc(c+d)+lc(a+c)+lc(b+d)-lc(a)-lc(b)-lc(c)-lc(d)-lc(a+b+c+d))
    p0=p(a,b,c,d); tot=0.0
    for i in range(0,a+b+1):
        j=a+b-i; k=a+c-i; l=d-a+i
        if k<0 or l<0: continue
        pi=p(i,j,k,l)
        if pi<=p0*(1+1e-9): tot+=pi
    return min(tot,1.0)

def locus(al):
    al=(al or "").strip()
    if al.startswith("A"): return "A"
    if al.startswith("B"): return "B"
    if al.startswith("C"): return "C"
    return "unrecorded"

for tgt,key,tgt_locus in (("NY-ESO-1  (target locus A)","liu","A"),
                          ("KF9       (target locus C)","kf9","C")):
    allb,win = surv_scaffolds(key); man = manifest(key)
    grp=defaultdict(lambda:[0,0])
    for s in allb:
        if s not in man: continue
        L=locus(man[s][0])
        grp[L][0]+=1
        if s in win: grp[L][1]+=1
    print(f"\n=== {tgt} ===   scaffolds attempted {sum(g[0] for g in grp.values())}")
    for L in ("A","B","C","unrecorded"):
        if L not in grp: continue
        n,w=grp[L]
        tag=" <= MATCHES TARGET LOCUS" if L==tgt_locus else ""
        print(f"  source locus {L:11} {w:3}/{n:3} = {100*w/n:5.1f}%{tag}")
    if "A" in grp and "C" in grp:
        (na,wa),(nc,wc)=grp["A"],grp["C"]
        p=fisher(wc,nc-wc,wa,na-wa)
        print(f"  C vs A sourced: {100*wc/nc:.1f}% vs {100*wa/na:.1f}%   Fisher p = {p:.5g}")
