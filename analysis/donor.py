import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _repo import B, ROOT, sibling
import csv, re
from collections import defaultdict
exec(open(sibling("src.py")).read().split("for tgt,key")[0])

res=defaultdict(dict); meta={}
for key in ("liu","kf9"):
    allb,win=surv_scaffolds(key)
    for r in csv.DictReader(open(B+"manifest_"+key+".tsv"),delimiter="\t"):
        if r["accepted"]!="1" or r["scaffold"] not in allb: continue
        s=r["source_target"]; meta[s]=r["source_allele"] or "—"
        d=res[s].setdefault(key,[0,0]); d[0]+=1
        if r["scaffold"] in win: d[1]+=1

rows=[]
for s,d in res.items():
    n1,w1=d.get("liu",[0,0]); n2,w2=d.get("kf9",[0,0])
    rows.append((s,meta[s],n1,w1,n2,w2))
rows.sort(key=lambda r:-(r[3]/r[2] if r[2] else 0))
print(f"{'donor target':16} {'source allele':14} {'NY-ESO-1':>16} {'KF9':>16}")
for s,a,n1,w1,n2,w2 in rows:
    f1=f"{w1}/{n1} = {100*w1/n1:4.1f}%" if n1 else "—"
    f2=f"{w2}/{n2} = {100*w2/n2:4.1f}%" if n2 else "—"
    print(f"{s:16} {a:14} {f1:>16} {f2:>16}")
t1=(sum(r[3] for r in rows), sum(r[2] for r in rows))
t2=(sum(r[5] for r in rows), sum(r[4] for r in rows))
print(f"{'ALL BORROWED':16} {'':14} {t1[0]}/{t1[1]} = {100*t1[0]/t1[1]:4.1f}%".ljust(64))
print(f"{'KF9 total':16} {'':14} {t2[0]}/{t2[1]} = {100*t2[0]/t2[1]:4.1f}%")
