import sys, csv, re, collections, statistics
man={r["scaffold"]:r for r in csv.DictReader(open(sys.argv[1]),delimiter="\t")}
rows=[]
with open(sys.argv[2]) as fh:
    hdr=fh.readline().split()[1:]
    for line in fh:
        f=line.split()
        if f and f[0]=="SCORE:": rows.append(dict(zip(hdr,f[1:])))
d=[]
for r in rows:
    m=re.match(r"(scaf\d+)", r["description"].replace("liu","",1))
    if not m or m.group(1) not in man: continue
    d.append((int(man[m.group(1)]["binder_len"]), float(r["plddt_binder"]),
              float(r["pae_interaction"])))
print("arm 3: %d designs\n" % len(d))
print("=== pLDDT pass rate (>92) by binder length ===")
print("  %-14s %6s %7s %8s %9s" % ("length","n","pass","rate","mean pLDDT"))
bins=[(50,70),(70,80),(80,90),(90,100),(100,120)]
for lo,hi in bins:
    v=[x for x in d if lo<=x[0]<hi]
    if not v: continue
    p=sum(1 for x in v if x[1]>92)
    print("  %-14s %6d %7d %7.1f%% %9.1f" % ("%d-%d"%(lo,hi-1),len(v),p,100*p/len(v),
          statistics.mean(x[1] for x in v)))
# Pearson r, length vs pLDDT
n=len(d); mx=statistics.mean(x[0] for x in d); my=statistics.mean(x[1] for x in d)
sx=(sum((x[0]-mx)**2 for x in d)/n)**.5; sy=(sum((x[1]-my)**2 for x in d)/n)**.5
r=sum((x[0]-mx)*(x[1]-my) for x in d)/(n*sx*sy)
print("\n  Pearson r (binder length vs pLDDT) = %+.3f" % r)
print("  arm 1 reference: 128 aa fixed, 61.7%% pass, n=60")
long_=[x for x in d if x[0]>=100]
p=sum(1 for x in long_ if x[1]>92)
print("\n=== arm 3's longest binders (>=100 aa), closest to arm 1's 128 ===")
print("  n=%d  pass=%d  %.1f%%   mean pLDDT %.1f" % (len(long_),p,100*p/len(long_),
      statistics.mean(x[1] for x in long_)))
