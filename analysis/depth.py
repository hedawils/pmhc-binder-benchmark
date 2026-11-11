import sys, re, collections
rows=[]
with open(sys.argv[1]) as fh:
    hdr=fh.readline().split()[1:]
    for line in fh:
        f=line.split()
        if f and f[0]=="SCORE:": rows.append(dict(zip(hdr,f[1:])))
by=collections.defaultdict(list)
for r in rows:
    m=re.match(r"(liuscaf\d+_0)_dldesign_(\d+)", r["description"])
    if not m: continue
    ok = (float(r["pae_interaction"])<7.0 and float(r["plddt_binder"])>92.0
          and float(r["binder_aligned_rmsd"])<2.0)
    by[m.group(1)].append((int(m.group(2)), ok))
print("arm 3 backbone hit rate at matched sequence depth")
print("  (arm 5 gets 4 seqs/backbone; arms 1b and 3 get 6)\n")
print("  %-10s %8s %8s %8s" % ("depth k", "bbones", "hits", "rate"))
for k in (1,2,3,4,5,6):
    elig=[v for v in by.values() if len(v)>=1]
    hits=0; n=0
    for v in by.values():
        sub=[ok for sid,ok in sorted(v) if sid < k]
        if not sub: continue
        n+=1; hits+= any(sub)
    print("  %-10d %8d %8d %7.1f%%" % (k,n,hits,100*hits/n if n else 0))
full=[v for v in by.values()]
print("\n  scored sequences per backbone: min %d  median %d  max %d"
      % (min(len(v) for v in full),
         sorted(len(v) for v in full)[len(full)//2],
         max(len(v) for v in full)))
