import sys, os, csv, re, collections, statistics
A3={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
scafdir, man_path, sc_path = sys.argv[1], sys.argv[2], sys.argv[3]
tail={}
for fn in os.listdir(scafdir):
    if fn.endswith(".pdb"):
        s=[]
        for l in open(os.path.join(scafdir,fn),errors="replace"):
            if l.startswith("ATOM") and l[12:16].strip()=="CA" and l[21]=="B":
                s.append(A3.get(l[17:20].strip(),"X"))
        tail[fn[:-4]]="".join(s)
man={r["scaffold"]:r for r in csv.DictReader(open(man_path),delimiter="\t")}
rows=[]
with open(sc_path) as fh:
    hdr=fh.readline().split()[1:]
    for line in fh:
        f=line.split()
        if f and f[0]=="SCORE:": rows.append(dict(zip(hdr,f[1:])))
by=collections.defaultdict(list)
for r in rows:
    m=re.match(r"(scaf\d+)_0_dldesign_(\d+)", r["description"].replace("liu","",1))
    if not m: continue
    ok=(float(r["pae_interaction"])<7.0 and float(r["plddt_binder"])>92.0
        and float(r["binder_aligned_rmsd"])<2.0)
    by[m.group(1)].append((int(m.group(2)), ok, float(r["pae_interaction"])))

def cohort(s):
    m=man[s]
    if m["source_peptide"] not in ("-",""):
        return "%s / %s / %s" % (m["source_target"], m["source_allele"], m["source_peptide"])
    t=tail.get(s,"")
    return "untyped / ? / %s" % (t[-8:] if t else "?")

g=collections.defaultdict(lambda: dict(n=0,h4=0,h6=0,pae=[]))
for s,v in by.items():
    if s not in man or man[s]["accepted"]!="1": continue
    c=cohort(s); e=g[c]; e["n"]+=1
    e["h6"] += any(ok for _,ok,_ in v)
    e["h4"] += any(ok for sid,ok,_ in v if sid<4)
    e["pae"].append(min(p for _,_,p in v))

tot=dict(n=0,h4=0,h6=0)
print("\n=== arm 3 provenance of EVERY hit (%d designs scored) ===" % len(rows))
print("  %-40s %4s %6s %6s %8s %s" % ("cohort: target / allele / peptide","n","hits@6","hits@4","rate@6","median best pAE"))
for c,e in sorted(g.items(), key=lambda kv:-kv[1]["h6"]):
    tot["n"]+=e["n"]; tot["h4"]+=e["h4"]; tot["h6"]+=e["h6"]
    print("  %-40s %4d %6d %6d %7.1f%%   %.1f"
          % (c,e["n"],e["h6"],e["h4"],100*e["h6"]/e["n"],statistics.median(e["pae"])))
print("  %-40s %4d %6d %6d %7.1f%%" % ("TOTAL",tot["n"],tot["h6"],tot["h4"],100*tot["h6"]/tot["n"]))
print()
print("  hits accounted for: %d of %d  (depth 6);  %d (depth 4, matched to arm 5)"
      % (tot["h6"],tot["h6"],tot["h4"]))
print()
print("=== share of all hits contributed by each cohort ===")
for c,e in sorted(g.items(), key=lambda kv:-kv[1]["h6"]):
    if e["h6"]:
        print("  %-40s %2d hits = %4.1f%% of all hits  (from %4.1f%% of scaffolds)"
              % (c,e["h6"],100*e["h6"]/tot["h6"],100*e["n"]/tot["n"]))
