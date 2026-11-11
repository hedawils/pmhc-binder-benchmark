import sys, collections
A3={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
    "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
    "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
# minimal mmCIF atom_site reader: find the loop header, then read rows
path=sys.argv[1]
cols=[]; rows=[]; inloop=False
for line in open(path, errors="replace"):
    s=line.strip()
    if s.startswith("_atom_site."):
        cols.append(s.split(".",1)[1]); inloop=True; continue
    if inloop:
        if s.startswith("#") or s.startswith("loop_") or s=="":
            if rows: break
            continue
        f=s.split()
        if len(f)==len(cols): rows.append(f)
i={c:n for n,c in enumerate(cols)}
need=("group_PDB","label_atom_id","label_comp_id","auth_asym_id","auth_seq_id","pdbx_PDB_model_num")
if not all(k in i for k in need): sys.exit("unexpected cif columns: %s" % cols[:12])
ch=collections.OrderedDict()
for r in rows:
    if r[i["group_PDB"]]!="ATOM": continue
    if r[i["pdbx_PDB_model_num"]]!="1": continue
    if r[i["label_atom_id"]]!="CA": continue
    ch.setdefault(r[i["auth_asym_id"]],[]).append((r[i["auth_seq_id"]], A3.get(r[i["label_comp_id"]],"X")))
print("\n=== 6ULR chains (model 1, CA only) ===")
for c,v in ch.items():
    seq="".join(x for _,x in v)
    print("  chain %-3s n=%-4d  first=%-5s last=%-5s" % (c,len(v),v[0][0],v[-1][0]))
    print("        %s" % (seq if len(seq)<=60 else seq[:60]+"..."))
