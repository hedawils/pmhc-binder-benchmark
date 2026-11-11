import json, sys, urllib.request, collections
SC = dict(G=0,A=1,S=2,C=2,T=3,V=3,P=3,L=4,I=4,N=4,D=4,M=4,Q=5,K=5,E=5,H=6,F=7,R=7,Y=8,W=10)
IDS = sys.argv[1:]
def get(u):
    try:
        with urllib.request.urlopen(u, timeout=20) as r: return json.load(r)
    except Exception: return None
rows=[]
for pid in IDS:
    e = get(f"https://data.rcsb.org/rest/v1/core/entry/{pid}")
    if not e: continue
    title = (e.get("struct",{}) or {}).get("title","?")
    ents = (e.get("rcsb_entry_container_identifiers",{}) or {}).get("polymer_entity_ids",[])
    res = (e.get("rcsb_entry_info",{}) or {}).get("resolution_combined",[None])
    res = res[0] if res else None
    pep=None; nchain=0; allele="?"
    for en in ents:
        d = get(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pid}/{en}")
        if not d: continue
        seq=((d.get("entity_poly",{}) or {}).get("pdbx_seq_one_letter_code_can","") or "").replace("\n","")
        desc=((d.get("rcsb_polymer_entity",{}) or {}).get("pdbx_description") or "")
        ch=((d.get("rcsb_polymer_entity_container_identifiers",{}) or {}).get("auth_asym_ids") or [])
        nchain+=len(ch)
        if 8<=len(seq)<=11 and set(seq)<=set(SC): pep=seq
        if "Cw" in desc or "C-" in desc or "class I" in desc.lower(): allele=desc[:30]
    if pep:
        rows.append((sum(SC[c] for c in pep), pep, pid, nchain, res, title[:58]))
rows.sort(reverse=True)
print("  %-6s %-11s %6s %7s %6s  %s" % ("entry","peptide","sc-atm","chains","res","title"))
for n,pep,pid,nc,res,title in rows:
    print("  %-6s %-11s %6d %7d %6s  %s" % (pid,pep,n,nc, ("%.1f"%res) if res else "—", title))
print()
print("  NY-ESO-1 SLLMWITQV = 39 sc-atoms, static peptide_focus 0.354")
