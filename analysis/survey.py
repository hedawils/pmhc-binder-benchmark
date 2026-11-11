import sys
from pathlib import Path
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1 as T

def seq(ch):
    return "".join(T.get(r.get_resname(), "X") for r in ch if r.id[0] == " ")

for p in sys.argv[1:]:
    path = Path(p)
    model = MMCIFParser(QUIET=True).get_structure(path.stem, path)[0]
    # resolution / title from the raw cif
    txt = path.read_text()
    res = next((l.split()[-1] for l in txt.splitlines()
                if l.startswith("_refine.ls_d_res_high")), "?")
    print(f"\n=== {path.stem.upper()}  resolution {res} ===")
    for ch in model:
        s = seq(ch)
        if not s:
            continue
        rs = [r.id[1] for r in ch if r.id[0] == " "]
        tag = "  <-- 9mer" if len(s) == 9 else ("  <-- 10mer" if len(s) == 10 else "")
        print(f"  {ch.id}: {len(s):>4} res  {rs[0]}..{rs[-1]}  {s[:50]}{tag}")
