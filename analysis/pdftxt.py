import sys
from pypdf import PdfReader
r = PdfReader(sys.argv[1])
out = []
for i, p in enumerate(r.pages, 1):
    try:
        t = p.extract_text() or ""
    except Exception as e:
        t = f"[extract failed: {e}]"
    out.append(f"\n===== PAGE {i} =====\n{t}")
print(f"pages: {len(r.pages)}", file=sys.stderr)
open(sys.argv[2], "w").write("".join(out))
