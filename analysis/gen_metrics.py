PA="#2d6a9f"; FR="#6a7179"; INK="#16202b"; MUT="#4a5568"; DIM="#8a908f"; RULE="#d4d6d2"; RUST="#b4653f"
FONT="'IBM Plex Sans', Arial, sans-serif"
CELL=f"padding:5px 10px; text-align:right; font-size:21px; line-height:1.25; color:{INK}"
MET =f"padding:5px 10px; text-align:left; width:330px; font-size:21px; line-height:1.25; color:{MUT}"
SEP =f"border-left:2px solid {RULE}"
PCOL=f"padding:5px 10px; text-align:right; font-size:19px; color:{DIM}; width:150px"

def sub(label,n,color):
    return (f'<th style="padding:5px 10px; text-align:right; font-size:19px; font-weight:600; '
            f'color:{color}; border-bottom:2px solid {RULE}; padding-bottom:5px{";"+SEP if n=="SEP" else ""}">'
            f'{label}</th>')

def hdr(a,b,na,nb,first_sep):
    def cell(lbl,n,col,sep):
        return (f'<th style="padding:5px 10px; text-align:right; font-size:19px; font-weight:600; '
                f'color:{col}; border-bottom:2px solid {RULE}{"; "+SEP if sep else ""}">'
                f'{lbl}<br><span style="font-weight:400; color:{DIM}">n = {n}</span></th>')
    return (cell("partial",na,PA,first_sep) + cell("free",nb,FR,False)
            + f'<th style="padding:5px 10px; text-align:right; font-size:19px; font-weight:600; '
              f'color:{DIM}; border-bottom:2px solid {RULE}">p</th>')

def row(label, ny, kf):
    """ny/kf = (partial, free, p, winner)  winner in {'p','f',None}"""
    out=[f'<td style="{MET}">{label}</td>']
    for i,(pv,fv,pval,win) in enumerate((ny,kf)):
        s = SEP if i==1 else ""
        pw = "font-weight:600; color:"+PA if win=="p" else ""
        fw = "font-weight:600; color:"+FR if win=="f" else ""
        out.append(f'<td style="{CELL}; {s}; {pw}">{pv}</td>')
        out.append(f'<td style="{CELL}; {fw}">{fv}</td>')
        out.append(f'<td style="{PCOL}">{pval}</td>')
    return "<tr>"+"".join(out)+"</tr>"

def band(txt):
    return (f'<tr><td colspan="7" style="padding:11px 10px 1px; text-align:left; font-size:18px; '
            f'font-weight:600; letter-spacing:.07em; text-transform:uppercase; color:{RUST}">{txt}</td></tr>')

T=f'border-collapse:collapse; width:1664px; font-family:{FONT}'
A=[f'<table style="{T}">',
   '<tr><th style="width:330px"></th>'
   f'<th colspan="3" style="padding:5px 10px; text-align:center; font-size:19px; font-weight:600; '
   f'letter-spacing:.06em; text-transform:uppercase; color:{DIM}">NY-ESO-1 · HLA-A*02:01</th>'
   f'<th colspan="3" style="padding:5px 10px; text-align:center; font-size:19px; font-weight:600; '
   f'letter-spacing:.06em; text-transform:uppercase; color:{DIM}; {SEP}">SARS-CoV-2 KF9 · HLA-C*12:02</th></tr>',
   f'<tr><th style="padding:5px 10px; text-align:left; font-size:19px; font-weight:600; color:{MUT}; '
   f'border-bottom:2px solid {RULE}">all designs</th>' + hdr("partial","free",852,800,False).replace(
       '<th style="padding:5px 10px; text-align:right; font-size:19px; font-weight:600; color:'+PA+'; border-bottom:2px solid '+RULE+'">',
       '<th style="padding:5px 10px; text-align:right; font-size:19px; font-weight:600; color:'+PA+'; border-bottom:2px solid '+RULE+'">')
   + hdr("partial","free",588,600,True) + '</tr>']
# the hdr() helper above emits 3 cells each; build header manually for the sep
A=[f'<table style="{T}">']
A.append('<tr><th style="width:330px"></th>'
   f'<th colspan="3" style="padding:5px 10px; text-align:center; font-size:19px; font-weight:600; letter-spacing:.06em; text-transform:uppercase; color:{DIM}">NY-ESO-1 · HLA-A*02:01</th>'
   f'<th colspan="3" style="padding:5px 10px; text-align:center; font-size:19px; font-weight:600; letter-spacing:.06em; text-transform:uppercase; color:{DIM}; {SEP}">SARS-CoV-2 KF9 · HLA-C*12:02</th></tr>')
def headrow(left,na,nb,nc,nd):
    def c(lbl,n,col,sep):
        return (f'<th style="padding:5px 10px; text-align:right; font-size:19px; font-weight:600; color:{col}; '
                f'border-bottom:2px solid {RULE}{"; "+SEP if sep else ""}">{lbl}'
                + (f'<br><span style="font-weight:400; color:{DIM}">n = {n}</span>' if n else '') + '</th>')
    p=lambda sep: (f'<th style="padding:5px 10px; text-align:right; font-size:19px; font-weight:600; color:{DIM}; '
                   f'border-bottom:2px solid {RULE}">p</th>')
    return ('<tr>'
            f'<th style="padding:5px 10px; text-align:left; font-size:19px; font-weight:600; color:{MUT}; border-bottom:2px solid {RULE}">{left}</th>'
            + c("partial",na,PA,False)+c("free",nb,FR,False)+p(False)
            + c("partial",nc,PA,True)+c("free",nd,FR,False)+p(False)+'</tr>')

A.append(headrow("all designs",852,800,588,600))
A.append(band("fold"))
A.append(row("binder pLDDT", ("90.91","92.15","7.6e-06","f"), ("90.79","92.96","7.1e-10","f")))
A.append(row("scRMSD (Å)",   ("0.726","0.784","0.008","p"),   ("0.784","0.764","ns",None)))
A.append(band("interface"))
A.append(row("designs touching the peptide", ("43.4%","30.0%","1.6e-08","p"), ("43.7%","43.5%","ns",None)))
A.append(row("buried peptide Å² <sup>†</sup>", ("297","46","9.5e-13","p"), ("107","103","ns",None)))
A.append(row("interface pAE <sup>†</sup>",    ("14.8","25.0","1.5e-14","p"), ("20.3","20.9","ns",None)))
A.append(row("buried MHC Å²",                 ("63","25","1.7e-17","p"),     ("89","64","ns",None)))
A.append("</table>")

B=[f'<table style="{T}; margin-top:6px">', headrow("four-gate survivors",67,15,21,25)]
B.append(row("binder pLDDT",  ("93.35","94.41","0.0003","f"), ("94.24","94.80","ns",None)))
B.append(row("interface pAE", ("4.54","5.37","ns",None),      ("4.56","5.22","0.015","p")))
B.append(row("peptide focus", ("0.319","0.280","0.006","p"),  ("0.220","0.222","ns",None)))
B.append(row("buried MHC Å²", ("769","905","0.001","f"),      ("865","863","ns",None)))
B.append("</table>")
open("metrics_tables.html","w").write("\n".join(A)+"\n"+"\n".join(B))
print("tables written, no <style> block")
