BLUE="#2D6A9F"; GREY="#C4C9C6"; MID="#8A908F"; RUST="#B4653F"
P=[]
def r(x,y,w,h,fill,rx=3,op=1.0):
    P.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}"'
             + (f' fill-opacity="{op}"' if op!=1.0 else '') + '/>')
def t(x,y,s,size=20,fill="#4a5568",weight=400,anchor="start"):
    P.append(f'<text x="{x}" y="{y}" font-family="IBM Plex Sans, Arial, sans-serif" '
             f'font-size="{size}" font-weight="{weight}" fill="{fill}" text-anchor="{anchor}">{s}</text>')

Y0=42; RW,G,C,ROWS=33,5,12,6
SIG={3,7,14,19,26,31,38,45,52,57,63}
for i in range(C*ROWS):
    r(G*0+(i%C)*(RW+G), Y0+(i//C)*(RW+G), RW, RW, BLUE if i in SIG else GREY, 3,
      1.0 if i in SIG else 0.85)
W1=C*(RW+G)-G; H1=ROWS*(RW+G)-G
t(0,Y0-14,"STAGE 1 · SHALLOW",20,RUST,600)

AX=W1+24; CY=Y0+H1/2
P.append(f'<path d="M {AX} {CY} L {AX+56} {CY}" stroke="{MID}" stroke-width="3" fill="none"/>')
P.append(f'<path d="M {AX+56} {CY} l -12 -8 v 16 z" fill="{MID}"/>')
t(AX+28,CY-18,"keep",18,MID,600,"middle")
t(AX+28,CY+34,"signal",18,MID,600,"middle")

X2=AX+86; COLW,COLG,BH,BG,NB=62,24,20,5,11
t(X2,Y0-14,"STAGE 2 · DEEP",20,RUST,600)
for c in range(4):
    cx=X2+c*(COLW+COLG)
    for k in range(NB):
        r(cx, Y0+k*(BH+BG), COLW, BH, BLUE, 3, 0.50+0.045*k)
W2=4*(COLW+COLG)-COLG
WIDTH=X2+W2; HEIGHT=Y0+max(H1, NB*(BH+BG)-BG)
svg=(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {WIDTH} {HEIGHT}" '
     f'style="width:{WIDTH}px; height:{HEIGHT}px" aria-label="Two-stage screen: stage one runs about '
     f'four designs against every donor scaffold and most show nothing; stage two runs many designs '
     f'against only the donors that scored.">' + "\n".join(P) + '</svg>')
open("donor_svg.txt","w").write(svg)
open("donor_dims.txt","w").write(f"{W1} {AX+86-W1} {W2} {WIDTH}")
print("svg",WIDTH,"x",HEIGHT," stage1",W1," gap",X2-W1," stage2",W2)
