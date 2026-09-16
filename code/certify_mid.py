# Step 4: certified check of Prop 5.1, NLO<=n<=NHI, delta in [8/3,4], self-contained constants.
import sys, time
from flint import arb, ctx
ctx.prec=100
from consts import Econst, PI, R2
NLO,NHI,S=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
I1=lambda x:x.bessel_i(1)
lo=arb(8)/3
pts=[lo+(arb(3)-lo)*i/S for i in range(S)]+[arb(3)+arb(i)/S for i in range(S+1)]
subs=[(pts[i],pts[i+1]) for i in range(len(pts)-1)]
def cosmax(a,b):
    u,v=abs(arb.cos_pi(a/2)),abs(arb.cos_pi(b/2))
    return u if u>v else (v if v>u else u+v)
def Plow(k,Bl,yl,yh): return (2*PI/k)*Bl.sqrt()/yh*I1((2*PI/k)*Bl.sqrt()*yl)
def Pup(k,Bh,yl,yh):  return (2*PI/k)*Bh.sqrt()/yl*I1((2*PI/k)*Bh.sqrt()*yh)
Ecache={}
fails=[];t0=time.time()
for n in range(NLO,NHI+1):
    for a,b in subs:
        yl=(2*n+a).sqrt(); yh=(2*n+b).sqrt()
        N=int(float((4*PI*n+2*PI*b).sqrt().upper()))+1
        key=(float(b),N)
        if key not in Ecache: Ecache[key]=Econst(b,N)
        E=Ecache[key]
        nk4=max(0,N//4-5)                       # k=24,28,...,<=N with 4|k
        T24=nk4*PI*(b.sqrt()*I1((2*PI/24)*b.sqrt()*yh)+b*(b-2).sqrt()*I1((2*PI/24)*(b-2).sqrt()*yh))/yl
        P16c=4*b*Pup(16,b-2,yl,yh)
        if n%8 in (2,6):
            T16=2*R2*cosmax(a,b)*Pup(16,b,yl,yh)+P16c
            m=R2*a*Plow(8,a-2,yl,yh)-T16-T24-E
        else:
            T16=4*Pup(16,b,yl,yh)+P16c
            m=R2*Plow(8,a,yl,yh)-2*b*Pup(8,b-2,yl,yh)-T16-T24-E
        if not (m>0): fails.append((n,float(a),float(b)))
print('n=%d..%d, %d subintervals: %.1fs, failures %d'%(NLO,NHI,len(subs),time.time()-t0,len(fails)), fails[:12])
