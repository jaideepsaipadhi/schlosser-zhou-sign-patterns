# Audit of the epsilon-difference bound (Q8, -1<delta<=-0.99): for n=3,7 mod 8,
# |c_delta(n) - [k=8 term] - [k=16 term]| / eps must be <= (k>=24 difference bounds + arc + chord + non-growing).
import sys, mpmath as mp
sys.path.insert(0,'/home/claude/review')
from bound_check import exact
from flint import arb
mp.mp.dps=60
src=open('/home/claude/repo/code/neg_certify.py').read()
exec(src.split('NLO,NHI=int(sys.argv[1])')[0].replace('import sys, time','import sys, time\nsys.argv=["x","0","0","0.99"]'))
e8=lambda d:(1 if d%8 in (1,7) else (-1 if d%8 in (3,5) else 0))
worst=0
for dstr in ('-0.995','-0.9999'):
    d=mp.mpf(dstr); s=-d; eps=1-s
    f=exact(e8,400,dstr)
    for n in range(280,401):
        if n%8 not in (3,7): continue
        y=mp.sqrt(2*n-s)
        Pm=lambda k,B:(2*mp.pi/k)*mp.sqrt(B)/y*mp.besseli(1,(2*mp.pi/k)*mp.sqrt(B)*y)
        main8=2*mp.cos(mp.pi*(n-s)/4)*Pm(8,s)
        A16=2*mp.cos(mp.pi*(n-s)/8)+2*mp.cos(mp.pi*(7*n+s)/8)
        diff=abs(f[n]-main8-A16*Pm(16,s))/eps
        # bound pieces from neg_certify (at n)
        N=int(float((4*PI*n).sqrt().upper()))+2; G=growing(N); yb=(2*n-arb(1)).sqrt()
        rest=arb(0)
        for k,c in G:
            if k<24: continue
            a=2*PI/k; Dk=arb('0.506')+1/(2*yb*yb)+arb('0.5026')*a*yb*(1+1/(a*arb('0.99').sqrt()*yb))
            rest+=c*((arb(k*k)/2+arb('0.01'))+Dk)*P(k,arb(1),n)
        arcs=sum((arb(c)/(k*k)*2*E1*PI.exp()*((arb(k*k)/2+arb('0.01'))*(PI/R2)*k/(N+1)+(2+PI/2)/arb('0.99')+(PI*PI/R2)/(k*(N+1))) for k,c in G),arb(0))
        chord=sum((arb(c)/(k*k)*(2*R2*k/(N+1))*E1*PI.exp()*((arb(k*k)/2+2*PI*N*N/arb(k*k)+arb('0.1'))*F1m+F1*ell) for k,c in G),arb(0))
        ng=2*R2*E1*KP*(KP.log()+arb('0.01')+arb('0.556')*N*N)
        bd=float((rest+arcs+chord+ng).mid())
        worst=max(worst,float(diff)/bd)
print('Q8 epsilon-difference audit (delta=-0.995,-0.9999; 280<=n<=400, n=3,7 mod 8): max ratio = %.3e (must be <=1)'%worst)
