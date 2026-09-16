import sys, mpmath as mp, math
sys.path.insert(0,'/home/claude/repo/code/conj20')
from bound_check import exact
import g7lib
from flint import arb
mp.mp.dps=60
def fr(x): return mp.mpf(x.numerator)/x.denominator
def g7model(d,n):
    D=mp.mpf(d); C=2*n-D/2; B=D/2; B2=D/2-2; N=int(mp.ceil(mp.sqrt(4*mp.pi*n)))+2
    P=lambda k,B: (2*mp.pi/k)*mp.sqrt(B/C)*mp.besseli(1,(2*mp.pi/k)*mp.sqrt(B*C)) if B>0 else 0
    M=0
    for k in (7,14,21):
        for h,sh,hp in g7lib.shifts(k):
            ph=-mp.pi*D*fr(sh)-2*mp.pi*h*n/k
            M+=mp.cos(ph)*P(k,B)
            if D>4: M+=-D*mp.cos(ph-2*mp.pi*hp/k)*P(k,B2)
    T=sum(g7lib.math.gcd(1,1)*0 for _ in [])
    T=sum(sum(1 for h in range(k) if math.gcd(h,k)==1)*(P(k,B)+(D*P(k,B2) if D>4 else 0)) for k in range(28,N+1,7))
    E=mp.mpf(g7lib.Econ(arb(str(d)),N).mid().str(20,radius=False))
    return M,T+E
e7=lambda d:(1 if d%7 else 0)
for d in ('2.5','3.5','4.6','4.8735'):
    f=exact(e7,420,d); worst=0
    for n in range(60,421):
        M,Bd=g7model(d,n); worst=max(worst,abs(f[n]-M)/Bd)
    print('G7 delta=%s: max |c-explicit|/(tail+E), 60<=n<=420 = %s (must be <=1)'%(d,mp.nstr(worst,4)))
