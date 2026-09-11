# Explicit Rademacher-type model for G7^delta: main terms (Schlosser-Zhou Thm 3) + correction terms (delta>4).
import mpmath as mp, math
from fractions import Fraction as F
mp.mp.dps=50
def dedekind(h,k):
    tot=F(0)
    for j in range(1,k):
        x=F(j,k); y=F(j*h,k)
        fx=x-math.floor(x)-F(1,2) if x.denominator!=1 else 0
        fy=y-math.floor(y)-F(1,2) if y.denominator!=1 else 0
        tot+=fx*fy
    return tot
SC={}
def omega_exp(h,k):   # s(h,k) - s(h,k/7) for 7|k
    key=(h,k)
    if key not in SC: SC[key]=dedekind(h,k)-dedekind(h,k//7)
    return SC[key]
def model(n,d,KMAX=70,corr=True,hsign=1):
    d=mp.mpf(d); tot=mp.mpf(0); C=2*n-d/2
    for k in range(7,KMAX+1,7):
        for h in range(k):
            if math.gcd(h,k)!=1: continue
            ph=mp.expjpi(-d*omega_exp(h,k))*mp.expjpi(-2*mp.mpf(h*n)/k)
            B=d/2
            tot+=ph*(2*mp.pi/k)*mp.sqrt(B/C)*mp.besseli(1,(2*mp.pi/k)*mp.sqrt(B*C))
            if corr and d>4:
                hp=pow(h,-1,k)*hsign
                B2=d/2-2
                tot+=ph*(-d)*mp.expjpi(2*mp.mpf(hp)/k)*(2*mp.pi/k)*mp.sqrt(B2/C)*mp.besseli(1,(2*mp.pi/k)*mp.sqrt(B2*C))
    return tot
if __name__=='__main__':
    import sys
    d=mp.mpf(sys.argv[1]); N=400
    b=[0]*(N+1)
    for dd in range(1,N+1):
        if dd%7:
            for m in range(dd,N+1,dd): b[m]-=dd
    f=[mp.mpf(1)]+[0]*N
    for n in range(1,N+1): f[n]=d*mp.fsum(b[k]*f[n-k] for k in range(1,n+1))/n
    for n in range(N-6,N+1):
        m1=model(n,d,corr=False); m2=model(n,d,hsign=1); m3=model(n,d,hsign=-1)
        print(n,n%7,mp.nstr(f[n],12),' main only',mp.nstr(mp.re(m1),12),' +corr(h\')',mp.nstr(mp.re(m2),12),' +corr(-h\')',mp.nstr(mp.re(m3),12))
