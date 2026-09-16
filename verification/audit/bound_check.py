# Empirical audit: the explicit expansions + certified error/tail bounds must bracket the exact coefficients.
import sys, mpmath as mp
sys.path.insert(0,'/home/claude/repo/code'); sys.path.insert(0,'/home/claude/repo/code/conj20')
mp.mp.dps=60
def exact(eps,N,d):
    b=[0]*(N+1)
    for dd in range(1,N+1):
        e=eps(dd)
        if e:
            for m in range(dd,N+1,dd): b[m]-=e*dd
    f=[mp.mpf(1)]+[mp.mpf(0)]*N; d=mp.mpf(d)
    for n in range(1,N+1): f[n]=d*mp.fsum(b[k]*f[n-k] for k in range(1,n+1))/n
    return f
# ---- Q8, delta in (2,4]: main k=8 + corr k=8 + main k=16 exact; rest bounded by T16corr + T>=24 + E
from flint import arb
import consts
def q8(d,n):
    D=mp.mpf(d); y=mp.sqrt(2*n+D); N=int(mp.ceil(mp.sqrt(4*mp.pi*n+2*mp.pi*D)))+1
    P=lambda k,B: (2*mp.pi/k)*mp.sqrt(B)/y*mp.besseli(1,(2*mp.pi/k)*mp.sqrt(B)*y) if B>0 else 0
    M=2*mp.cos(3*mp.pi*n/4)*P(8,D)-2*D*mp.cos((3*n-1)*mp.pi/4)*P(8,D-2)+4*mp.cos(mp.pi*n/8)*mp.cos(mp.pi*(n+D)/2)*P(16,D)
    T=4*D*P(16,D-2)+sum((k/2)*(P(k,D)+D*P(k,D-2)) for k in range(24,N+1,4))
    E=mp.mpf(consts.Econst(arb(str(d)),N).mid().str(20,radius=False))
    return M,T+E
e8=lambda d:(1 if d%8 in (1,7) else (-1 if d%8 in (3,5) else 0))
for d in ('2.7','3.2','3.9'):
    f=exact(e8,420,d); worst=0
    for n in range(60,421):
        M,B=q8(d,n); worst=max(worst,abs(f[n]-M)/B)
    print('Q8 delta=%s: max |c-explicit| / (tail+E) over 60<=n<=420 = %s  (must be <=1)'%(d,mp.nstr(worst,4)))
