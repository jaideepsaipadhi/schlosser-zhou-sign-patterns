import sys, mpmath as mp, math
sys.path.insert(0,'/home/claude/repo/code/conj24'); sys.path.insert(0,'/home/claude/repo/code/conj23')
from bound_check import exact
from flint import arb
mp.mp.dps=60
# Q12, 2<=delta<=3 (criterion ingredients of q12_pos.py)
exec(open('/home/claude/repo/code/conj24/q12_pos.py').read().split('NLO,NHI,S=')[0].replace('import sys, time',''))
def q12model(d,n):
    D=mp.mpf(d); C=2*(n+D); N=int(mp.ceil(mp.sqrt(4*mp.pi*(n+D))))+2
    P_=lambda k,B:(2*mp.pi/k)*mp.sqrt(B/C)*mp.besseli(1,(2*mp.pi/k)*mp.sqrt(B*C)) if B>0 else 0
    M=2*mp.cos(5*mp.pi*n/6)*P_(12,2*D)-2*D*mp.cos(mp.pi*(5*n-1)/6)*P_(12,2*(D-1))+D*(D+1)*mp.cos(mp.pi*(5*n-2)/6)*P_(12,2*(D-2))
    M+=(2*mp.cos(5*mp.pi*n/12)+2*mp.cos(7*mp.pi*n/12))*P_(24,2*D)
    T=4*D*P_(24,2*(D-1))+2*D*(D+1)*P_(24,2*(D-2))
    T+=6*(P_(36,2*D)+D*P_(36,2*(D-1))+D*(D+1)/2*P_(36,2*(D-2)))+sum((k//2)*(P_(k,2*D)+D*P_(k,2*(D-1))+D*(D+1)/2*P_(k,2*(D-2))) for k in range(40,N+1,4))
    E=mp.mpf(Econ(arb(str(d)),N).mid().str(20,radius=False))
    return M,T+E
e12=lambda d:(1 if d%12 in (1,11) else (-1 if d%12 in (5,7) else 0))
for d in ('2.2','2.6','3'):
    f=exact(e12,400,d); w=max(abs(f[n]-q12model(d,n)[0])/q12model(d,n)[1] for n in range(60,401))
    print('Q12 delta=%s: max |c-explicit|/(tail+E) = %s'%(d,mp.nstr(w,4)))
# G11
import gplib
def fr(x): return mp.mpf(x.numerator)/x.denominator
def g11model(d,n):
    D=mp.mpf(d); C=2*n-5*D/6; B=5*D/6; N=int(mp.ceil(mp.sqrt(4*mp.pi*n)))+2
    P_=lambda k,B:(2*mp.pi/k)*mp.sqrt(B/C)*mp.besseli(1,(2*mp.pi/k)*mp.sqrt(B*C))
    M=sum(mp.cos(-mp.pi*D*fr(sh)-2*mp.pi*h*n/k)*P_(k,B) for k in (11,22) for h,sh,hp in gplib.shifts(k))
    T=sum(sum(1 for h in range(k) if math.gcd(h,k)==1)*P_(k,B) for k in range(33,N+1,11))
    E=mp.mpf(gplib.Econ(arb(str(d)),N).mid().str(20,radius=False))
    return M,T+E
e11=lambda d:(1 if d%11 else 0)
for d in ('1.76','1.9','2'):
    f=exact(e11,400,d); w=max(abs(f[n]-g11model(d,n)[0])/g11model(d,n)[1] for n in range(60,401))
    print('G11 delta=%s: max |c-explicit|/(tail+E) = %s'%(d,mp.nstr(w,4)))
