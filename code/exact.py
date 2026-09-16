# exact Gauss-sum data in Z[zeta_{16k}], reduction mod Phi_{16k}(x)=x^{8k}+1 (16k power of 2)
from fractions import Fraction
import cmath, math
def psi(a,j):
    r=a-4
    if (j-r)%8: return 0
    return -1 if ((j-r)//8)%2 else 1
def G(a,h,k,l):
    M=16*k; half=M//2
    v=[0]*half
    for s in range(M):
        p=psi(a,s)
        if not p: continue
        e=(h*s*s+l*s)%M
        if e>=half: v[e-half]-=p
        else: v[e]+=p
    return tuple(v)
def cnum(v,k):
    M=16*k
    return sum(c*cmath.exp(2j*math.pi*i/M) for i,c in enumerate(v))
def iszero(v): return all(c==0 for c in v)
def mul_root(v,t,k):  # multiply by zeta^t
    M=16*k; half=M//2; out=[0]*half
    for i,c in enumerate(v):
        if c:
            e=(i+t)%M
            if e>=half: out[e-half]-=c
            else: out[e]+=c
    return tuple(out)
def add(u,v): return tuple(a+b for a,b in zip(u,v))
def sub(u,v): return tuple(a-b for a,b in zip(u,v))
def Ghat(a,h,k,e):
    l=int(round(math.isqrt(e)))
    if l*l!=e: return None
    v=G(a,h,k,l)
    if l: v=add(v,G(a,h,k,-l))
    return v
for k in (8,16):
    for h in range(1,k,2):
        sup={}
        for a in (1,3):
            sup[a]=[l*l for l in range(0,3*k) if not iszero(Ghat(a,h,k,l*l))][:4]
        e1,e3=sup[1][0],sup[3][0]
        g1=Ghat(1,h,k,e1); g3=Ghat(3,h,k,e3)
        # find t with g1 == zeta^t g3 exactly
        M=16*k
        ts=[t for t in range(M) if mul_root(g3,t,k)==g1]
        print('h/k=%d/%d supp1'%(h,k),sup[1],'supp3',sup[3],'e3-e1=',e3-e1,' g1=zeta_%d^t*g3, t='%M,ts, ' |g1|=%.6f'%abs(cnum(g1,k)))
print('--- correction ratios rho = ghat3(100)/ghat3(36), and ghat1(196)/ghat1(4)')
for k in (8,16):
    M=16*k
    for h in range(1,k,2):
        for a,(e0,e1_) in ((3,(36,100)),(1,(4,196))):
            g0=Ghat(a,h,k,e0); gn=Ghat(a,h,k,e1_)
            if iszero(g0): continue
            ts=[t for t in range(M) if mul_root(g0,t,k)==gn]
            ts2=[t for t in range(M) if mul_root(g0,t,k)==tuple(-x for x in gn)]
            print(h,k,'a=%d'%a,'gn = zeta^t g0 t=',ts,' or -zeta^t: ',ts2)
