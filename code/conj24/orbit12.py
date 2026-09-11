# Exact projective orbit of (S1,S5) for Q12 under the Weil representation (level 96), in Z[zeta_192]
import numpy as np, cmath, math, pickle
from math import gcd
D=64; M=192   # Phi_192(x) = x^64 - x^32 + 1
def red(v):   # reduce polynomial coefficient list (any length) mod x^64 - x^32 + 1
    v=list(v)
    for i in range(len(v)-1,D-1,-1):
        c=v[i]
        if c:
            v[i]=0; v[i-32]+=c; v[i-64]-=c      # x^i = x^{i-32} - x^{i-64}
    return np.array((v+[0]*D)[:D],dtype=object)
def mono(e):
    e%=M; v=[0]*(e+1); v[e]=1; return red(v)
MONO=[mono(e) for e in range(M)]
def mulmono(a,e):
    e%=M
    if e==0: return a.copy()
    v=[0]*(D+e); 
    for i,c in enumerate(a):
        if c: v[i+e]+=c
    return red(v)
def mul(a,b):
    v=[0]*(2*D)
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y: v[i+j]+=x*y
    return red(v)
Z=np.zeros(D,dtype=object); ONE=MONO[0]
NC=48
def actT(v): return [mulmono(v[c],2*c*c) for c in range(NC)]
def actS(v):
    out=[Z.copy() for _ in range(NC)]
    for c in range(NC):
        if any(v[c]):
            for d in range(NC): out[d]=out[d]+mulmono(v[c],4*c*d)
    return out
def content(p):
    g=0
    for v in p:
        for a in v:
            for x in a:
                if x: g=gcd(g,int(x))
    return g
def normalize(p):
    g=content(p); return tuple([a//g for a in v] for v in p) if g>1 else p
Z192=[cmath.exp(2j*math.pi*i/M) for i in range(D)]
def cnum(a): return sum(int(x)*Z192[i] for i,x in enumerate(a) if x)
def key(p):
    vals=[cnum(a) for v in p for a in v]
    i=next(i for i,x in enumerate(vals) if abs(x)>1e-9)
    return tuple((round((x/vals[i]).real,6)+0.0,round((x/vals[i]).imag,6)+0.0) for x in vals)
def prop(p,q):
    P=[a for v in p for a in v]; Q=[a for v in q for a in v]
    i=next(i for i,a in enumerate(P) if any(a))
    return all(np.array_equal(mul(P[i],Q[j]),mul(P[j],Q[i])) for j in range(len(P)))
def vec(d):
    v=[Z.copy() for _ in range(NC)]
    for c,s in d.items():
        v[c%NC]=v[c%NC]+s*ONE; v[(-c)%NC]=v[(-c)%NC]+s*ONE
    return v
v1=vec({10:1,14:-1}); v5=vec({2:1,22:-1})
if __name__=='__main__':
    start=normalize((v1,v5)); orbit=[start]; idx={key(start):0}; edges=[]; q=[0]
    while q:
        i=q.pop(); p=orbit[i]
        for nm,f in (('S',actS),('T',actT)):
            np_=normalize((f(p[0]),f(p[1]))); k=key(np_)
            if k in idx:
                j=idx[k]; assert prop(np_,orbit[j]),'collision'
            else:
                j=len(orbit); orbit.append(np_); idx[k]=j; q.append(j)
            edges.append((i,nm,j))
    print('orbit size',len(orbit),'edges',len(edges))
    pickle.dump((orbit,edges),open('orbit12.pkl','wb'))
