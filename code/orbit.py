# Exact projective orbit of the pair (v1,v3) under right multiplication by E (S) and D (T),
# vectors symmetric (v_c = v_{-c}) over Z[zeta_128], stored as integer arrays mod x^64+1.
import numpy as np, cmath, math, pickle, sys
from math import gcd
from functools import reduce
D_=64
def zmul_mono(a,e):           # a * zeta^e
    e%=128
    out=np.zeros(D_,dtype=object)
    for i in range(D_):
        if a[i]:
            j=i+e
            if j>=128: j-=128
            if j<D_: out[j]+=a[i]
            else: out[j-D_]-=a[i]
    return out
def zmul(a,b):
    out=np.zeros(D_,dtype=object)
    for i in range(D_):
        if b[i]: out=out+b[i]*zmul_mono(a,i)
    return out
ZERO=np.zeros(D_,dtype=object)
def vec(dic):
    v=[ZERO.copy() for _ in range(32)]
    for c,x in dic.items():
        v[c%32]=v[c%32]+x*zmul_mono(np.array([1]+[0]*63,dtype=object),0)
    return v
def sym(v):  # symmetrize: w_c = v_c + v_{-c}  (so that sum_c v_c th_c = sum_c w_c th_c /... ) -> use canonical: w_c=w_{-c}
    return [v[c]+v[(-c)%32] if c not in (0,16) else 2*v[c] for c in range(32)]
def actT(v): return [zmul_mono(v[c],2*c*c) for c in range(32)]
def actS(v):
    out=[ZERO.copy() for _ in range(32)]
    for c in range(32):
        if any(v[c]):
            for d in range(32):
                out[d]=out[d]+zmul_mono(v[c],4*c*d)     # E_{cd}=zeta_32^{cd}=zeta_128^{4cd}
    return out
def content(pair):
    g=0
    for v in pair:
        for a in v:
            for x in a:
                if x: g=gcd(g,int(x))
    return g
def normalize(pair):
    g=content(pair)
    return tuple([a//g for a in v] for v in pair) if g>1 else pair
def cnum(a): return sum(int(x)*cmath.exp(2j*math.pi*i/128) for i,x in enumerate(a) if x)
def key(pair):
    vals=[cnum(a) for v in pair for a in v]
    i=next(i for i,x in enumerate(vals) if abs(x)>1e-9)
    return tuple((round((x/vals[i]).real,6)+0.0,round((x/vals[i]).imag,6)+0.0) for x in vals)
def prop(p,q):   # exact proportionality test: all cross products equal
    P=[a for v in p for a in v]; Q=[a for v in q for a in v]
    i=next(i for i,a in enumerate(P) if any(a))
    for j in range(len(P)):
        if not np.array_equal(zmul(P[i],Q[j]),zmul(P[j],Q[i])): return False
    return True
one=np.array([1]+[0]*63,dtype=object)
v1=[ZERO.copy() for _ in range(32)]; v3=[ZERO.copy() for _ in range(32)]
for c,s in ((6,1),(26,1),(10,-1),(22,-1)): v1[c]=v1[c]+s*one
for c,s in ((2,1),(30,1),(14,-1),(18,-1)): v3[c]=v3[c]+s*one
start=normalize((v1,v3))
orbit=[start]; idx={key(start):0}; edges=[]
q=[0]
while q:
    i=q.pop()
    p=orbit[i]
    for name,f in (('S',actS),('T',actT)):
        np_=normalize((f(p[0]),f(p[1])))
        k=key(np_)
        if k in idx:
            j=idx[k]
            assert prop(np_,orbit[j]), 'key collision without exact proportionality'
        else:
            j=len(orbit); orbit.append(np_); idx[k]=j; q.append(j)
        edges.append((i,name,j))
print('orbit size',len(orbit),' edges',len(edges))
pickle.dump((orbit,edges),open('orbit.pkl','wb'))
