import pickle, numpy as np
from orbit12 import *
orbit,edges=pickle.load(open('orbit12.pkl','rb'))
def conj(a):
    out=Z.copy()
    for i,x in enumerate(a):
        if x: out=out+x*mulmono(ONE,-i)
    return out
Tn={i:j for i,n,j in edges if n=='T'}
seen=set(); cyc=[]
for i in range(len(orbit)):
    if i in seen: continue
    c=[i]; seen.add(i); j=Tn[i]
    while j!=i: c.append(j); seen.add(j); j=Tn[j]
    cyc.append(c)
print('T-cycles:',len(cyc))
for c in cyc:
    w1,w5=orbit[c[0]]
    J1=min(k for k in range(25) if any(w1[k])); J5=min(k for k in range(25) if any(w5[k]))
    s1=[k for k in range(25) if any(w1[k])]; s5=[k for k in range(25) if any(w5[k])]
    mods=set(tuple(mul(v[k],conj(v[k]))) for v,s in ((w1,s1),(w5,s5)) for k in s)
    g='GROW(+)' if J1<J5 else ('GROW(-)' if J1>J5 else 'flat')
    print('len %2d  J1=%2d J5=%2d  exp=(J1^2-J5^2)/96=%6.3f  %-8s supp1=%s supp5=%s  all moduli equal: %s'%(len(c),J1,J5,(J1*J1-J5*J5)/96,g,s1,s5,len(mods)==1))
