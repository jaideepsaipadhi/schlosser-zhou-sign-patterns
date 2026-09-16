import pickle, numpy as np
exec(open('orbit.py').read().split('one=np.array')[0])
orbit,edges=pickle.load(open('orbit.pkl','rb'))
def conj(a):
    out=np.zeros(D_,dtype=object)
    for i in range(D_):
        if a[i]: out=out+a[i]*zmul_mono(np.array([1]+[0]*63,dtype=object),-i)
    return out
def nz(a): return any(a)
def lead(v): return min(c for c in range(17) if nz(v[c]))
Tnext={i:j for i,n,j in edges if n=='T'}
# T-cycles
seen=set(); cycles=[]
for i in range(len(orbit)):
    if i in seen: continue
    cyc=[i]; seen.add(i); j=Tnext[i]
    while j!=i: cyc.append(j); seen.add(j); j=Tnext[j]
    cycles.append(cyc)
print('number of T-cycles (cusp types):',len(cycles))
summary={}
for cyc in cycles:
    w1,w3=orbit[cyc[0]]
    J1,J3=lead(w1),lead(w3)
    supp1=sorted(c for c in range(17) if nz(w1[c])); supp3=sorted(c for c in range(17) if nz(w3[c]))
    # exact modulus-squared values on supports
    m1={c:tuple(zmul(w1[c],conj(w1[c]))) for c in supp1}; m3={c:tuple(zmul(w3[c],conj(w3[c]))) for c in supp3}
    const1=len(set(m1.values()))==1; const3=len(set(m3.values()))==1
    same=(m1[supp1[0]]==m3[supp3[0]])
    tag='GROWING' if J1<J3 else ('zero-growth' if J1==J3 else 'decaying')
    print('cycle len %2d  J1=%2d J3=%2d  %-11s supp1=%s supp3=%s  |w| const on supp: %s %s  |w1|=|w3|: %s'%(len(cyc),J1,J3,tag,supp1,supp3,const1,const3,same))
