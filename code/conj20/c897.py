# Exact c_897(delta) for G7^delta and certified isolation of its root delta_c in the bracket.
import time, pickle
from flint import fmpq_poly, fmpq, fmpz_poly, arb, ctx
t0=time.time()
NMAX=897
b=[0]*(NMAX+1)
for d in range(1,NMAX+1):
    if d%7:
        for m in range(d,NMAX+1,d): b[m]-=d
X=fmpq_poly([0,1]); c=[fmpq_poly([1])]
for n in range(1,NMAX+1):
    acc=fmpq_poly([0])
    for k in range(1,n+1):
        acc+=b[k]*c[n-k]
    c.append(acc*X/n)
    if n%100==0: print(n,'%.0fs'%(time.time()-t0),flush=True)
p=c[NMAX]; pickle.dump([[ (int(x.p),int(x.q)) for x in cc.coeffs()] for cc in c],open('call.pkl','wb'))
lo=fmpq(487350758538673426330055236816,10**29); hi=fmpq(487350758538673426343362426758,10**29)
print('sign at lo,hi:', p(lo)>0, p(hi)>0)
# exact root count in (lo,hi) via Descartes after Moebius map
p1=p(fmpq_poly([lo,hi-lo])); cf=p1.coeffs(); q=fmpq_poly(cf[::-1])(fmpq_poly([1,1]))
s=[x for x in q.coeffs() if x!=0]; V=sum(1 for i in range(len(s)-1) if (s[i]>0)!=(s[i+1]>0))
print('Descartes sign variations on (lo,hi):',V,' (1 => exactly one root)  %.0fs'%(time.time()-t0))
