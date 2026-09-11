# Q12^delta, -1<delta<0: exact coefficient polynomials; find, for n<=NMAX, the roots in (-1,0) and the
# resulting maximal interval (-1, d1] on which the conjectured pattern +++++------+ holds for n<=NMAX.
import sys
from flint import fmpq_poly, fmpq, fmpz_poly
NMAX=int(sys.argv[1])
eps=lambda d:(1 if d%12 in (1,11) else (-1 if d%12 in (5,7) else 0))
b=[0]*(NMAX+1)
for d in range(1,NMAX+1):
    e=eps(d)
    if e:
        for m in range(d,NMAX+1,d): b[m]-=e*d
X=fmpq_poly([0,1]); c=[fmpq_poly([1])]
for n in range(1,NMAX+1):
    acc=fmpq_poly([0])
    for k in range(1,n+1):
        if b[k]: acc+=b[k]*c[n-k]
    c.append(acc*X/n)
pat='+++++------+'
print('c_{-1/4}(6) =', c[6](fmpq(-1,4)), ' (pattern requires sign -)')
best=None
for n in range(1,NMAX+1):
    p=c[n]; want=1 if pat[n%12]=='+' else -1
    num=fmpz_poly([int(x*p.denom()) for x in p.coeffs()])
    rts=[r for r,m in num.complex_roots() if abs(r.imag)<1e-30 and -1<float(r.real.mid())<0]
    for r in rts:
        x=float(r.real.mid())
        # sign just below the root (towards -1) should be 'want'; the pattern fails just above x if sign flips
        if best is None or x<best[0]: best=(x,n)
print('smallest root in (-1,0) among c_n, n<=%d: delta=%.15f at n=%d'%(NMAX,best[0],best[1]))
