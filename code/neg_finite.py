# Exact finite check: for n<=NMAX and all delta in [-1,-0.99], c_delta(n) has sign sigma'(n mod 8),
# sigma' = (+,+,+,-,-,-,-,+); where c_{-1}(n)=0 we divide by (delta+1) and check the quotient on the closed interval.
import sys,time
from flint import fmpq_poly, fmpq
exec(open('c3.py').read().split("A=fmpq(8,3)")[0].replace("NMAX=int(sys.argv[1])","NMAX=int(sys.argv[1])"))
def sgn(v): return (v>0)-(v<0)
def descartes(p):  # sign variations of p on (0,1) via t=1/(1+s)... roots in (0,1)
    cf=[x for x in p.coeffs()]
    d=len(cf)-1
    rev=fmpq_poly(cf[::-1])
    q=rev(fmpq_poly([1,1]))
    s=[sgn(x) for x in q.coeffs() if x!=0]
    return sum(1 for i in range(len(s)-1) if s[i]!=s[i+1])
def roots_in(p,a,bb,depth=0):
    # number of distinct-sign-change certification on open (a,bb); returns list of subintervals possibly containing roots
    p1=p(fmpq_poly([a,bb-a]))
    v=descartes(p1)
    if v==0: return []
    if v==1: return [(a,bb,'exactly one root')]
    if depth>40: return [(a,bb,'unresolved')]
    m=(a+bb)/2
    out=[]
    if p(m)==0: out.append((m,m,'root at midpoint'))
    return out+roots_in(p,a,m,depth+1)+roots_in(p,m,bb,depth+1)

A=fmpq(-1); B=fmpq(-99,100)
sig='+++----+'
def sgn(v): return (v>0)-(v<0)
bad=[]; divided=[]
for n in range(NMAX+1):
    p=c[n]; want=1 if sig[n%8]=='+' else -1
    if p(A)==0:
        q,r=divmod(p,fmpq_poly([1,1])); assert r==0; p=q; divided.append(n)
        assert p(A)!=0
    issue=[]
    for e in (A,B):
        if sgn(p(e))!=want: issue.append(('end',float(e),sgn(p(e))))
    rr=roots_in(p,A,B)
    if rr: issue.append(('interior',rr))
    if issue: bad.append((n,issue))
print('n<=%d on [-1,-0.99]: divided by (delta+1) at %d indices (residues %s); problems: %s'%(NMAX,len(divided),sorted(set(d%8 for d in divided)),bad[:10]))
