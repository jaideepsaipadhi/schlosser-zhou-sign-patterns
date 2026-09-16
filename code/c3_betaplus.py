import flint, time, sys
from flint import fmpq_poly, fmpq, fmpz_poly
NMAX=int(sys.argv[1])
b=[0]*(NMAX+1)
for d in range(1,NMAX+1):
    r=d%8; s=-d if r in (1,7) else (d if r in (3,5) else 0)
    if s:
        for m in range(d,NMAX+1,d): b[m]+=s
X=fmpq_poly([0,1])
t=time.time()
c=[fmpq_poly([1])]
for n in range(1,NMAX+1):
    acc=fmpq_poly([0])
    for k in range(1,n+1):
        if b[k]: acc+=b[k]*c[n-k]
    c.append(acc*X/n)
print('polys built %.1fs'%(time.time()-t),flush=True)
A=fmpq(2664479110226973,10**15); B=fmpq(8,3)
sigma='+-++-+--'
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
bad=[]; zeros=[]
for n in range(NMAX+1):
    p=c[n]; want=1 if sigma[n%8]=='+' else -1
    ea,eb=p(A),p(B)
    issue=[]
    if sgn(ea)!=want: issue.append(('at 8/3',sgn(ea)))
    if sgn(eb)!=want:
        if eb==0: zeros.append(n)
        else: issue.append(('at 4',sgn(eb)))
    r=roots_in(p,A,B)
    if r: issue.append(('interior',r))
    if issue: bad.append((n,issue))
print('checked n=0..%d in %.1fs'%(NMAX,time.time()-t))
print('exact zeros at delta=4:',zeros)
print('problems:',bad)
