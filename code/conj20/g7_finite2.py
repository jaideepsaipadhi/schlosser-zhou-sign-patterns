# Exact check, from pickled exact polynomials: for N1<=n<=N2, sign of c_n(delta) on the closed interval [A,B]
# equals pattern (dividing by (delta-A) at endpoint zeros, factor >=0 on [A,B]).
import sys, pickle, time
from flint import fmpq_poly, fmpq
src=open('/home/claude/c20/c3.py').read() if False else None
def sgn(v): return (v>0)-(v<0)
def descartes(p):
    cf=p.coeffs(); q=fmpq_poly(cf[::-1])(fmpq_poly([1,1]))
    s=[sgn(x) for x in q.coeffs() if x!=0]
    return sum(1 for i in range(len(s)-1) if s[i]!=s[i+1])
def roots_in(p,a,b,depth=0):
    v=descartes(p(fmpq_poly([a,b-a])))
    if v==0: return []
    if v==1: return [(a,b)]
    if depth>60: return [(a,b,'unres')]
    m=(a+b)/2
    return ([(m,m)] if p(m)==0 else [])+roots_in(p,a,m,depth+1)+roots_in(p,m,b,depth+1)
N1,N2=int(sys.argv[1]),int(sys.argv[2])
A=fmpq(*map(int,sys.argv[3].split('/'))); B=fmpq(*map(int,sys.argv[4].split('/'))); pat=sys.argv[5]
raw=pickle.load(open('call.pkl','rb'))
t0=time.time(); bad=[]; div=[]
for n in range(N1,N2+1):
    p=fmpq_poly([fmpq(a,b) for a,b in raw[n]])
    want=1 if pat[n%len(pat)]=='+' else -1
    fac=1
    while p(A)==0 and not p.is_zero():
        q,r=divmod(p,fmpq_poly([-A,1])); p=q; div.append(n)
    while p(B)==0 and not p.is_zero():
        q,r=divmod(p,fmpq_poly([-B,1])); p=q; div.append(n); fac=-fac
    if p.is_zero(): bad.append((n,'zero poly')); continue
    iss=[]
    for e in (A,B):
        if sgn(p(e))!=want*fac: iss.append(('end',float(e)))
    rr=roots_in(p,A,B)
    if rr: iss.append(('roots',len(rr)))
    if iss: bad.append((n,iss))
print('n=%d..%d on [%s,%s]: divided at A for %d indices; problems %s; %.0fs'%(N1,N2,float(A),float(B),len(div),bad[:6],time.time()-t0))
