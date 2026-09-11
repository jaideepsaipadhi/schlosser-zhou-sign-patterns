# Exact finite check for G11^delta on [gamma,2] (n<=NMAX), gamma = root of c_21 in (1.5,2).
import sys, time
from flint import fmpq_poly, fmpq, fmpz_poly
NMAX=int(sys.argv[1]); pat='+--+++-+--+'
b=[0]*(NMAX+1)
for d in range(1,NMAX+1):
    if d%11:
        for m in range(d,NMAX+1,d): b[m]-=d
X=fmpq_poly([0,1]); c=[fmpq_poly([1])]
for n in range(1,NMAX+1):
    acc=fmpq_poly([0])
    for k in range(1,n+1): acc+=b[k]*c[n-k]
    c.append(acc*X/n)
def sgn(v): return (v>0)-(v<0)
def desc(p):
    q=fmpq_poly(p.coeffs()[::-1])(fmpq_poly([1,1])); s=[sgn(x) for x in q.coeffs() if x!=0]
    return sum(1 for i in range(len(s)-1) if s[i]!=s[i+1])
def roots_in(p,a,bb,dep=0):
    v=desc(p(fmpq_poly([a,bb-a])))
    if v==0: return []
    if v==1 or dep>60: return [(a,bb)]
    m=(a+bb)/2; return ([(m,m)] if p(m)==0 else [])+roots_in(p,a,m,dep+1)+roots_in(p,m,bb,dep+1)
# gamma: isolate the root of c_21 in (1.5,2)
num=fmpz_poly([int(x*c[21].denom()) for x in c[21].coeffs()])
fac=[f for f,e in num.factor()[1] if any(1.5<float(z.real.mid())<2 and abs(z.imag)<1e-30 for z,mm in f.complex_roots())]
P18=fac[0]; print('factor containing gamma: degree',P18.degree(), ' roots in (1.5,2):',[z.real.str(20) for z,mm in P18.complex_roots() if abs(z.imag)<1e-30 and 1.5<float(z.real.mid())<2])
G_lo=fmpq(17584535519419,10**13); G_hi=fmpq(17584535519420,10**13)
P18q=fmpq_poly([fmpq(int(x)) for x in P18.coeffs()])
print('P18 sign change in (G_lo,G_hi):', sgn(P18q(G_lo))!=sgn(P18q(G_hi)), ' roots of P18 in [G_hi,2]:', roots_in(P18q,G_hi,fmpq(2)))
A=G_lo; B=fmpq(2); bad=[]; notes=[]
for n in range(NMAX+1):
    p=c[n]; want=1 if pat[n%11]=='+' else -1; f=1
    if n==21:
        q,r=divmod(p,P18q); assert r==0; p=q; f=sgn(P18q(fmpq(19,10))); notes.append('c21/P18')
    while p(B)==0 and not p.is_zero():
        q,r=divmod(p,fmpq_poly([-B,1])); p=q; f=-f; notes.append((n,2))
    if p.is_zero(): continue
    iss=[]
    for e in (A,B):
        if sgn(p(e))*f!=want: iss.append(('end',float(e)))
    rr=roots_in(p,A,B)
    if rr: iss.append(('roots',len(rr)))
    if iss: bad.append((n,iss))
print('n<=%d on [gamma_lo,2]: notes %s; problems %s'%(NMAX,notes[:8],bad[:6]))
