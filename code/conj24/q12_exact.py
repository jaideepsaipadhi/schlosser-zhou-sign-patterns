# Exact Gauss-sum ratios (mod Phi_{24k}) and certified branch integers for Q12 at the growing cusps
# for delta>0 (h=5,7 mod 12) and delta<0 (h=1,11 mod 12), k=12,24.
import math
from flint import fmpz_poly, arb, acb, ctx
ctx.prec=160
PI=arb.pi(); I=acb(0,1)
def psi(a,j):
    r=a-6
    if (j-r)%12: return 0
    return -1 if ((j-r)//12)%2 else 1
def Gvec(a,h,k,l):
    M=24*k; c=[0]*M
    for s in range(M):
        p=psi(a,s)
        if p:
            c[(h*s*s+l*s)%M]+=p
            if l: c[(h*s*s-l*s)%M]+=p
    return c
def Gp(a,h,k,l): return fmpz_poly(Gvec(a,h,k,l))
def Gacb(a,h,k,l):
    M=24*k; s=acb(0)
    for i,x in enumerate(Gvec(a,h,k,l)):
        if x: s+=x*acb.exp_pi_i(acb(arb(2*i)/M))
    return s
def root_ratio(pn,pd,k):          # find t with pn == zeta^t pd mod Phi
    M=24*k; Phi=fmpz_poly.cyclotomic(M)
    for t in range(M):
        mono=fmpz_poly([0]*t+[1])
        if (pn-mono*pd)%Phi==0: return t
    return None
Bm=8000
b=[0]*(Bm+1)
for d in range(1,Bm+1):
    e=1 if d%12 in (1,11) else (-1 if d%12 in (5,7) else 0)
    if e:
        for m in range(d,Bm+1,d): b[m]-=e*d
def logQ(tau):
    s=acb(0)
    for m in range(1,Bm+1):
        if b[m]: s+=(2*PI*I*tau*m).exp()*b[m]/m
    rr=float(abs((2*PI*I*tau).exp()).upper()); LB=1+math.log(Bm)
    tail=arb(rr)**Bm*(LB*arb(rr)/(1-arb(rr))+arb(rr)/(Bm*(1-arb(rr))**2))
    return s+acb(arb(0,tail.upper()),arb(0,tail.upper()))
def analyse(h,k,sign):
    na,da=(1,5) if sign>0 else (5,1)     # R^{sign}: numerator lead J=2, denominator lead J=10
    M=24*k
    t=root_ratio(Gp(na,h,k,2),Gp(da,h,k,10),k)          # c = zeta^t
    t1=root_ratio(Gp(da,h,k,14),Gp(da,h,k,10),k)        # rho1 = zeta^t1 (first denominator correction)
    # branch: sign*log Q12 = sign*(-2 pi i tau) + Log c + 96*pi/(48 z) + Log(1+u) - Log(1+v) + 2 pi i m
    z=arb(1); tau=acb(arb(h)/k,z/k**2); X=(-PI/(48*z)).exp()
    Gn={l:Gacb(na,h,k,l) for l in range(2,70)}; Gd={l:Gacb(da,h,k,l) for l in range(10,70)}
    u=sum((Gn[l]/Gn[2]*X**(l*l-4) for l in range(3,70)),acb(0)); v=sum((Gd[l]/Gd[10]*X**(l*l-100) for l in range(11,70)),acb(0))
    tr=4*M*X**(70*70-100); u+=acb(arb(0,tr.upper()),arb(0,tr.upper())); v+=acb(arb(0,tr.upper()),arb(0,tr.upper()))
    c=acb.exp_pi_i(acb(arb(2*t)/M))
    L=sign*(-2*PI*I*tau)+c.log()+2*PI/z+(1+u).log()-(1+v).log()
    d=(sign*logQ(tau)-L)/(2*PI*I); m=round(float(d.real.mid())); ok=abs(d-m)<arb('0.25')
    phase=(c.arg()+2*PI*m-sign*2*PI*arb(h)/k)/PI       # D-phase/pi  (D=|delta|):  D*(Arg c + 2 pi m - sign*2 pi h/k)
    return t,t1,m,ok,phase
for sign,hs in ((1,[(5,12),(7,12),(5,24),(7,24),(17,24),(19,24)]),(-1,[(1,12),(11,12),(1,24),(11,24),(13,24),(23,24)])):
    for h,k in hs:
        t,t1,m,ok,ph=analyse(h,k,sign)
        print('sign %+d  h/k=%2d/%d  c=zeta_%d^%d  rho1=zeta^%s  branch m=%d certified=%s  |delta|-phase/pi=%s'%(sign,h,k,24*k,t,t1,m,ok,ph.str(8)))
