# Common certified quantities for G7^delta = ((q;q)/(q^7;q^7))^delta (arb ball arithmetic).
PP=11
import math
from fractions import Fraction as F
from flint import arb, acb, ctx
ctx.prec=120
PI=arb.pi(); E1=arb(1).exp(); R2=arb(2).sqrt()
def dedekind(h,k):
    tot=F(0)
    for j in range(1,k):
        x=F(j,k); y=F(j*h,k)
        fx=x-math.floor(x)-F(1,2) if x.denominator!=1 else 0
        fy=y-math.floor(y)-F(1,2) if y.denominator!=1 else 0
        tot+=fx*fy
    return tot
SH={}
def shifts(k):          # list of (h, s(h,k)-s(h,k/7), h') for 7|k,  hh' = 1 mod k
    if k not in SH:
        SH[k]=[(h,dedekind(h,k)-dedekind(h,k//PP),pow(h,-1,k)) for h in range(k) if math.gcd(h,k)==1]
    return SH[k]
def fr(x): return arb(x.numerator)/x.denominator
def amp_main(k,n,d):    # sum_h Re( e^{-pi i d s_h} e(-hn/k) )
    s=arb(0)
    for h,sh,hp in shifts(k): s+=arb.cos_pi(-d*fr(sh)-arb(2*h*n)/k)
    return s
def amp_corr(k,n,d):    # sum_h Re( e^{-pi i d s_h} e(-hn/k) * (-d) e(-h'/k) )
    s=arb(0)
    for h,sh,hp in shifts(k): s+=arb.cos_pi(-d*fr(sh)-arb(2*h*n)/k-arb(2*hp)/k)
    return -d*s
def f_prod(x,T=400):
    v=arb(1)
    for m in range(1,T): v*=(1-x**m)
    return 1/v            # x<=1/2: tail negligible vs radius; add explicit tail
K1=arb(PP).sqrt()*(-PI*(PP-1)/(12*PP)).exp()*f_prod((-2*PI/PP).exp())*f_prod((-2*PI).exp())*(1+arb(10)**-40)
def F7m(d,extract):     # f(t)^d f(t^7)^d - 1 - [extract] d t  at t=e^{-2pi}
    t=(-2*PI).exp()
    return f_prod(t)**d*f_prod(t**PP)**d-1-(d*t if extract else 0)
def Econ(d,N):          # error constant, d = upper end of delta-interval (E increasing in d)
    ext = bool(d>arb(24)/(PP-1))
    t=(-2*PI).exp()
    S=arb(PP-1)/(PP*PP)
    a=R2*PI*S*E1*(PI*d*(PP-1)/12).exp()*(1+(d*t if ext else 0))
    b=2*R2*S*E1*(PI*d*(PP-1)/12).exp()*F7m(d,ext)
    c=2*R2*E1*K1**d
    return a+b+c
Uenv=lambda x:x.exp()/(2*PI*x).sqrt()
def kappa(x0): return (2*PI*x0).sqrt()*(-x0).exp()*x0.bessel_i(1)
def P(k,B,C):   # ENVELOPE upper bound for (2pi/k) sqrt(B/C) I1((2pi/k) sqrt(BC)); valid for all n (U(x)>=I1(x))
    if not (B>0): return arb(0)
    x=(2*PI/k)*(B*C).sqrt()
    return (2*PI/k)*(B/C).sqrt()*Uenv(x)
def Plow(k,B,Clo,Chi):   # lower bound valid for all n>=n1 (C>=Clo at n1): kappa(x1)*U(x)
    x=(2*PI/k)*(B*Clo).sqrt()
    return (2*PI/k)*(B/Chi).sqrt()*kappa(x)*Uenv(x)
def amp_main_mv(k,n,dball):
    # mean-value (centred) enclosure: amp(m) + amp'(J)*(J-m)
    m=arb(dball.mid()); r=arb(dball.rad())
    a0=amp_main(k,n,m)
    der=arb(0)
    for h,sh,hp in shifts(k): der+=-arb.sin_pi(-dball*fr(sh)-arb(2*h*n)/k)*(-PI*fr(sh))
    return a0+der*arb(0,r)
