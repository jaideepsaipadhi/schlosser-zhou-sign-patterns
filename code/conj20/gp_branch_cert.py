# Certified check of the transformation formula for G_p (p=7,11) at all cusps h/k with p|k<=KMAX:
# log G_p(tau) = (1/2)log(p/p) + pi i(-s(h,k)+s(h,k/p)) + pi(p^2-p)/(12pz) - (p-1)pi z/(12k^2) + log Ghat,
# log Ghat = log f(e(p h'_p/k) e^{-2pi p/z}) - log f(e(sigma*h'/k) e^{-2pi/z}),  with log f(x) = -sum log(1-x^m)
# for sigma = -1 (h h' = -1 mod k convention) vs sigma = +1; the difference must be 0 (branch m=0).
import sys, math
from fractions import Fraction as F
from flint import arb, acb, ctx
ctx.prec=200
PI=arb.pi(); I=acb(0,1)
def ded(h,k):
    t=F(0)
    for j in range(1,k):
        x=F(j,k); y=F(j*h,k)
        fx=x-math.floor(x)-F(1,2) if x.denominator!=1 else 0
        fy=y-math.floor(y)-F(1,2) if y.denominator!=1 else 0
        t+=fx*fy
    return t
def logf(x,M=400):     # -sum_{m} log(1-x^m), |x|<=e^{-2pi/p}
    s=acb(0)
    for m in range(1,M): s-= (1-x**m).log()
    return s
def logG_series(p,tau,B=6000):
    s=acb(0)
    for m in range(1,B):
        # log G_p = sum_{p not| m} log(1-q^m) = -sum_n sigma'(n) q^n / n
        pass
    q=(2*PI*I*tau).exp(); s=acb(0)
    for m in range(1,B):
        if m%p: s+=(1-(2*PI*I*tau*m).exp()).log()
    r=abs(q).upper(); tail=2*arb(r)**B/(1-arb(r))            # |sum_{m>=B} log(1-q^m)| <= 2|q|^B/(1-|q|)
    return s+acb(arb(0,tail.upper()),arb(0,tail.upper()))
p=int(sys.argv[1]); KS=[int(x) for x in sys.argv[2].split(',')]
for sigma in (-1,1):
    worst=0; bad=[]
    if sigma>0: KS_=KS[:1]
    else: KS_=KS
    for k in KS_:
        for h in range(k):
            if math.gcd(h,k)!=1: continue
            z=arb(1); tau=acb(arb(h)/k,z/k**2)
            hp=pow(h,-1,k)*sigma % k
            hpp=pow(h,-1,k//p)                       # (hp/d) h'_d = 1 mod k/d with d=p
            x=acb.exp_pi_i(acb(arb(2*p*hpp)/k))*(-2*PI*p/z).exp()
            qt=acb.exp_pi_i(acb(arb(2*hp)/k))*(-2*PI/z).exp()
            phase=PI*I*(-(lambda v:arb(v.numerator)/v.denominator)(ded(h,k))+(lambda v:arb(v.numerator)/v.denominator)(ded(h,k//p)))
            rhs=phase+PI*(p*p-p)/(12*p*z)-(p-1)*PI*z/(12*k**2)+logf(x)-logf(qt)
            lhs=logG_series(p,tau)
            d=(lhs-rhs)/(2*PI*I)
            err=float(abs(d).upper()); worst=max(worst,err)
            if not abs(d)<arb('0.5'): bad.append((h,k,d.str(4)))
    print('p=%d sigma=%+d: max |(LHS-RHS)/2pi i| <= %.2e ; |.|<1/2 at all cusps (=> branch integer 0): %s'%(p,sigma,worst,len(bad)==0))
