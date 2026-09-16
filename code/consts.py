# Self-contained error constant (no He-Li input)
from flint import arb, ctx
PI=arb.pi(); E1=arb(1).exp(); R2=arb(2).sqrt()
KNG=arb('23.6475')                      # Step 2 (certified upper bound 23.64743...)
def Ftilde_minus(d):
    """F~_d(t)-1-d t at t=e^{-2pi}; F~=(1-U)^{-d}(1-V)^{-d}."""
    t=(-2*PI).exp()
    U=arb(0); V=arb(0)
    for J in range(3,400):
        if J%16 in (2,14) and J>2: U+=t**(arb(J*J-4)/64)
        if J%16 in (6,10) and J>6: V+=t**(arb(J*J-36)/64)
    tail=2*t**(arb(399*399-36)/64)     # remaining terms, super-geometric
    U+=arb(0,tail.upper()); V+=arb(0,tail.upper())
    return (1-U)**(-d)*(1-V)**(-d)-1-d*t
def Econst(d,N):
    """E(delta) for the Farey order N (any N >= sqrt(4 pi n + 2 pi delta)); increasing in d, decreasing in N."""
    t=(-2*PI).exp()
    # chord length <= 2*sqrt2*k/(N+1); small arcs <= (pi/2)*sqrt2*k/(N+1); sum_{growing} 1/k <= 1/2 + N/8
    S=(arb(1)/2+arb(N)/8)/(N+1)
    a=2*(PI/2)*R2*S*E1*(PI*d).exp()*(1+d*t)
    b=2*R2*S*E1*(PI*d).exp()*Ftilde_minus(d)
    c=2*R2*E1*KNG**d
    return a+b+c
if __name__=='__main__':
    ctx.prec=100
    for d in ('2.6667','3','4'):
        print(d, Econst(arb(d),81).str(5))
