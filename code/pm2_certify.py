# Certified large-n signs for delta = +2 and delta = -2 (self-contained constants).
import sys
from flint import arb, ctx
ctx.prec=100
PI=arb.pi(); E1=arb(1).exp(); R2=arb(2).sqrt(); K=arb('23.6475')
I1=lambda x:x.bessel_i(1)
t=(-2*PI).exp(); U=arb(0); V=arb(0)
for J in range(3,400):
    if J%16 in (2,14) and J>2: U+=t**(arb(J*J-4)/64)
    if J%16 in (6,10) and J>6: V+=t**(arb(J*J-36)/64)
F2m=(1-U)**(-2)*(1-V)**(-2)-1                 # F~_2(e^{-2pi}) - 1 (t-term kept: correction not growing at exponent 2)
def run(sign,NLO,NHI):
    s=arb(2); fails=[]; sig={}
    for n in range(NLO,NHI+1):
        y=(2*n+sign*s).sqrt()                  # +: sqrt(2n+delta), -: sqrt(2n-delta')
        P=lambda k: (2*PI/k)*s.sqrt()/y*I1((2*PI/k)*s.sqrt()*y)
        N=int(float((4*PI*n+2*PI*2).sqrt().upper()))+2
        G=[(8,2),(16,4)]+[(k,k//2) for k in range(24,N+1) if k%4==0]
        Sg=sum((arb(c)/k for k,c in G),arb(0))/(N+1)
        E=R2*PI*Sg*E1*(2*PI).exp()+2*R2*Sg*E1*(2*PI).exp()*F2m+2*R2*E1*K**2
        T24=sum((c*P(k) for k,c in G if k>=24),arb(0))
        if sign>0:
            a8=2*arb.cos_pi(arb(3*n)/4); a16=4*arb.cos_pi(arb(n)/8)*arb.cos_pi(arb(n+2)/2)
        else:
            a8=2*arb.cos_pi(arb(n-2)/4); a16=2*arb.cos_pi(arb(n-2)/8)+2*arb.cos_pi(arb(7*n+2)/8)
        if abs(a8)>0.1:
            m=abs(a8)*P(8)-4*P(16)-T24-E; lead=a8
        else:
            m=abs(a16)*P(16)-T24-E; lead=a16
        if not (m>0): fails.append(n)
        else: sig[n%16]=('+' if lead>0 else '-')
    return fails,sig
for sign in (1,-1):
    f,sg=run(sign,int(sys.argv[1]),int(sys.argv[2]))
    print('delta=%+d: failures %d, last %s; dominant-term sign pattern mod 16: %s'%(2*sign,len(f),f[-3:],''.join(sg.get(r,'?') for r in range(16))))
