# Negative delta, -1 < delta <= -0.99 (delta' = -delta in [0.99,1)): certified large-n criterion.
# Nondegenerate residues: standard main-term domination, uniform in delta'.
# Degenerate residues n = 3,7 mod 8: epsilon-free difference-quotient criterion (see neg_step2.md).
import sys, time
from flint import arb, ctx
ctx.prec=100
PI=arb.pi(); E1=arb(1).exp(); R2=arb(2).sqrt()
KP=arb('23.6475')
I1=lambda x: x.bessel_i(1)
t=(-2*PI).exp()
U=arb(0); V=arb(0)
for J in range(3,400):
    if J%16 in (2,14) and J>2: U+=t**(arb(J*J-4)/64)
    if J%16 in (6,10) and J>6: V+=t**(arb(J*J-36)/64)
F1=(1-U)**(-1)*(1-V)**(-1); F1m=F1-1; ell=-(1-U).log()-(1-V).log()
def P(k,s,n):
    y=(2*n-s).sqrt(); a=2*PI/k
    return a*s.sqrt()/y*I1(a*s.sqrt()*y)
def growing(N):          # (k, count) with count >= number of inverse-growing h
    out=[(8,2),(16,4)]
    out+=[(k,k//2) for k in range(24,N+1) if k%4==0]
    return out
S_LO=sys.argv[3] if len(sys.argv)>3 else '0.99'
s_lo=arb(S_LO); s_hi=arb(1)
epsmax=1-s_lo
x=PI*epsmax/4; SINC=1-x*x/6          # sin(pi eps/4)/(pi eps/4) >= 1 - x^2/6
CMIN={'0.99':arb('0.7016'),'0.772':arb('0.5698')}[S_LO]
NLO,NHI=int(sys.argv[1]),int(sys.argv[2])
fails=[]; t0=time.time()
for n in range(NLO,NHI+1):
    N=int(float((4*PI*n).sqrt().upper()))+2
    G=growing(N)
    Sg=sum((arb(c)/k for k,c in G),arb(0))           # sum_{growing} 1/k
    Epr=(R2*PI*Sg/(N+1))*E1*PI.exp()+(2*R2*Sg/(N+1))*E1*PI.exp()*F1m+2*R2*E1*KP
    r=n%8
    if r in (3,7):
        y=(2*n-s_hi).sqrt()
        main=(PI/2)*SINC*P(8,s_lo,n)
        rest=(PI/2)*P(16,s_hi,n)
        for k,c in G:
            if k<24: continue
            a=2*PI/k
            D=1/(2*s_lo)+1/(2*y*y)+a*y/(2*s_lo.sqrt())*(1+1/(a*s_lo.sqrt()*y))
            rest+=c*((arb(k*k)/2+arb('0.01'))+D)*P(k,s_hi,n)
        # derivative bounds of the error (independent of epsilon)
        arcs=sum((arb(c)/(k*k)*2*E1*PI.exp()*((arb(k*k)/2+arb('0.01'))*(PI/R2)*k/(N+1)+(2+PI/2)/s_lo+(PI*PI/R2)/(k*(N+1))) for k,c in G),arb(0))
        chord=sum((arb(c)/(k*k)*(2*R2*k/(N+1))*E1*PI.exp()*((arb(k*k)/2+2*PI*N*N/arb(k*k)+arb('0.1'))*F1m+F1*ell) for k,c in G),arb(0))
        ng=2*R2*E1*KP*(KP.log()+arb('0.01')+arb('0.556')*N*N)
        m=main-rest-arcs-chord-ng
    else:
        main=2*CMIN*P(8,s_lo,n)
        T=4*P(16,s_hi,n)+sum((c*P(k,s_hi,n) for k,c in G if k>=24),arb(0))
        m=main-T-Epr
    if not (m>0): fails.append(n)
print('n=%d..%d: %.1fs, failures: %d'%(NLO,NHI,time.time()-t0,len(fails)), fails[-20:])
