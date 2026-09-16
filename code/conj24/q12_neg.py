# Conj 24, -1<delta<=delta1 (s=-delta in [S_LO,1)): certified large-n criterion. Degenerate residues n=5,11 mod 12
# via the epsilon-free difference-quotient bound; others via main-term domination on s-subintervals.
import sys, time
from flint import arb, ctx
ctx.prec=100
PI=arb.pi(); E1=arb(1).exp(); R2=arb(2).sqrt()
KP=arb('547.23'); LOGK=arb('550.75').log()
I1=lambda x:x.bessel_i(1)
t=(-2*PI).exp(); U=arb(0); V=arb(0)
for J in range(3,600):
    if J%24 in (2,22) and J>2: U+=t**(arb(J*J-4)/96)
    if J%24 in (10,14) and J>10: V+=t**(arb(J*J-100)/96)
F1=(1-U)**(-1)*(1-V)**(-1); F1m=F1-1; ell=-(1-U).log()-(1-V).log()
S_LO=arb(sys.argv[3]); NLO,NHI=int(sys.argv[1]),int(sys.argv[2]); NSUB=int(sys.argv[4]) if len(sys.argv)>4 else 16
def P(k,s,n):          # (2pi/k) sqrt(B/C) I1((2pi/k) sqrt(BC)), B=2s, C=2(n-s)
    B=2*s; C=2*(n-s); a=2*PI/k
    return a*(B/C).sqrt()*I1(a*(B*C).sqrt())
def growing(N): return [(12,2),(24,4),(36,6)]+[(k,k//2) for k in range(40,N+1,4)]
epsmax=1-S_LO; x=PI*epsmax/3; SINC=1-x*x/6
pts=[S_LO+(1-S_LO)*i/NSUB for i in range(NSUB+1)]
fails=[]; t0=time.time()
for n in range(NLO,NHI+1):
    N=int(float((4*PI*n).sqrt().upper()))+2
    G=growing(N); Sg=sum((arb(c)/k for k,c in G),arb(0))
    r=n%12
    if r in (5,11):
        s1=arb(1); y=(2*(n-s1)).sqrt()
        main=(2*PI/3)*SINC*P(12,S_LO,n)
        rest=(2*PI/3)*P(24,s1,n)
        cnt40=sum(k//2 for k in range(40,N+1,4))
        for k,c in ((36,6),(40,cnt40)):
            a=2*PI/k; xx=a*(2*S_LO*2*(n-1)).sqrt()
            D=1/(2*S_LO)+1/(2*(n-s1))+(1+1/xx)*a*(n/S_LO).sqrt()
            kk=k if k==36 else N                                # |theta|<=k^2/2: use N for the lumped k>=40 part
            rest+=c*((arb(kk)**2/2+arb('0.01'))+D)*P(k,s1,n)
        arcs=sum((arb(c)/(k*k)*2*E1*(2*PI).exp()*((arb(k*k)/2+arb('0.01'))*(PI/R2)*k/(N+1)+(2+PI/2)/S_LO+(2*PI*PI/R2)/(k*(N+1))) for k,c in G),arb(0))
        chord=sum((arb(c)/(k*k)*(2*R2*k/(N+1))*E1*(2*PI).exp()*((arb(k*k)/2+4*PI*N*N/arb(k*k)+arb('0.2'))*F1m+F1*ell) for k,c in G),arb(0))
        ng=2*R2*E1*KP*(LOGK+arb('0.02')+arb('0.556')*N*N)
        m=main-rest-arcs-chord-ng
        if not (m>0): fails.append((n,'deg'))
    else:
        E=R2*PI*(Sg/(N+1))*E1*(2*PI).exp()+2*R2*(Sg/(N+1))*E1*(2*PI).exp()*F1m+2*R2*E1*KP
        cnt40=sum(k//2 for k in range(40,N+1,4))
        T=4*P(24,arb(1),n)+6*P(36,arb(1),n)+cnt40*P(40,arb(1),n)
        for i in range(NSUB):
            sa,sb=pts[i],pts[i+1]
            ca=abs(arb.cos_pi((n-2*sa)/6)); cb=abs(arb.cos_pi((n-2*sb)/6))
            cmin=ca if ca<cb else cb
            m=2*cmin*P(12,sa,n)-T-E
            if not (m>0): fails.append((n,float(sa))); break
print('n=%d..%d (s in [%s,1)): %.1fs, failures %d, last %s'%(NLO,NHI,S_LO.str(6),time.time()-t0,len(fails),fails[-4:]))
