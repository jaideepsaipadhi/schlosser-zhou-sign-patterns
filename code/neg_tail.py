# Tail n>20000 for -1<delta<=-0.99: every competitor/main ratio is C*y^p*e^{-g y} (g>0), decreasing for y>=200;
# verify the sum of ratios at y=200 (n=20000) is < 1, using envelopes of Lemma 5.2.
from flint import arb, ctx
ctx.prec=100
exec(open('neg_certify.py').read().split('def P(k,s,n)')[0].split("I1=lambda")[0])
PI=arb.pi(); E1=arb(1).exp(); R2=arb(2).sqrt(); KP=arb('23.6475')
t=(-2*PI).exp(); U=arb(0); V=arb(0)
for J in range(3,400):
    if J%16 in (2,14) and J>2: U+=t**(arb(J*J-4)/64)
    if J%16 in (6,10) and J>6: V+=t**(arb(J*J-36)/64)
F1=(1-U)**(-1)*(1-V)**(-1); F1m=F1-1; ell=-(1-U).log()-(1-V).log()
y=arb(200)                     # y = sqrt(2n - s) >= 199.99 for n>20000; use 199.9 for lower bounds
yl=arb('199.9')
N=(4*PI).sqrt()*y/(2).__float__()**0.5+3   # N <= sqrt(4 pi n)+2 <= sqrt(2 pi) y + 3
N=(2*PI).sqrt()*y+3
a8=(PI/4)*arb('0.99').sqrt(); x0=a8*yl
kap=(2*PI*x0).sqrt()*(-x0).exp()*x0.bessel_i(1)
Uenv=lambda x: x.exp()/(2*PI*x).sqrt()
Main=(PI/2)*arb('0.99999')*(PI/4)*arb('0.99').sqrt()/y*kap*Uenv(a8*yl)      # lower bound, degenerate case
Main2=2*arb('0.7016')*(PI/4)*arb('0.99').sqrt()/y*kap*Uenv(a8*yl)           # nondegenerate
Pup=lambda k: (2*PI/k)/yl*Uenv((2*PI/k)*y)
a24=2*PI/24
D24=arb('0.506')+1/(2*yl*yl)+arb('0.5026')*a24*y*(1+1/(a24*arb('0.99').sqrt()*yl))
Sg=arb('0.5')+N/8
k16=(PI/2)*Pup(16)
k24=(N/4)*(N/2)*(N*N/2+arb('0.01')+D24)*Pup(24)
arcs=Sg*2*E1*PI.exp()*((N*N/2+1)*(PI/R2)/1+(2+PI/2)/arb('0.99')+1)
chord=(2*R2*Sg/(N+1))*E1*PI.exp()*((N*N/2+2*PI*N*N+1)*F1m+F1*ell)*N
ng=2*R2*E1*KP*(KP.log()+1+arb('0.556')*N*N)
tot=(k16+k24+arcs+chord+ng)/Main
Epr=(R2*PI*Sg/(N+1)+2*R2*Sg/(N+1)*F1m)*E1*PI.exp()+2*R2*E1*KP
tot2=(4*Pup(16)+(N/4)*(N/2)*Pup(24)+Epr)/Main2
print('degenerate residues: sum of ratios at y=200:',tot.str(4),' <1:',tot<1)
print('nondegenerate residues: sum of ratios at y=200:',tot2.str(4),' <1:',tot2<1)
print('exponential gaps: vs k=16:',(a8-PI/8).str(3),' vs k=24:',(a8-a24).str(3),' vs O(poly):',a8.str(3))
