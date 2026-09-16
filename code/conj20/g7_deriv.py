# Positivity of d/d delta c_delta(n) for n = r mod 7 (default r=1), all n >= n1, delta in [a,b] (4 < a).
# Dominant: amp7'(delta) P7(delta/2); competitors bounded; ratios are C y^p e^{-g y} (decreasing in n).
import sys
from g7lib import *
def D(k,B,C):  # bound for |d/ds log P_k(B(s),C(s))|, B=s/2 (or s/2-2), C=2n-s/2
    x=(2*PI/k)*(B*C).sqrt()
    return 1/(4*B)+1/(4*C)+(1+1/x)*(2*PI/k)*C.sqrt()/(4*B.sqrt())
def check(a,b,n1,r):
    A=arb(a); Bb=arb(b); db=arb(((A+Bb)/2).mid(),float(((Bb-A)/2).mid())*1.000001+1e-45)
    n=n1+((r-n1)%7); N=int(float((4*PI*n).sqrt().upper()))+2
    Clo=2*n-Bb/2; Chi=2*n-A/2
    der=arb(0); amp=arb(0)
    for h,sh,hp in shifts(7):
        der+=arb.sin_pi(-db*fr(sh)-arb(2*h*n)/7)*PI*fr(sh)
        amp+=arb.cos_pi(-db*fr(sh)-arb(2*h*n)/7)
    Blo=A/2; Bhi=Bb/2
    P7lo=Plow(7,Blo,Clo,Chi)
    P7hi=(2*PI/7)*(Bhi/Clo).sqrt()*((2*PI/7)*(Bhi*Chi).sqrt()).bessel_i(1)
    # the |amp| D P7 term grows relative to the dominant term; it is only needed for n<2500 (dominance covers n>=2500):
    # bound its ratio by its value at n=2500.
    n2=2500; C2lo=2*n2-Bb/2; C2hi=2*n2-A/2
    rat=abs(amp)*D(7,Blo,C2lo)*P(7,Bhi,C2lo)*(C2hi/C2lo).sqrt()*((2*PI/7)*((Bhi*C2hi).sqrt()-(Bhi*C2lo).sqrt())).exp()/Plow(7,Blo,C2lo,C2hi)
    dom=arb(der.lower())*P7lo - rat*P7lo
    comp=arb(0)
    B2=Bhi-2
    # k=7 correction
    comp+=(6+6*Bb*PI*arb(5)/14+6*Bb*D(7,A/2-2,Clo))*P(7,B2,Clo)*(Chi/Clo).sqrt()
    for k,ph in ((14,6),(21,12)):
        comp+=ph*((PI*k/10+D(k,Blo,Clo))*P(k,Bhi,Clo)+(1+Bb*PI*k/10+Bb*D(k,A/2-2,Clo))*P(k,B2,Clo))*(Chi/Clo).sqrt()
    comp+=(arb(N)/7)*(arb(6)*N/7)*((PI*N/10+D(28,Blo,Clo))*P(28,Bhi,Clo)+(1+Bb*PI*N/10+Bb*D(28,A/2-2,Clo))*P(28,B2,Clo))*(Chi/Clo).sqrt()
    # error derivative
    t=(-2*PI).exp(); F=f_prod(t)**Bb*f_prod(t**7)**Bb; l7=(f_prod(t)*f_prod(t**7)).log()
    S=arb(6)/49
    arcs=2*S*E1*(PI*Bb/2).exp()*( (PI*N/10)*(PI/R2)*N + (2+PI/2)*(1/(A/2)+Bb*t/(A/2-2))*(1+(arb(N)/7).log()) + 10 )
    chords=2*R2*S*E1*(PI*Bb/2).exp()*((PI*N/10+PI*N*N+arb('0.1'))*(F-1-Bb*t)+F*l7-t)*N
    ng=2*R2*E1*K1**Bb*(arb('0.556')*N*N+arb('0.449')*N*N+2+K1.log())
    comp+=arcs+chords+ng
    return bool(dom>comp),(dom/comp).str(3)
if __name__=='__main__':
    print(check(sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]) if len(sys.argv)>4 else 1))
