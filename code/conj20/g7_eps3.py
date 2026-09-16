# epsilon-free criterion near delta=3 for residues r in {2,4,5}: c_3(n)=0 there; for 0<|delta-3|<=eta and n>=n1,
# sign c_delta(n) = sign(amp7(r,delta)) = sign((delta-3)*amp7'(r,3)).
import sys
from g7lib import *
from g7_deriv import D
def check(side,eta,n1,r):
    E3=arb(3); et=arb(eta)
    lo,hi=(E3,E3+et) if side>0 else (E3-et,E3)
    db=arb(((lo+hi)/2).mid(),float(((hi-lo)/2).mid())*1.000001+1e-45)
    n=n1+((r-n1)%7); N=int(float((4*PI*n).sqrt().upper()))+2
    Clo=2*n-hi/2; Chi=2*n-lo/2
    d1=arb(0); d2=arb(0)
    for h,sh,hp in shifts(7):
        d1+=arb.sin_pi(-E3*fr(sh)-arb(2*h*n)/7)*PI*fr(sh)                 # amp'(3)
        d2+=(PI*fr(sh))**2                                                 # |amp''| bound
    lower=abs(d1)-et*d2                                                    # |amp(delta)|/eps >= |amp'(3)| - eps*sup|amp''|
    Blo=lo/2; Bhi=hi/2
    dom=lower*Plow(7,Blo,Clo,Chi)
    comp=arb(0)
    for k,ph in ((14,6),(21,12)):
        comp+=ph*(PI*k/10+D(k,Blo,Clo))*P(k,Bhi,Clo)*(Chi/Clo).sqrt()
    comp+=(arb(N)/7)*(arb(6)*N/7)*(PI*N/10+D(28,Blo,Clo))*P(28,Bhi,Clo)*(Chi/Clo).sqrt()
    t=(-2*PI).exp(); F=f_prod(t)**hi*f_prod(t**7)**hi; l7=(f_prod(t)*f_prod(t**7)).log()
    S=arb(6)/49
    arcs=2*S*E1*(PI*hi/2).exp()*((PI*N/10)*(PI/R2)*N+(2+PI/2)/(lo/2)*(1+(arb(N)/7).log())+10)
    chords=2*R2*S*E1*(PI*hi/2).exp()*((PI*N/10+PI*N*N+arb('0.1'))*(F-1)+F*l7)*N
    ng=2*R2*E1*K1**hi*(arb('0.556')*N*N+arb('0.449')*N*N+2+K1.log())
    comp+=arcs+chords+ng
    sgn=(1 if d1>0 else -1)*side
    return bool(dom>comp),sgn,(dom/comp).str(3)
if __name__=='__main__':
    for side,pat in ((1,'+-++---'),(-1,'+--+++-')):
        for r in (2,4,5):
            ok,sg,mg=check(side,'1e-9',898,r)
            print('side',side,'r',r,'criterion',ok,'margin',mg,'sign',sg,'pattern',pat[r])
