# Dominance certification: for delta in [a,b] and ALL n >= n1 with n = r mod 7, sign(c_delta(n)) = sign(k=7 main term).
# Checked at n1 with exact Bessel values; all competitor/dominant ratios are of the form C*y^p*exp(-g*y), g>0 (decreasing).
import sys
from g7lib import *
def check(a,b,n1,r,pat):
    A=arb(str(a)); Bb=arb(str(b)); big=bool(Bb>4)
    d=arb(a).union(arb(b)) if hasattr(arb,'union') else None
    m=(A+Bb)/2; rad=(Bb-A)/2
    dball=arb(m.mid(),float(rad.mid())*1.0000001+1e-40)
    n=n1+((r-n1)%7)
    Clo=2*n-Bb/2; Chi=2*n-A/2
    amp=amp_main_mv(7,r,dball)
    want=1 if pat[r]=='+' else -1
    if not ((amp>0 and want>0) or (amp<0 and want<0)): return False,'amp sign',amp.str(4)
    amin=abs(amp).lower() if hasattr(abs(amp),'lower') else abs(amp)
    amin=arb(abs(amp).lower())
    B=A/2
    dom=amin*Plow(7,B,Clo,Chi)
    Bh=Bb/2; B2=Bh-2
    comp=arb(0)
    if big: comp+=6*Bb*P(7,B2,Clo)*(Chi/Clo).sqrt()
    N=int(float((4*PI*n).sqrt().upper()))+2
    for k,ph in ((14,6),(21,12)):
        comp+=ph*(P(k,Bh,Clo)+(Bb*P(k,B2,Clo) if big else 0))*(Chi/Clo).sqrt()
    comp+=(arb(N)/7)*(arb(6)/7)*2*PI*(Bh/Clo).sqrt()*((2*PI/28)*(Bh*Chi).sqrt()).bessel_i(1)*(1+Bb)
    comp+=Econ(Bb,N)
    return bool(dom>comp),'margin',(dom/comp).str(3)
if __name__=='__main__':
    a,b,n1=float(sys.argv[1]),float(sys.argv[2]),int(sys.argv[3]); pat=sys.argv[4]
    print([ (r,)+check(a,b,n1,r,pat) for r in range(7)])
