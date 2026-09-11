# Tail lemma (n>20000, y>=200), self-contained constants.
from flint import arb, ctx
ctx.prec=100
from consts import Econst, PI, R2
y=arb(200); lo=arb(8)/3; hi=arb(4)
a8c_min=(PI/4)*(lo-2).sqrt()
x0=a8c_min*y/2
kap=(2*PI*x0).sqrt()*(-x0).exp()*x0.bessel_i(1)
U=lambda x:x.exp()/(2*PI*x).sqrt()
L=lambda x:kap*U(x)
Emax=Econst(hi,501)                          # N >= sqrt(4 pi 20000) > 501 ; E decreasing in N
nk4=((2*PI).sqrt()*y+2)/4                    # #{24<=k<=N, 4|k} <= N/4, N<=sqrt(2pi) y + 1
dh=hi
a16c=(2*PI/16)*(dh-2).sqrt(); a24=(2*PI/24)*dh.sqrt(); a24c=(2*PI/24)*(dh-2).sqrt()
# residues 2,6
D=R2*lo*(PI/4)*(lo-2).sqrt()*L(a8c_min*y)/y
R16main=arb('0.375')
R16c=4*dh*(2*PI/16)*(dh-2).sqrt()*U(a16c*y)/y/D
R24=nk4*PI*(dh.sqrt()*U(a24*y)+dh*(dh-2).sqrt()*U(a24c*y))/y/D
RE=Emax/D
print('res 2,6: 0.375 + %s + %s + %s  <1: %s'%(R16c.str(3),R24.str(3),RE.str(3),(R16main+R16c+R24+RE)<1))
print('gaps', (a8c_min-a16c).str(3), (a8c_min-a24).str(3))
a8_min=(PI/4)*lo.sqrt()
D2=R2*(PI/4)*lo.sqrt()*L(a8_min*y)/y
gapC=(PI/4)*(hi.sqrt()-(hi-2).sqrt())
RC8=R2*(dh*(dh-2)).sqrt()*(1/kap)*(dh.sqrt()/(dh-2).sqrt()).sqrt()*(-gapC*y).exp()
RT16=(4*(2*PI/16)*dh.sqrt()*U((2*PI/16)*dh.sqrt()*y)+4*dh*(2*PI/16)*(dh-2).sqrt()*U(a16c*y))/y/D2
RT24=nk4*PI*(dh.sqrt()*U(a24*y)+dh*(dh-2).sqrt()*U(a24c*y))/y/D2
RE2=Emax/D2
print('other: %s + %s + %s + %s <1: %s  kappa=%s'%(RC8.str(3),RT16.str(3),RT24.str(3),RE2.str(3),(RC8+RT16+RT24+RE2)<1,kap.str(5)))
