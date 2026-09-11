# Tail n>20000 for Conj 24, 2<=delta<=3: all competitor/dominant ratios are C*y^p*e^{-g y}, g>0, y=sqrt(2(n+delta)).
from flint import arb, ctx
ctx.prec=100
exec(open('q12_pos.py').read().split('NLO,NHI,S=')[0])
yl=(2*arb(20000)+4).sqrt(); yh=(2*arb(20000)+6).sqrt()
N=(2*PI).sqrt()*yh+3
U_=lambda x:x.exp()/(2*PI*x).sqrt()
def env_low(k,B):   # lower bound of P_k(B) at y in [yl,yh] (C=y^2)
    x=(2*PI/k)*B.sqrt()*yl; kap=(2*PI*x).sqrt()*(-x).exp()*x.bessel_i(1)
    return (2*PI/k)*B.sqrt()/yh*kap*U_(x)
def env_up(k,B): return (2*PI/k)*B.sqrt()/yl*U_((2*PI/k)*B.sqrt()*yh)
ok=True
for jj in range(40):
    lo=arb(2)+arb(jj)/40; hi=arb(2)+arb(jj+1)/40
    E=Econ(hi,int(float(N.upper()))+1)
    cnt=(N/4)*(N/2)
    # degenerate residues: dominant delta*P(12,2(delta-1)), worst: delta=2 (smallest B)  -- competitor uppers at delta=3
    D=lo*env_low(12,2*(lo-1))
    r_c2=hi*(hi+1)*env_up(12,2*(hi-2))/D
    r_24=(4*env_up(24,2*hi)+4*hi*env_up(24,2*(hi-1))+2*hi*(hi+1)*env_up(24,2*(hi-2)))/D
    r_k=(6+cnt)*(env_up(36,2*hi)+hi*env_up(36,2*(hi-1))+hi*(hi+1)/2*env_up(36,2*(hi-2)))/D
    r_E=E/D
    _=('degenerate: ratios corr2 %s k24 %s k>=36 %s E %s  sum<1: %s'%(r_c2.str(3),r_24.str(3),r_k.str(3),r_E.str(3),(r_c2+r_24+r_k+r_E)<1))
    # nondegenerate: dominant 1*P(12,2delta) (|a0|>=1), competitors: corr1 2*delta*P(12,2(delta-1)) etc.
    D2=env_low(12,2*lo)
    s=(2*hi*env_up(12,2*(hi-1))+hi*(hi+1)*env_up(12,2*(hi-2))+4*env_up(24,2*hi)+4*hi*env_up(24,2*(hi-1))+2*hi*(hi+1)*env_up(24,2*(hi-2)))/D2
    s+=(6+cnt)*(env_up(36,2*hi)+hi*env_up(36,2*(hi-1))+hi*(hi+1)/2*env_up(36,2*(hi-2)))/D2+E/D2
    _=('nondegenerate: sum of ratios %s <1: %s'%(s.str(3),s<1))
    _=('gaps (must be >0): deg vs corr2 %s, deg vs k24 %s, nondeg vs corr1 %s'%(((PI/6)*(2*(lo-1)).sqrt()-(PI/6)*(2*(hi-2)).sqrt()).str(3),((PI/6)*(2*(lo-1)).sqrt()-(PI/12)*(2*hi).sqrt()).str(3),((PI/6)*(2*lo).sqrt()-(PI/6)*(2*(hi-1)).sqrt()).str(3)))
    ok=ok and (r_c2+r_24+r_k+r_E)<1 and s<1
    if jj in (0,39): print(float(lo),'deg sum',(r_c2+r_24+r_k+r_E).str(3),'nondeg sum',s.str(3))
print('tail n>20000 certified on all 40 subintervals of [2,3]:',ok)
