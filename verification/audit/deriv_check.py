# Audit of the derivative error bound for G7 (window near delta_c): compare exact d/d delta c_n with the
# derivative of the explicit k=7,14,21 terms; the difference must be <= derivative of k>=28 terms + E'(s,N).
import sys, mpmath as mp, math
sys.path.insert(0,'/home/claude/repo/code/conj20')
import g7lib
from flint import arb
mp.mp.dps=60
def exact_fg(N,d):
    b=[0]*(N+1)
    for dd in range(1,N+1):
        if dd%7:
            for m in range(dd,N+1,dd): b[m]-=dd
    d=mp.mpf(d); f=[mp.mpf(1)]+[0]*N; g=[mp.mpf(0)]+[0]*N
    for n in range(1,N+1):
        s1=mp.fsum(b[k]*f[n-k] for k in range(1,n+1)); s2=mp.fsum(b[k]*g[n-k] for k in range(1,n+1))
        f[n]=d*s1/n; g[n]=(s1+d*s2)/n
    return f,g
def fr(x): return mp.mpf(x.numerator)/x.denominator
def explicit(d,n,kmax):
    D=mp.mpf(d); C=2*n-D/2; B=D/2; B2=D/2-2
    P=lambda k,B: (2*mp.pi/k)*mp.sqrt(B/C)*mp.besseli(1,(2*mp.pi/k)*mp.sqrt(B*C)) if B>0 else 0
    M=0
    for k in range(7,kmax+1,7):
        for h,sh,hp in g7lib.shifts(k):
            ph=-mp.pi*D*fr(sh)-2*mp.pi*h*n/k
            M+=mp.cos(ph)*P(k,B)-D*mp.cos(ph-2*mp.pi*hp/k)*P(k,B2)
    return M
d0='4.8735'; h=mp.mpf('1e-20')
f,g=exact_fg(400,d0)
exec(open('/home/claude/repo/code/conj20/g7_deriv.py').read().split('def check')[0])
worst=0
for n in range(100,401,7):
    der_expl=(explicit(mp.mpf(d0)+h,n,21)-explicit(mp.mpf(d0)-h,n,21))/(2*h)
    N=int(mp.ceil(mp.sqrt(4*mp.pi*n)))+2
    # derivative bound for the remaining explicit k>=28 terms + error derivative (as in g7_deriv.py)
    Bb=arb(d0); A=arb(d0); Clo=2*n-Bb/2
    t=(-2*PI).exp(); F=f_prod(t)**Bb*f_prod(t**7)**Bb; l7=(f_prod(t)*f_prod(t**7)).log(); S=arb(6)/49
    arcs=2*S*E1*(PI*Bb/2).exp()*((PI*N/10)*(PI/R2)*N+(2+PI/2)*(1/(A/2)+Bb*t/(A/2-2))*(1+(arb(N)/7).log())+10)
    chords=2*R2*S*E1*(PI*Bb/2).exp()*((PI*N/10+PI*N*N+arb('0.1'))*(F-1-Bb*t)+F*l7-t)*N
    ng=2*R2*E1*K1**Bb*(arb('0.556')*N*N+arb('0.449')*N*N+2+K1.log())
    rest=(arb(N)/7)*(arb(6)*N/7)*((PI*N/10+D(28,A/2,Clo))*P(28,Bb/2,Clo)+(1+Bb*PI*N/10+Bb*D(28,A/2-2,Clo))*P(28,Bb/2-2,Clo))
    bound=float((arcs+chords+ng+rest).mid())
    worst=max(worst,float(abs(g[n]-der_expl))/bound)
print('G7 d/d delta audit at delta=%s, 100<=n<=400: max |exact - explicit derivative| / bound = %.3e (must be <=1)'%(d0,worst))
