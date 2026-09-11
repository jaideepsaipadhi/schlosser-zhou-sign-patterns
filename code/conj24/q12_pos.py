# Conj 24, 2<=delta<=3: certified sign criterion for NLO<=n<=NHI (Q12 analogue of Prop 4.3 / 5.1).
import sys, time
from flint import arb, ctx
ctx.prec=100
PI=arb.pi(); E1=arb(1).exp(); R2=arb(2).sqrt()
K12=arb('550.75')
I1=lambda x:x.bessel_i(1)
t=(-2*PI).exp(); U=arb(0); V=arb(0)
for J in range(3,600):
    if J%24 in (2,22) and J>2: U+=t**(arb(J*J-4)/96)
    if J%24 in (10,14) and J>10: V+=t**(arb(J*J-100)/96)
def Frem(d):   # F~_d(t)-1-d t-d(d+1)/2 t^2  at t=e^{-2pi}
    return (1-U)**(-d)*(1-V)**(-d)-1-d*t-d*(d+1)/2*t*t
def Econ(d,N):
    S=(arb(1)/2+arb(N)/8)/(N+1)
    a=R2*PI*S*E1*(2*PI*d).exp()*(1+d*t+d*(d+1)/2*t*t)
    b=2*R2*S*E1*(2*PI*d).exp()*Frem(d)
    c=2*R2*E1*K12**d
    return a+b+c
def P(k,B,C):
    if not (B>0): return arb(0)
    return (2*PI/k)*(B/C).sqrt()*I1((2*PI/k)*(B*C).sqrt())
def Pl(k,Bl,Cl,Ch): return (2*PI/k)*(Bl/Ch).sqrt()*I1((2*PI/k)*(Bl*Cl).sqrt()) if Bl>0 else arb(0)
def Pu(k,Bh,Cl,Ch): return (2*PI/k)*(Bh/Cl).sqrt()*I1((2*PI/k)*(Bh*Ch).sqrt()) if Bh>0 else arb(0)
def cospi(x): return arb.cos_pi(arb(x))
NLO,NHI,S=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])
pts=[arb(2)+arb(i)/S for i in range(S+1)]
subs=[(pts[i],pts[i+1]) for i in range(S)]
fails=[]; t0=time.time(); Ec={}
for n in range(NLO,NHI+1):
    for a,b in subs:
        Cl=2*(n+a); Ch=2*(n+b)
        N=int(float((4*PI*(n+b)).sqrt().upper()))+2
        key=(float(b),N)
        if key not in Ec: Ec[key]=Econ(b,N)
        E=Ec[key]
        a0=abs(2*cospi(arb(5*n)/6)); a1=abs(2*cospi(arb(5*n-1)/6))    # times delta
        c2u=b*(b+1)*Pu(12,2*(b-2),Cl,Ch)                                 # |corr2| <= d(d+1) P(12,2(d-2))
        T24=4*Pu(24,2*b,Cl,Ch)+4*b*Pu(24,2*(b-1),Cl,Ch)+2*b*(b+1)*Pu(24,2*(b-2),Cl,Ch)
        Tk=arb(0)
        Tk+=6*(Pu(36,2*b,Cl,Ch)+b*Pu(36,2*(b-1),Cl,Ch)+b*(b+1)/2*Pu(36,2*(b-2),Cl,Ch))
        cntsum=sum(k//2 for k in range(40,N+1,4))                 # P_k decreasing in k: bound by k=40 term
        Tk+=cntsum*(Pu(40,2*b,Cl,Ch)+b*Pu(40,2*(b-1),Cl,Ch)+b*(b+1)/2*Pu(40,2*(b-2),Cl,Ch))
        if n%12 in (3,9):
            m=a*a1*Pl(12,2*(a-1),Cl,Ch)-c2u-T24-Tk-E
        else:
            m=a0*Pl(12,2*a,Cl,Ch)-b*a1*Pu(12,2*(b-1),Cl,Ch)-c2u-T24-Tk-E
        if not (m>0): fails.append((n,float(a)))
print('n=%d..%d, %d subintervals: %.1fs, failures %d, last %s'%(NLO,NHI,S,time.time()-t0,len(fails),fails[-4:]))
