# Numerical model of Q12^delta near its growing cusps (Poisson form), with branch fixed by the log-series.
import mpmath as mp, math, sys
mp.mp.dps=40
PI=mp.pi
def psi(a,j):
    r=a-6
    if (j-r)%12: return 0
    return -1 if ((j-r)//12)%2 else 1
Bm=6000
b=[0]*(Bm+1)
for d in range(1,Bm+1):
    e=1 if d%12 in (1,11) else (-1 if d%12 in (5,7) else 0)
    if e:
        for m in range(d,Bm+1,d): b[m]-=e*d
def logQ(q): return mp.fsum(b[m]*q**m/m for m in range(1,Bm+1) if b[m])
def Ghat(a,h,k,e):
    l=math.isqrt(e); M=24*k
    def G(l): return mp.fsum(psi(a,s)*mp.expjpi(mp.mpf(2*((h*s*s+l*s)%M))/M) for s in range(M) if psi(a,s))
    return G(l)+(G(-l) if l else 0)
def local(h,k,delta,sign=+1,EMAX=900):
    # returns list of (A, C) with Q12^delta ~ sum C * exp(pi*A/z) * e(-delta*tau)-part handled outside
    a1,a5=(1,5) if sign>0 else (5,1)          # numerator/denominator of R^{+-1}
    lead1,lead5=(4,100)
    W1={e:Ghat(1,h,k,e) for e in range(0,EMAX) if math.isqrt(e)**2==e}
    W5={e:Ghat(5,h,k,e) for e in range(0,EMAX) if math.isqrt(e)**2==e}
    W1={e:c for e,c in W1.items() if abs(c)>1e-20}; W5={e:c for e,c in W5.items() if abs(c)>1e-20}
    e1=min(W1); e5=min(W5)
    num,den=(W1,W5) if sign>0 else (W5,W1)
    en,ed=min(num),min(den)
    # series in X: num/den normalized
    sn={e-en:c/num[en] for e,c in num.items()}; sd={e-ed:c/den[ed] for e,c in den.items()}
    maxE=96*4
    def lg(s):
        v={e:c for e,c in s.items() if 0<e<=maxE}; out={}; p={0:mp.mpf(1)}
        for m in range(1,12):
            np_={}
            for x,cx in p.items():
                for y,cy in v.items():
                    if x+y<=maxE: np_[x+y]=np_.get(x+y,0)+cx*cy
            p=np_
            for e,c in p.items(): out[e]=out.get(e,0)+(-1)**(m+1)*c/m
        return out
    L={e:lg(sn).get(e,0)-lg(sd).get(e,0) for e in set(lg(sn))|set(lg(sd))}
    c0=num[en]/den[ed]
    # branch: log Q12 = -2 pi i tau*sign... sign*log Q12 = sign*(-2 pi i tau) + Log c0 + (en-ed)*log X + L
    z=mp.mpf('0.3'); tau=mp.mpf(h)/k+1j*z/k**2; X=mp.exp(-PI/(48*z)); q=mp.exp(2j*PI*tau)
    lhs=sign*logQ(q); rhs=sign*(-2j*PI*tau)+mp.log(c0)+(en-ed)*mp.log(X)+mp.fsum(c*X**e for e,c in L.items())
    m=mp.nint(mp.im(lhs-rhs)/(2*PI)); assert abs(lhs-rhs-2j*PI*m)<0.05
    D=abs(delta); lc=mp.log(c0)+2j*PI*m
    E={0:mp.mpf(1)}; term={0:mp.mpf(1)}
    for j in range(1,40):
        nt={}
        for x,cx in term.items():
            for y,cy in L.items():
                if x+y<=maxE: nt[x+y]=nt.get(x+y,0)+cx*cy*D/j
        term=nt
        for e,c in term.items(): E[e]=E.get(e,0)+c
    out=[]
    for j,c in sorted(E.items()):
        tot=D*(en-ed)+j
        if tot<0 and abs(c)>1e-25: out.append((-tot/48, mp.exp(D*lc)*c))   # X^tot = exp(pi*(-tot/48)/z)
    return out
def predict(n,delta,cusps,sign=+1):
    D=abs(delta); tot=0; C=2*(n+sign*D)
    for (h,k),terms in cusps.items():
        ph=mp.expjpi(-2*n*mp.mpf(h)/k)*mp.expjpi(-2*sign*D*mp.mpf(h)/k)
        for A,Cf in terms:
            tot+=Cf*ph*(2*PI/k)*mp.sqrt(A/C)*mp.besseli(1,(2*PI/k)*mp.sqrt(A*C))
    return tot
