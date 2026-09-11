# Step 3: certified branch constants m_{h,k} in
#   log Q8(tau) = -pi i tau + Log(c) + (e1-e3) log X + Log(1+u) - Log(1+v) + 2 pi i m,
# tau = h/k + i z/k^2 with z = 1 (so Im tau' = 1), X = e^{-pi/(32 z)} real.
from flint import arb, acb, ctx
from math import gcd, log
ctx.prec=160
PI=arb.pi(); I=acb(0,1)
def psi(a,j):
    r=a-4
    if (j-r)%8: return 0
    return -1 if ((j-r)//8)%2 else 1
def G(a,h,k,l):
    M=16*k; s_=acb(0)
    for s in range(M):
        p=psi(a,s)
        if p: s_+=p*acb.exp_pi_i(acb(arb(2*((h*s*s+l*s)%M))/M))
    return s_
B=8000
b=[0]*(B+1)
for d in range(1,B+1):
    r=d%8; s=-d if r in (1,7) else (d if r in (3,5) else 0)
    if s:
        for m in range(d,B+1,d): b[m]+=s
def logQ8(tau):
    q=(2*PI*I*tau).exp(); s=acb(0)
    for m in range(1,B+1):
        if b[m]: s+=(2*PI*I*tau*m).exp()*b[m]/m
    r=abs(q).upper(); rr=float(r)
    LB=1+log(B)
    tail=arb(rr)**B*(LB*arb(rr)/(1-arb(rr))+arb(rr)/(B*(1-arb(rr))**2))   # |b_m/m| <= 1+ln m
    return s+acb(arb(0,tail.upper()),arb(0,tail.upper()))
for (h,k) in [(3,8),(5,8),(3,16),(5,16),(11,16),(13,16)]:
    z=arb(1); tau=acb(arb(h)/k, z/k**2)
    X=(-PI/(32*z)).exp()
    LMAX=64
    Gh={}
    for a in (1,3):
        for l in range(0,LMAX+1):
            g=G(a,h,k,l)+(G(a,h,k,-l) if l else 0)
            Gh[(a,l)]=g
    # leading terms (Step 1: l=2 for a=1, l=6 for a=3)
    g1,g3=Gh[(1,2)],Gh[(3,6)]
    u=sum((Gh[(1,l)]/g1*X**(l*l-4) for l in range(3,LMAX+1)),acb(0))
    v=sum((Gh[(3,l)]/g3*X**(l*l-36) for l in range(7,LMAX+1)),acb(0))
    # truncation: |Ghat| <= 2*16k (trivial bound M), X^(l^2) for l>LMAX
    trunc=2*16*k*X**((LMAX+1)**2-36)*2
    u+=acb(arb(0,trunc.upper()),arb(0,trunc.upper())); v+=acb(arb(0,trunc.upper()),arb(0,trunc.upper()))
    assert abs(u)<0.5 and abs(v)<0.5
    c=g1/g3
    Lloc=-PI*I*tau+c.log()+(4-36)*X.log()+(1+u).log()-(1+v).log()
    Lex=logQ8(tau)
    diff=(Lex-Lloc)/(2*PI*I)
    m=round(float(diff.real.mid()))
    ok=abs(diff-m)<arb('0.25')
    print('h/k=%d/%d  Arg c/pi=%s  (Lex-Lloc)/(2 pi i)=%s  -> m=%d certified:%s'%(h,k,(c.arg()/PI).str(8),diff.str(6),m,ok))
