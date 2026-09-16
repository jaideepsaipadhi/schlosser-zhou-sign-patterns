# Step 1 (negative delta): exact Gauss-sum ratios and certified branch constants at the
# cusps growing for Q8^{-d'}: h/k in {1/8,7/8,1/16,7/16,9/16,15/16}.
exec(open('exact.py').read().split('for k in (8,16):')[0])
from flint import arb, acb, ctx
ctx.prec=160
PI=arb.pi(); I=acb(0,1)
B=8000
b=[0]*(B+1)
for d in range(1,B+1):
    r=d%8; s=-d if r in (1,7) else (d if r in (3,5) else 0)
    if s:
        for m in range(d,B+1,d): b[m]+=s
from math import log
def logQ8(tau):
    s=acb(0)
    for m in range(1,B+1):
        if b[m]: s+=(2*PI*I*tau*m).exp()*b[m]/m
    rr=float(abs((2*PI*I*tau).exp()).upper()); LB=1+log(B)
    tail=arb(rr)**B*(LB*arb(rr)/(1-arb(rr))+arb(rr)/(B*(1-arb(rr))**2))
    return s+acb(arb(0,tail.upper()),arb(0,tail.upper()))
def zacb(v,k):
    M=16*k; s=acb(0)
    for i,c in enumerate(v):
        if c: s+=c*acb.exp_pi_i(acb(arb(2*i)/M))
    return s
res={}
for (h,k) in [(1,8),(7,8),(1,16),(7,16),(9,16),(15,16)]:
    M=16*k
    g3=Ghat(3,h,k,4); g1=Ghat(1,h,k,36)
    assert not iszero(g3) and not iszero(g1)
    # support check: a=3 leads at l=2, a=1 at l=6
    s3=[l for l in range(0,20) if not iszero(Ghat(3,h,k,l*l))]; s1=[l for l in range(0,20) if not iszero(Ghat(1,h,k,l*l))]
    t=[t for t in range(M) if mul_root(g1,t,k)==g3]      # c' = g3/g1 = zeta^t
    assert len(t)==1; t=t[0]
    # certified branch: -log Q8 = pi i tau + Log c' + pi/z + Log(1+u') - Log(1+v') + 2 pi i m'
    z=arb(1); tau=acb(arb(h)/k,z/k**2); X=(-PI/(32*z)).exp()
    Gc={}
    for a in (1,3):
        for l in range(0,65):
            Gc[(a,l)]=zacb(Ghat(a,h,k,l*l),k) if l else zacb(G(a,h,k,0),k)
    u=sum((Gc[(3,l)]/Gc[(3,2)]*X**(l*l-4) for l in range(3,65)),acb(0))
    v=sum((Gc[(1,l)]/Gc[(1,6)]*X**(l*l-36) for l in range(7,65)),acb(0))
    tr=2*M*X**(65*65-36)*2
    u+=acb(arb(0,tr.upper()),arb(0,tr.upper())); v+=acb(arb(0,tr.upper()),arb(0,tr.upper()))
    assert abs(u)<0.5 and abs(v)<0.5
    cp=acb.exp_pi_i(acb(arb(2*t)/M))
    L=PI*I*tau+cp.log()+PI/z+(1+u).log()-(1+v).log()
    d=(-logQ8(tau)-L)/(2*PI*I)
    m=round(float(d.real.mid())); ok=abs(d-m)<arb('0.25')
    argc=cp.arg()/PI
    phase=argc+2*m+arb(h)/k       # delta'-phase / pi  :  delta'*(Arg c' + 2 pi m' + pi h/k)
    res[(h,k)]=phase
    print('h/k=%d/%d supp(S3)=%s supp(S1)=%s  c\'=zeta_%d^%d  Arg/pi=%s  m\'=%d (certified %s)  phase/pi=%s'%(h,k,s3[:3],s1[:3],M,t,argc.str(6),m,ok,phase.str(6)))
