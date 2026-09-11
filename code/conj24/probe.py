# Sign-pattern probe for P(q)^delta, P = prod (1-q^d)^{eps(d)}; CDQ ball-arithmetic recurrence.
import sys
from flint import arb, arb_series, ctx
def coeffs(eps, N, delta, prec):
    ctx.prec=prec; ctx.cap=N+1
    b=[0]*(N+1)
    for d in range(1,N+1):
        e=eps(d)
        if e:
            for m in range(d,N+1,d): b[m]-=e*d
    dl=arb(delta); f=[None]*(N+1); f[0]=arb(1); BS=2048; LEAF=64
    def solve(l,r):
        if r-l<=LEAF:
            for n in range(l,r):
                if n==0: continue
                s=f[n] if f[n] is not None else arb(0)
                for j in range(l,n):
                    if b[n-j]: s+=b[n-j]*f[j]
                f[n]=dl*s/n
            return
        m=(l+r)//2; solve(l,m)
        for s0 in range(l,m,BS):
            s1=min(s0+BS,m); A=arb_series(f[s0:s1])
            for t0 in range(m,r,BS):
                t1=min(t0+BS,r); d0=t0-(s1-1); d1=(t1-1)-s0
                C=(A*arb_series(b[d0:d1+1])).coeffs()
                for n in range(t0,t1):
                    i=n-s0-d0
                    if 0<=i<len(C): f[n]=C[i] if f[n] is None else f[n]+C[i]
        solve(m,r)
    solve(0,N+1)
    return f
def signs(f):
    out=[]
    for x in f:
        out.append('+' if x>0 else ('-' if x<0 else ('0' if (x.is_exact() and x==0) else '?')))
    return ''.join(out)
def check(sg, pat):
    L=len(pat); bad=[n for n,s in enumerate(sg) if s!='?' and s!='0' and pat[n%L]!='*' and s!=pat[n%L]]
    zeros=[n for n,s in enumerate(sg) if s=='0']; und=[n for n,s in enumerate(sg) if s=='?']
    return bad,zeros,und
G=lambda p:(lambda d:1 if d%p else 0)
Q12=lambda d:(1 if d%12 in (1,11) else (-1 if d%12 in (5,7) else 0))
