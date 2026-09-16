#!/usr/bin/env python3
"""
Lighter certified check (independent of the circle-method analysis):
coefficients of Q8(q)^delta via the recurrence n f_n = delta * sum_{k<=n} b_k f_{n-k},
evaluated by divide-and-conquer (CDQ) with ball-arithmetic polynomial products.
Ball radii grow like the majorant prod(1-q^m)^(-delta) ~ exp(pi*sqrt(2*delta*n/3)), much
slower than for the Newton exponential, so far fewer bits (and ~1/4 the memory) are needed.

Usage: python3 q8_lite.py N DELTA PREC [WINDOW] [BLOCK]
  e.g. python3 q8_lite.py 40000 2.66448 1200            (quick test)
       python3 q8_lite.py 860000 2.66448 2600 400 8192   (the real run)
Reports certified signs of f_n for n = 2,6 (mod 8) in the last WINDOW indices and any
certified violation of the pattern (+,-,+,+,-,+,-,-).
"""
import sys, time
from flint import arb, arb_series, ctx
N=int(sys.argv[1]); D=sys.argv[2]; P=int(sys.argv[3]); W=int(sys.argv[4]) if len(sys.argv)>4 else 200
ctx.prec=P
delta=arb(D)
b=[0]*(N+1)
for e in range(1,N+1):
    r=e%8; s=-e if r in (1,7) else (e if r in (3,5) else 0)
    if s:
        for m in range(e,N+1,e): b[m]+=s
# memory-lean: f[n] holds the partial sum sum_k b_k f_{n-k} until it is finalised
f=[None]*(N+1); f[0]=arb(1)
ctx.cap=N+1
LEAF=64
BS=int(sys.argv[5]) if len(sys.argv)>5 else 8192
def solve(l,r):
    if r-l<=LEAF:
        for n in range(l,r):
            if n==0: continue
            s=f[n] if f[n] is not None else arb(0)
            for j in range(l,n):
                if b[n-j]: s+=b[n-j]*f[j]
            f[n]=delta*s/n
        return
    m=(l+r)//2
    solve(l,m)
    # blocked convolution: contribution of f[l..m) to f[m..r), peak memory O(BS)
    for s0 in range(l,m,BS):
        s1=min(s0+BS,m)
        A=arb_series(f[s0:s1])
        for t0 in range(m,r,BS):
            t1=min(t0+BS,r)
            d0=t0-(s1-1); d1=(t1-1)-s0
            Bs=arb_series(b[d0:d1+1])
            C=(A*Bs).coeffs()
            for n in range(t0,t1):
                i=n-s0-d0
                if 0<=i<len(C):
                    f[n]=C[i] if f[n] is None else f[n]+C[i]
            del Bs,C
        del A
    solve(m,r)
t=time.time()
solve(0,N+1)
print('done in %.1fs'%(time.time()-t),flush=True)
pat='+-++-+--'
bad=[];und=[]
for n in range(N+1):
    x=f[n]
    if x>0: s='+'
    elif x<0: s='-'
    else:
        if not (x.is_exact() and x==0): und.append(n)
        continue
    if s!=pat[n%8]: bad.append(n)
print('undecided:',len(und),und[:20]); print('certified violations:',len(bad),bad[:40])
for n in range(max(0,N-W),N+1):
    if n%8 in (2,6):
        x=f[n]; print(n,n%16,'+' if x>0 else ('-' if x<0 else '?'), x.mid().str(6), 'rad',x.rad().str(3))
