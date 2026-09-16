#!/usr/bin/env python3
"""
Rigorous (ball-arithmetic) computation of coefficients of Q8(q)^delta,
Q8 = (q,q^7;q^8)_inf / (q^3,q^5;q^8)_inf  (Schlosser-Zhou Conjecture 21).

Q8^delta = exp(delta * sum_k b_k q^k / k),  b_k = sum_{d|k} (+d if d=3,5 mod 8; -d if d=1,7 mod 8).
Uses arb_series.exp (Newton iteration, quasi-linear) with every coefficient a certified ball.

Usage:  python3 q8_signcheck.py N DELTA PREC_BITS [test]
  e.g.  python3 q8_signcheck.py 3000 2.66448 400 test   (quick self-check vs naive recurrence)
        python3 q8_signcheck.py 860000 2.66448 3600      (the real run)
Prints sign and ball radius of f_n for n = 2,6 mod 8 in the last ~2000 indices,
and flags any coefficient whose certified sign contradicts the conjectured pattern
(+ - + + - + - -) for n mod 8 = 0..7.
"""
import sys, time
import flint
from flint import arb, arb_series, ctx

N=int(sys.argv[1]); DELTA=sys.argv[2]; PREC=int(sys.argv[3]); TEST=len(sys.argv)>4
ctx.prec=PREC
ctx.cap=N+1
delta=arb(DELTA)          # exact decimal -> ball (tiny radius)

t=time.time()
b=[0]*(N+1)
for d in range(1,N+1):
    r=d%8
    s=-d if r in (1,7) else (d if r in (3,5) else 0)
    if s:
        for m in range(d,N+1,d): b[m]+=s
L=arb_series([arb(0)]+[delta*arb(b[k])/k for k in range(1,N+1)], prec=N+1)
F=L.exp()
c=F.coeffs()
print('series done in %.1fs'%(time.time()-t), flush=True)

pattern='+-++-+--'
bad=[]; undecided=[]
for n in range(len(c)):
    x=c[n]
    if x>0: s='+'
    elif x<0: s='-'
    else:
        if x.is_exact() and x==0: s='0'
        else: undecided.append(n); continue
    if s!='0' and s!=pattern[n%8]: bad.append(n)
print('undecided signs (ball contains 0):',undecided[:50],'count',len(undecided))
print('certified pattern violations:',bad[:50],'count',len(bad))
for n in range(max(0,N-200),N+1):
    if n%8 in (2,6):
        x=c[n]; print(n, n%16, 'sign', '+' if x>0 else ('-' if x<0 else '?'), 'mid~', x.mid().str(5), 'rad', x.rad().str(3))

if TEST:
    from fractions import Fraction
    import mpmath as mp
    mp.mp.dps=PREC//3
    dd=mp.mpf(DELTA); f=[mp.mpf(1)]+[0]*N
    for n in range(1,N+1): f[n]=dd*mp.fsum(b[k]*f[n-k] for k in range(1,n+1))/n
    err=max(abs(mp.mpf(c[n].mid().str(PREC//4, radius=False))-f[n])/abs(f[n]) for n in range(1,N+1) if f[n]!=0)
    print('self-check: max relative deviation vs naive recurrence =', err)
