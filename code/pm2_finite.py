# Exact integer coefficients of Q8^{+2}, Q8^{-2} up to NMAX and comparison with the length-16 patterns of [SZ, Conj. 21].
NMAX=400
def series(e):   # product over m of (1-q^m)^{eps_m * e}, e=+-2, exact integers
    f=[0]*(NMAX+1); f[0]=1
    for m in range(1,NMAX+1):
        r=m%8; s=1 if r in (1,7) else (-1 if r in (3,5) else 0)
        p=s*e
        if p==0: continue
        for _ in range(abs(p)):
            if p>0:   # multiply by (1-q^m)
                for n in range(NMAX,m-1,-1): f[n]-=f[n-m]
            else:     # divide by (1-q^m)
                for n in range(m,NMAX+1): f[n]+=f[n-m]
    return f
for e,pat in ((2,'+-++-+--+--+-++-'),(-2,'+++++----+++----')):
    f=series(e)
    bad=[n for n in range(NMAX+1) if f[n]==0 or ('+' if f[n]>0 else '-')!=pat[n%16]]
    print('delta=%+d: coefficients n<=%d violating/zero vs pattern %s: %s'%(e,NMAX,pat,bad[:10]))
