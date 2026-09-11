from flint import fmpq_poly, fmpq
N=400
b=[0]*(N+1)
for e in range(1,N+1):
    r=e%8; s=-e if r in (1,7) else (e if r in (3,5) else 0)
    if s:
        for m in range(e,N+1,e): b[m]+=s
X=fmpq_poly([0,1]); c=[fmpq_poly([1])]
for n in range(1,N+1):
    acc=fmpq_poly([0])
    for k in range(1,n+1):
        if b[k]: acc+=b[k]*c[n-k]
    c.append(acc*X/n)
m1=fmpq(-1)
val=[c[n](m1) for n in range(N+1)]
der=[c[n].derivative()(m1) for n in range(N+1)]
for r in range(8):
    zs=[n for n in range(r,N+1,8) if val[n]==0]
    nz=[n for n in range(r,N+1,8) if val[n]!=0]
    sv=set((v>0)-(v<0) for v in [val[n] for n in nz])
    sd=set((d>0)-(d<0) for d in [der[n] for n in zs])
    print('res',r,' c_{-1}=0 at',len(zs),'of',len(range(r,N+1,8)),' nonzero idx',nz[:5],' sign(c_-1) on nonzero',sv,' sign(dc/dδ) on zeros',sd, ' first zero-der', [n for n in zs if der[n]==0][:5])
