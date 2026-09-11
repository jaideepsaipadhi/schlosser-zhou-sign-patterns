exec(open('neg_finite.py').read().split("A=fmpq(-1); B=fmpq(-99,100)")[0])
A=fmpq(-1); B=fmpq(-772,1000)
sig='+++----+'
quad=fmpq_poly([-6,-7,1])      # delta^2-7delta-6, roots (7 +- sqrt73)/2
bad=[]; divided=[]
for n in range(NMAX+1):
    p=c[n]; want=1 if sig[n%8]=='+' else -1
    if n==4:
        q,r=divmod(p,quad); assert r==0; p=q*fmpq_poly([1,1]) if False else q
        # c(4) = quad*q ; quad>=0 on [-1,b], =0 only at b; so sign c(4) = sign q there
    if p(A)==0:
        q,r=divmod(p,fmpq_poly([1,1])); assert r==0; p=q; divided.append(n)
    issue=[]
    for e in (A,B):
        if sgn(p(e))!=want: issue.append(('end',float(e),sgn(p(e))))
    rr=roots_in(p,A,B)
    if rr: issue.append(('interior',rr))
    if issue: bad.append((n,issue))
print('n<=%d on [-1,-0.772]: divided at %d indices; problems: %s'%(NMAX,len(divided),bad[:10]))
