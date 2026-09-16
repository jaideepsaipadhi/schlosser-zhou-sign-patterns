# Exact finite check for Q12^delta on a closed delta-interval [A,B] for n<=NMAX (pattern given), with
# exact division by linear factors at endpoint zeros (sign of the factor accounted for).
import sys
from flint import fmpq_poly, fmpq
src=open('c3.py').read()
exec(src[src.index('def sgn(v)'):src.index('bad=[]')])
NMAX=int(sys.argv[1]); A=fmpq(*map(int,sys.argv[2].split('/'))); B=fmpq(*map(int,sys.argv[3].split('/'))); pat=sys.argv[4]
eps=lambda d:(1 if d%12 in (1,11) else (-1 if d%12 in (5,7) else 0))
b=[0]*(NMAX+1)
for d in range(1,NMAX+1):
    e=eps(d)
    if e:
        for m in range(d,NMAX+1,d): b[m]-=e*d
X=fmpq_poly([0,1]); c=[fmpq_poly([1])]
for n in range(1,NMAX+1):
    acc=fmpq_poly([0])
    for k in range(1,n+1):
        if b[k]: acc+=b[k]*c[n-k]
    c.append(acc*X/n)
L=len(pat); bad=[]; notes=[]
for n in range(NMAX+1):
    p=c[n]; want=1 if pat[n%L]=='+' else -1; fac=1
    for e,s in ((A,1),(B,-1)):          # factor (x-A)>=0 on [A,B]; (x-B)<=0
        while p(e)==0 and not p.is_zero():
            q,r=divmod(p,fmpq_poly([-e,1])); assert r==0; p=q; fac*=s; notes.append((n,float(e)))
    target=want*fac
    issue=[]
    for e in (A,B):
        if sgn(p(e))!=target: issue.append(('end',float(e)))
    rr=roots_in(p,A,B)
    if rr: issue.append(('interior',rr))
    if issue: bad.append((n,issue))
print('n<=%d on [%s,%s], pattern %s: endpoint zeros divided %s; problems: %s'%(NMAX,A,B,pat,notes[:10],bad[:8]))
