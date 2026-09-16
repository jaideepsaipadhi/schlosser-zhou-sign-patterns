# delta in {1,3,5}: large-n dominance at n>=n1 (ratios decreasing in n) + exact integer coefficients for n<n1.
import sys
from g7lib import *
def dom_check(d,r,n1,use_corr):
    n=n1+((r-n1)%7); N=int(float((4*PI*n).sqrt().upper()))+2
    D=arb(d); C=2*n-D/2; B=D/2; B2=D/2-2
    if not use_corr:
        a=amp_main(7,r,D); dom=abs(a)*Plow(7,B,C,C); sg=a
        comp=(6*D*P(7,B2,C) if d>4 else 0)+6*(P(14,B,C)+(D*P(14,B2,C) if d>4 else 0))+12*(P(21,B,C)+(D*P(21,B2,C) if d>4 else 0))
    else:
        a=amp_corr(7,r,D); dom=abs(a)*Plow(7,B2,C,C); sg=a
        comp=6*D*P(14,B2,C)+12*(P(21,B,C)+D*P(21,B2,C))
    comp+=(arb(N)/7)*(arb(6)*N/7)*(P(28,B,C)+(D*P(28,B2,C) if d>4 else 0))+Econ(D,N)
    return bool(dom>comp),(1 if sg>0 else -1),(dom/comp).str(3)
def ints(d,NMAX):
    f=[0]*(NMAX+1); f[0]=1
    for m in range(1,NMAX+1):
        for _ in range(d):
            for n in range(NMAX,m-1,-1): f[n]-=f[n-m]            # times (1-q^m)^d
        if m%7==0:
            for _ in range(d):
                for n in range(m,NMAX+1): f[n]+=f[n-m]        # divide by (1-q^m)^d  (q^7 factors)
    return f
cases={1:('+--00+0',[3,4,6]),3:('+-0+00-',[2,4,5]),5:('+-++---',[])}
n1=int(sys.argv[1])
for d,(pat,zeros) in cases.items():
    res=[]
    for r in range(7):
        if r in zeros: continue
        uc=(d==5 and r in (1,2,6))
        ok,sg,mg=dom_check(d,r,n1,uc); want=1 if pat[r]=='+' else -1
        res.append((r,ok and sg==want,mg))
    f=ints(d,n1)
    bad=[n for n in range(n1+1) if (pat[n%7]=='0' and f[n]!=0) or (pat[n%7]!='0' and f[n]!=0 and (1 if f[n]>0 else -1)!=(1 if pat[n%7]=='+' else -1))]
    z=[n for n in range(n1+1) if f[n]==0 and pat[n%7]!='0']
    print('delta=%d: large-n residues %s ; exact n<=%d violations %s, extra zeros %s'%(d,res,n1,bad[:6],z[:6]))
