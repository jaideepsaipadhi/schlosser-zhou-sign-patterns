from gplib import *
import sys
def dom(d,r,n1):
    n=n1+((r-n1)%11); N=int(float((4*PI*n).sqrt().upper()))+2
    D=arb(d); C=2*n-5*D/6; B=5*D/6; B2=B-2; big=d>2.4
    a=arb(0)
    for h,sh,hp in shifts(11): a+=arb.cos_pi(-D*fr(sh)-arb(2*h*r)/11)
    d0=abs(a)*Plow(11,B,C,C)
    comp=(10*D*P(11,B2,C) if big else 0)+10*(P(22,B,C)+(D*P(22,B2,C) if big else 0))+20*(P(33,B,C)+(D*P(33,B2,C) if big else 0))
    comp+=(arb(N)/11)*(arb(10)*N/11)*(P(44,B,C)+(D*P(44,B2,C) if big else 0))+Econ(D,N)
    return bool(d0>comp),(1 if a>0 else -1),(d0/comp).str(3)
def ints(d,NM):
    f=[0]*(NM+1); f[0]=1
    for m in range(1,NM+1):
        for _ in range(d):
            for n in range(NM,m-1,-1): f[n]-=f[n-m]
        if m%11==0:
            for _ in range(d):
                for n in range(m,NM+1): f[n]+=f[n-m]
    return f
n1=int(sys.argv[1])
for d,pat,zeros in ((1,'+--0-+0+000',[3,6,8,9,10]),(3,'+-0+-0-000+',[2,5,7,8,9])):
    res=[(r,)+dom(d,r,n1) for r in range(11) if r not in zeros]
    ok=all(x[1] and x[2]==(1 if pat[x[0]]=='+' else -1) for x in res)
    f=ints(d,n1)
    bad=[n for n in range(n1+1) if (pat[n%11]=='0' and f[n]!=0) or (pat[n%11]!='0' and f[n]!=0 and (1 if f[n]>0 else -1)!=(1 if pat[n%11]=='+' else -1))]
    print('delta=%d large n>=%d all ok: %s margins %s; exact n<=%d violations %s extra zeros %s'%(d,n1,ok,[x[3] for x in res],n1,bad[:5],[n for n in range(n1+1) if f[n]==0 and pat[n%11]!='0'][:6]))
