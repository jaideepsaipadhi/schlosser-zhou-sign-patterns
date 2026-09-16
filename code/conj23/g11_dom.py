# Dominance for G11^delta (delta <= 2.4, no correction): for delta in [a,b], all n>=n1 with n=r mod 11.
import sys
from gplib import *
def amp_mv(r,dball):
    m=arb(dball.mid()); rr=arb(dball.rad())
    a0=arb(0); der=arb(0)
    for h,sh,hp in shifts(11):
        a0+=arb.cos_pi(-m*fr(sh)-arb(2*h*r)/11)
        der+=arb.sin_pi(-dball*fr(sh)-arb(2*h*r)/11)*PI*fr(sh)
    return a0+der*arb(0,rr)
def check(a,b,n1,r,pat):
    A=arb(a); Bb=arb(b)
    db=arb(((A+Bb)/2).mid(),float(((Bb-A)/2).mid())*1.000001+1e-45)
    n=n1+((r-n1)%11); N=int(float((4*PI*n).sqrt().upper()))+2
    Clo=2*n-5*Bb/6; Chi=2*n-5*A/6
    amp=amp_mv(r,db); want=1 if pat[r]=='+' else -1
    if not ((amp>0 and want>0) or (amp<0 and want<0)): return False,'amp',amp.str(3)
    Blo=5*A/6; Bhi=5*Bb/6
    dom=arb(abs(amp).lower())*Plow(11,Blo,Clo,Chi)
    comp=10*P(22,Bhi,Clo)*(Chi/Clo).sqrt()+20*P(33,Bhi,Clo)*(Chi/Clo).sqrt()
    comp+=(arb(N)/11)*(arb(10)*N/11)*P(44,Bhi,Clo)*(Chi/Clo).sqrt()
    comp+=Econ(Bb,N)
    return bool(dom>comp),'margin',(dom/comp).str(3)
if __name__=='__main__':
    import mpmath as mp
    mp.mp.dps=40
    a,b=mp.mpf(sys.argv[1]),mp.mpf(sys.argv[2]); n1=int(sys.argv[3]); pat=sys.argv[4]
    out={'ok':0,'fail':[]}
    def tile(a,b,dep=0):
        if all(check(mp.nstr(a,40),mp.nstr(b,40),n1,r,pat)[0] for r in range(11)): out['ok']+=1; return
        if dep>40: out['fail'].append((a,b)); return
        m=(a+b)/2; tile(a,m,dep+1); tile(m,b,dep+1)
    tile(a,b)
    print('[%s,%s] n>=%d: %d certified, failures %s'%(sys.argv[1],sys.argv[2],n1,out['ok'],out['fail'][:3]))
    print('sample margins at n1:',[check(sys.argv[1],sys.argv[2],n1,r,pat)[2] for r in range(11)])
