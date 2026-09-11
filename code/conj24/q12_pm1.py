# Conj 24 at delta = +1 and -1: certified large-n signs (zero residues vanish identically by Andrews-Bressoud).
from flint import arb, ctx
ctx.prec=100
PI=arb.pi(); E1=arb(1).exp(); R2=arb(2).sqrt()
I1=lambda x:x.bessel_i(1)
t=(-2*PI).exp(); U=arb(0); V=arb(0)
for J in range(3,600):
    if J%24 in (2,22) and J>2: U+=t**(arb(J*J-4)/96)
    if J%24 in (10,14) and J>10: V+=t**(arb(J*J-100)/96)
F1m=(1-U)**(-1)*(1-V)**(-1)-1
def P(k,C): a=2*PI/k; return a*(2/C).sqrt()*I1(a*(2*C).sqrt())
out={}
for sign,K,amp,zero in ((1,arb('550.75'),lambda n:2*arb.cos_pi(arb(5*n)/6),(3,9)),(-1,arb('547.23'),lambda n:2*arb.cos_pi(arb(n-2)/6),(5,11))):
    fails=[]
    for n in range(20,20001):
        if n%12 in zero: continue
        C=arb(2*(n+sign)); N=int(float((4*PI*(n+1)).sqrt().upper()))+2
        Sg=(arb('0.5')+arb(N)/8)/(N+1)
        E=R2*PI*Sg*E1*(2*PI).exp()+2*R2*Sg*E1*(2*PI).exp()*F1m+2*R2*E1*K
        cnt=sum(k//2 for k in range(40,N+1,4))
        m=abs(amp(n))*P(12,C)-4*P(24,C)-6*P(36,C)-cnt*P(40,C)-E
        if not (m>0): fails.append(n)
    print('delta=%+d: nonzero residues, failures up to 20000: %d, last %s'%(sign,len(fails),fails[-3:]))
