from q12_model import *
import sys
d=mp.mpf(sys.argv[1]); sign=1 if d>0 else -1
cusps={}
for k in (12,24,36):
    for h in range(k):
        if math.gcd(h,k)!=1: continue
        if (sign>0 and h%12 in (5,7)) or (sign<0 and h%12 in (1,11)):
            cusps[(h,k)]=local(h,k,d,sign)
for hk,t in list(cusps.items())[:2]: print(hk,[(mp.nstr(A,5),mp.nstr(C,5)) for A,C in t])
N=500
e=1
bb=b
f=[mp.mpf(1)]+[0]*N
for n in range(1,N+1): f[n]=d*mp.fsum(bb[k]*f[n-k] for k in range(1,n+1))/n
for n in range(N-12,N+1):
    p=predict(n,d,cusps,sign)
    print(n,n%12,mp.nstr(f[n],10),mp.nstr(mp.re(p),10),mp.nstr(mp.im(p),3))
