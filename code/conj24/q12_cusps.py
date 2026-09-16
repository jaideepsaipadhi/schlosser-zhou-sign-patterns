# Exact typing of all cusps h/k (k<=KMAX) for Q12 via Poisson/Gauss sums, zero tests mod Phi_{24k}.
import sys
from math import gcd
from flint import fmpz_poly
def psi(a,j):
    r=a-6
    if (j-r)%12: return 0
    return -1 if ((j-r)//12)%2 else 1
def Gpair(a,h,k,l,Phi):
    M=24*k; c=[0]*M
    for s in range(M):
        p=psi(a,s)
        if p:
            c[(h*s*s+l*s)%M]+=p
            if l: c[(h*s*s-l*s)%M]+=p
    return fmpz_poly(c)%Phi
KMAX=int(sys.argv[1]); res={}
for k in range(1,KMAX+1):
    Phi=fmpz_poly.cyclotomic(24*k)
    for h in range(k):
        if gcd(h,k)!=1: continue
        J={}
        for a in (1,5):
            l=0
            while Gpair(a,h,k,l,Phi)==0: l+=1
            J[a]=l
        res[(h,k)]=(J[1],J[5])
gp=sorted(c for c,v in res.items() if v[0]<v[1]); gm=sorted(c for c,v in res.items() if v[0]>v[1])
print('growing for delta>0:',gp)
print('growing for delta<0:',gm)
print('(+) <=> 12|k, h=5,7 mod 12:',all(((k%12==0) and h%12 in (5,7))==((h,k) in gp) for (h,k) in res))
print('(-) <=> 12|k, h=1,11 mod 12:',all(((k%12==0) and h%12 in (1,11))==((h,k) in gm) for (h,k) in res))
print('types:',sorted(set(res.values())))
