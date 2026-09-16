# Exact leading exponents J1,J3 at every cusp h/k, k<=KMAX, from Lemma 2.1 Gauss sums,
# with exact zero tests in Z[zeta_M] (reduction modulo the cyclotomic polynomial Phi_M).
import sys
from math import gcd
from flint import fmpz_poly
KMAX=int(sys.argv[1])
def psi(a,j):
    r=a-4
    if (j-r)%8: return 0
    return -1 if ((j-r)//8)%2 else 1
res={}
for k in range(1,KMAX+1):
    M=16*k; Phi=fmpz_poly.cyclotomic(M)
    for h in range(k):
        if gcd(h,k)!=1: continue
        J={}
        for a in (1,3):
            def Gpair(l):
                c=[0]*M
                for s in range(M):
                    p=psi(a,s)
                    if p:
                        c[(h*s*s+l*s)%M]+=p
                        if l: c[(h*s*s-l*s)%M]+=p
                return fmpz_poly(c)%Phi
            l=0
            while Gpair(l)==0: l+=1
            J[a]=l
        res[(h,k)]=(J[1],J[3])
grow=sorted((h,k) for (h,k),(j1,j3) in res.items() if j1>j3)
print('KMAX=%d: cusps checked %d'%(KMAX,len(res)))
print('growing cusps:',grow)
ok=all((k%8==0 and h%8 in (1,7))==(res[(h,k)][0]>res[(h,k)][1]) for (h,k) in res)
print('inverse-growing <=> (8|k and h = 1,7 mod 8):',ok)
print('all of type (6,2):',all(res[c]==(6,2) for c in grow))
print('distinct (J1,J3) values:',sorted(set(res.values())))
