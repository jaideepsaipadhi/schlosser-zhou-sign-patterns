# Main-term amplitude of G7^delta at k=7 from Schlosser-Zhou Thm 3: A_7(n)=(1/7) sum_h e^{-pi i delta s(h,7)} e(-hn/7)
from fractions import Fraction as F
import math, cmath
def s(h,k):
    tot=F(0)
    for j in range(1,k):
        x=F(j,k); y=F(j*h,k)
        fx=x-math.floor(x)-F(1,2) if x.denominator!=1 else 0
        fy=y-math.floor(y)-F(1,2) if y.denominator!=1 else 0
        tot+=fx*fy
    return tot
print('s(h,7):',[str(s(h,7)) for h in range(1,7)])
def A(n,d): return sum(cmath.exp(-1j*math.pi*d*float(s(h,7)))*cmath.exp(-2j*math.pi*h*n/7) for h in range(1,7))/7
for d in (2,2.5,3,3.5,4,4.5,5):
    print(d,[round(A(r,d).real,4) for r in range(7)], ' imag max',max(abs(A(r,d).imag) for r in range(7)))
def A14(n,d):
    return sum(cmath.exp(-1j*math.pi*d*float(s(h,14)-s(h,2)))*cmath.exp(-2j*math.pi*h*n/14) for h in range(14) if math.gcd(h,14)==1)/14
for d in (4.5,4.9,5):
    print('A14 d=%s:'%d,[ (r,round(A14(r,d).real,4)) for r in (1,8,2,9,6,13)])
