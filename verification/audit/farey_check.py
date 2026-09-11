# Numerical audit of the Farey/Ford chord lemma (Lemma 4.2 of the joint paper) for many N.
from fractions import Fraction as F
import math, cmath
worst={'len':0,'re':0,'arc':0,'re1z':9}
for N in list(range(1,60))+[97,150,211]:
    fr=sorted(set(F(h,k) for k in range(1,N+1) for h in range(0,k) if math.gcd(h,k)==1))
    L=len(fr)
    for i,q in enumerate(fr):
        h,k=q.numerator,q.denominator
        p=fr[i-1] if i>0 else fr[-1]-1; r=fr[i+1] if i+1<L else fr[0]+1
        k1,k2=p.denominator,r.denominator
        z1=complex(k*k,k*k1)/(k*k+k1*k1); z2=complex(k*k,-k*k2)/(k*k+k2*k2)
        worst['len']=max(worst['len'],abs(z1-z2)/(2*math.sqrt(2)*k/(N+1)))
        worst['re']=max(worst['re'],max(z1.real,z2.real)/(2*k*k/(N+1)**2))
        for t in [j/50 for j in range(51)]:
            z=z1+(z2-z1)*t; worst['re1z']=min(worst['re1z'],(1/z).real)
        for z in (z1,z2):   # small arc from 0 to z on |z-1/2|=1/2
            ang=abs(cmath.phase((z-0.5)/(-0.5)))
            worst['arc']=max(worst['arc'],0.5*ang/((math.pi/2)*math.sqrt(2)*k/(N+1)))
print('max ratios (must be <=1): chord length %.4f, Re z %.4f, arc length %.4f ; min Re(1/z) on chords %.6f (must be >=1)'%(worst['len'],worst['re'],worst['arc'],worst['re1z']))
