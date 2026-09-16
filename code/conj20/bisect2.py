# Refine delta_c: signs of c_n(delta) for n<=NM, n=1 mod 7, certified in ball arithmetic.
import sys
from probe import *
import mpmath as mp
mp.mp.dps=50
G7=G(7); NM=int(sys.argv[1]); PREC=int(sys.argv[2])
lo=mp.mpf(sys.argv[3]); hi=mp.mpf(sys.argv[4]); steps=int(sys.argv[5])
def run(d):
    f=coeffs(G7,NM,mp.nstr(d,45),PREC); sg=signs(f)
    und=[n for n in range(len(sg)) if sg[n]=='?']
    bad=[n for n in range(len(sg)) if n%7==1 and sg[n]=='+']
    other=[n for n in range(len(sg)) if sg[n] not in '0?' and n%7!=1 and sg[n]!='+-++---'[n%7]]
    return bad,und,other,f
for it in range(steps):
    mid=(lo+hi)/2; bad,und,oth,f=run(mid)
    assert not und, ('undecided',und[:5])
    if bad: hi=mid; fb=bad[:4]
    else: lo=mid
    print(it,mp.nstr(lo,30),mp.nstr(hi-lo,3),'viol at mid:',bad[:4],'other-residue violations:',oth[:3],flush=True)
b1,u1,o1,f1=run(lo); b2,u2,o2,f2=run(hi)
print('FINAL lo=%s: violations %s ; hi=%s: violations %s'%(mp.nstr(lo,30),b1[:5],mp.nstr(hi,30),b2[:5]))
print('c_883 at lo:',f1[883].str(5),' at hi:',f2[883].str(5))
