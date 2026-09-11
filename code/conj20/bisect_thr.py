import sys
from probe import *
from flint import arb
import mpmath as mp
mp.mp.dps=40
G7=G(7); NM=int(sys.argv[1])
def viol(dstr):
    f=coeffs(G7,NM,dstr,2600); sg=signs(f)
    und=[n for n in range(len(sg)) if sg[n]=='?']
    bad=[n for n in range(len(sg)) if n%7==1 and sg[n]=='+']
    return bad,und
lo=mp.mpf(sys.argv[3]) if len(sys.argv)>3 else mp.mpf('4.873507'); hi=mp.mpf('4.873507585386755088902325')
for it in range(int(sys.argv[2])):
    mid=(lo+hi)/2; bad,und=viol(mp.nstr(mid,35))
    assert not und
    if bad: hi=mid; last=bad[:3]
    else: lo=mid
    print(it, mp.nstr(lo,25), mp.nstr(hi-lo,3), 'first violations at hi:', last if bad else '', flush=True)
