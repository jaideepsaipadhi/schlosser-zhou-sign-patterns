# Point values of c_delta(n) (G7) for all n<=NM at delta_lo and delta_hi, certified in ball arithmetic.
import sys
from probe import *
G7=G(7); NM=int(sys.argv[1]); PREC=int(sys.argv[2])
for name,d in (('delta_lo','4.87350758538673426330055236816'),('delta_hi','4.87350758538673426343362426758')):
    f=coeffs(G7,NM,d,PREC); sg=signs(f)
    und=[n for n in range(len(sg)) if sg[n]=='?']
    bad=[n for n in range(len(sg)) if sg[n] not in '0?' and sg[n]!='+-++---'[n%7]]
    zer=[n for n in range(len(sg)) if sg[n]=='0']
    print('%s=%s  n<=%d: undecided %d, zeros %s, pattern violations %s'%(name,d,NM,len(und),zer,bad))
