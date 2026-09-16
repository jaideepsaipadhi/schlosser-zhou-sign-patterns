from q12_model import *
for ds in ('2','2.25','2.5','2.75','3'):
    d=mp.mpf(ds)
    T={h:local(h,12,d,1) for h in (5,7)}
    def amp(idx,r):
        s=0
        for h in (5,7):
            A,C=T[h][idx]
            s+=C*mp.expjpi(-2*r*mp.mpf(h)/12)*mp.expjpi(-2*d*mp.mpf(h)/12)
        return s
    main=[mp.re(amp(0,r)) for r in range(12)]
    c1=[mp.re(amp(1,r)) for r in (3,9)]
    c1all=[mp.re(amp(1,r)) for r in range(12)]
    print(ds,'A:',[mp.nstr(T[5][i][0],3) for i in range(len(T[5]))],' main amp r=0..11:',[mp.nstr(x,3) for x in main])
    print('    corr1 amp r=0..11:',[mp.nstr(x,3) for x in c1all])
