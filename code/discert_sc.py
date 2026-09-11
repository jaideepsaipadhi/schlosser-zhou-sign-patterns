from flint import arb, ctx
ctx.prec=200
from consts import Econst, PI
d=arb('2.66448'); I1=lambda x:x.bessel_i(1)
def P(k,B,y): return (2*PI/k)*B.sqrt()/y*I1((2*PI/k)*B.sqrt()*y)
for n in (1000010,):
    y=(2*n+d).sqrt(); N=int(float((4*PI*n+2*PI*d).sqrt().upper()))+1
    main16=4*arb.cos_pi(arb(n)/8)*arb.cos_pi((n+d)/2)*P(16,d,y)
    C8=-2*d*arb.cos_pi(arb(3*n-1)/4)*P(8,d-2,y)
    rest=4*d*P(16,d-2,y)+sum((arb(k)/2)*(P(k,d,y)+d*P(k,d-2,y)) for k in range(24,N+1,4))+Econst(d,N)
    ub=main16+C8+rest
    print(n,'upper bound on c_delta(n):',ub.str(5),' negative:',ub<0,' relative:',(ub/abs(main16)).str(4))
