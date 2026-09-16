from flint import arb, ctx
ctx.prec=100
exec(open('q12_neg.py').read().split("S_LO=arb(sys.argv[3])")[0])
S_LO=arb('0.6441134808'); n=arb(20000)
yl=(2*(n-1)).sqrt(); yh=(2*n).sqrt(); N=(2*PI).sqrt()*yh/ (arb(2).sqrt()) *arb(2).sqrt()+3   # N<=sqrt(4 pi n)+2
N=(4*PI*n).sqrt()+3
Uenv=lambda x:x.exp()/(2*PI*x).sqrt()
def low(k,s): x=(2*PI/k)*(2*s).sqrt()*yl; kap=(2*PI*x).sqrt()*(-x).exp()*x.bessel_i(1); return (2*PI/k)*(2*s).sqrt()/yh*kap*Uenv(x)
def up(k,s): return (2*PI/k)*(2*s).sqrt()/yl*Uenv((2*PI/k)*(2*s).sqrt()*yh)
cnt=(N/4)*(N/2)
D=(2*PI/3)*(1-(PI*(1-S_LO)/3)**2/6)*low(12,S_LO)
rest=(2*PI/3)*up(24,arb(1))+(6+cnt)*(N*N/2+1+1+(2*PI/36)*(n/S_LO).sqrt()*2)*up(36,arb(1))
Gs=arb('0.5')+N/8
arcs=Gs*2*E1*(2*PI).exp()*((N*N/2+1)*(PI/R2)+(2+PI/2)/S_LO+1)
chord=Gs*2*R2*E1*(2*PI).exp()*((N*N/2+4*PI*N*N+1)*F1m+F1*ell)
ng=2*R2*E1*KP*(LOGK+1+arb('0.556')*N*N)
print('degenerate ratio sum at n=20000:',((rest+arcs+chord+ng)/D).str(3))
cmin=arb('0.1478')      # min over nondegenerate residues and s in [S_LO,1] of |cos(pi(n-2s)/6)|
E=R2*PI*Gs*E1*(2*PI).exp()+2*R2*Gs*E1*(2*PI).exp()*F1m+2*R2*E1*KP
D2=2*cmin*low(12,S_LO)
print('nondegenerate ratio sum:',((4*up(24,arb(1))+(6+cnt)*up(36,arb(1))+E)/D2).str(3))
print('gaps: main vs k24 =',((2*PI/12)*(2*S_LO).sqrt()-(2*PI/24)*arb(2).sqrt()).str(3),'(>0; all other ratios are poly(y)e^{-gy})')
import math
print('check cmin:',min(abs(math.cos(math.pi*(r-2*(0.6441134808+i*(1-0.6441134808)/20000))/6)) for r in (0,1,2,3,4,6,7,8,9,10) for i in range(20001)))
