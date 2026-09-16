from flint import arb, ctx
ctx.prec=100
PI=arb.pi(); E1=arb(1).exp(); R2=arb(2).sqrt(); K=arb('23.6475')
t=(-2*PI).exp(); U=arb(0); V=arb(0)
for J in range(3,400):
    if J%16 in (2,14) and J>2: U+=t**(arb(J*J-4)/64)
    if J%16 in (6,10) and J>6: V+=t**(arb(J*J-36)/64)
F2m=(1-U)**(-2)*(1-V)**(-2)-1
yl=arb('199.9'); y=arb('200.1'); N=(2*PI).sqrt()*y+3
x0=(PI/8)*R2*yl
kap=(2*PI*x0).sqrt()*(-x0).exp()*x0.bessel_i(1)
Uenv=lambda x:x.exp()/(2*PI*x).sqrt()
low=lambda k:(2*PI/k)*R2/y*kap*Uenv((2*PI/k)*R2*yl)
up=lambda k:(2*PI/k)*R2/yl*Uenv((2*PI/k)*R2*y)
Sg=(arb('0.5')+N/8)/(N+1)
E=R2*PI*Sg*E1*(2*PI).exp()+2*R2*Sg*E1*(2*PI).exp()*F2m+2*R2*E1*K**2
T24=(N/4)*(N/2)*up(24)
r1=(4*up(16)+T24+E)/(R2*low(8)); r2=(T24+E)/(2*R2*low(16))
print('delta=+-2 tail at y~200: k=8-dominated residues ratio',r1.str(3),' k=16-dominated ratio',r2.str(3), ' both <1:', r1<1 and r2<1)
