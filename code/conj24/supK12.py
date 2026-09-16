# Certified sup |R^{sgn}| on Im tau' >= 1 for the non-growing Q12 cusp types (sgn=+1: R=S1/S5, sgn=-1: inverse).
import pickle, heapq, sys, time
from flint import arb, acb, ctx
from orbit12 import *
ctx.prec=64
SGN=int(sys.argv[1]) if len(sys.argv)>1 else 1
orbit,edges=pickle.load(open('orbit12.pkl','rb'))
Tn={i:j for i,n,j in edges if n=='T'}
seen=set(); cycles=[]
for i in range(len(orbit)):
    if i in seen: continue
    c=[i]; seen.add(i); j=Tn[i]
    while j!=i: c.append(j); seen.add(j); j=Tn[j]
    cycles.append(c)
PI=arb.pi()
def zval(a):
    re=arb(0); im=arb(0)
    for i,x in enumerate(a):
        x=int(x)
        if x: re+=x*arb.cos_pi(arb(2*i)/192); im+=x*arb.sin_pi(arb(2*i)/192)
    return acb(re,im)
JMAX=60
def coeffs(w):
    out=[]
    for J in range(JMAX+1):
        c=w[J%48]
        out.append(None if not any(c) else (zval(c) if J==0 else 2*zval(c)))
    return out
def Fval(A,tau):
    s=acb(0)
    for J,a in enumerate(A):
        if a is not None: s+=a*(2*PI*acb(0,1)*tau*J*J/96).exp()
    return s
def tailJ(A,y):
    amax=max(abs(a) for a in A if a is not None); t=(-2*PI*y/96).exp()
    return amax*2*t**((JMAX+1)**2)
Y=arb(6); res=[]
for cyc in cycles:
    w1,w5=orbit[cyc[0]]
    if SGN<0: w1,w5=w5,w1
    J1=min(c for c in range(25) if any(w1[c])); J5=min(c for c in range(25) if any(w5[c]))
    if J1<J5: continue
    L=len(cyc); A1=coeffs(w1); A5=coeffs(w5)
    t=(-2*PI*Y).exp()
    def eps(A,J0):
        s=arb(0)
        for J,a in enumerate(A):
            if a is not None and J>J0: s+=abs(a)/abs(A[J0])*t**(arb(J*J-J0*J0)/96)
        return s+tailJ(A,Y)/abs(A[J0])
    e1,e5=eps(A1,J1),eps(A5,J5)
    tailb=abs(A1[J1])/abs(A5[J5])*t**(arb(J1*J1-J5*J5)/96)*(1+e1)/(1-e5)
    def box_up(x0,y0,hx,hy):
        tau=acb(arb(x0+hx/2,hx/2),arb(y0+hy/2,hy/2))
        num=Fval(A1,tau); den=Fval(A5,tau); ad=abs(den)-tailJ(A5,arb(y0))
        if not (ad>0): return float('inf')
        return float(((abs(num)+tailJ(A1,arb(y0)))/ad).upper())
    def pt(x,y):
        tau=acb(x,y); return float(abs(Fval(A1,tau)/Fval(A5,tau)).lower())
    best=0.0
    for i in range(4*L+1):
        for j in range(11): best=max(best,pt(L*i/(4*L),1+5*j/10))
    heap=[]; nx,ny=8*L,20
    for i in range(nx):
        for j in range(ny):
            hx=L/nx; hy=5/ny; heapq.heappush(heap,(-box_up(i*hx,1+j*hy,hx,hy),i*hx,1+j*hy,hx,hy))
    while -heap[0][0]>best*1.02 and len(heap)<300000:
        u,x0,y0,hx,hy=heapq.heappop(heap); best=max(best,pt(x0+hx/2,y0+hy/2))
        for dx in (0,1):
            for dy in (0,1):
                heapq.heappush(heap,(-box_up(x0+dx*hx/2,y0+dy*hy/2,hx/2,hy/2),x0+dx*hx/2,y0+dy*hy/2,hx/2,hy/2))
    K=max(-heap[0][0],float(tailb.upper())); res.append(K)
    print('type (%d,%d) L=%2d lower=%.5f cert_upper=%.5f tail=%.5f boxes=%d'%(J1,J5,L,best,-heap[0][0],float(tailb.upper()),len(heap)),flush=True)
print('K12 (sgn=%+d) certified upper bound = %.6f'%(SGN,max(res)))
