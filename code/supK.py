# Step 2: certified sup of |R| on Im tau' >= 1 for the 9 bounded cusp types.
import pickle, heapq, sys, time
import numpy as np
from flint import arb, acb, ctx
ctx.prec=64
exec(open('orbit.py').read().split('one=np.array')[0])
exec(open('cusp_types.py').read().split("print('number")[0].split("orbit,edges=pickle.load")[0])
orbit,edges=pickle.load(open('orbit.pkl','rb'))
Tnext={i:j for i,n,j in edges if n=='T'}
seen=set(); cycles=[]
for i in range(len(orbit)):
    if i in seen: continue
    cyc=[i]; seen.add(i); j=Tnext[i]
    while j!=i: cyc.append(j); seen.add(j); j=Tnext[j]
    cycles.append(cyc)
PI=arb.pi()
def zval(a):
    re=arb(0); im=arb(0)
    for i,x in enumerate(a):
        x=int(x)
        if x:
            t=arb(2*i)/128; re+=x*arb.cos_pi(t); im+=x*arb.sin_pi(t)
    return acb(re,im)
JMAX=40
def coeffs(w):   # A_J for J=0..JMAX  (A_0=w_0, A_J=2 w_{J mod 32})
    out=[]
    for J in range(JMAX+1):
        c=w[J%32]
        out.append(None if not any(c) else (zval(c) if J==0 else 2*zval(c)))
    return out
def Fval(A,tau):
    s=acb(0)
    for J,a in enumerate(A):
        if a is not None: s+=a*(2*PI*acb(0,1)*tau*J*J/64).exp()
    return s
def tailJ(A,y):   # |sum_{J>JMAX}| <= Amax * sum_{J>JMAX} e^{-2pi y J^2/64}
    amax=max(abs(a) for a in A if a is not None)
    t=(-2*PI*y/64).exp()
    return amax*2*t**((JMAX+1)**2)          # geometric domination, t^{J^2} decays super-geometrically
Y=arb(6)
results=[]
t0=time.time()
for cyc in cycles:
    w1,w3=orbit[cyc[0]]
    J1=min(c for c in range(17) if any(w1[c])); J3=min(c for c in range(17) if any(w3[c]))
    if J1<J3: continue
    L=len(cyc)
    A1=coeffs(w1); A3=coeffs(w3)
    # ---- tail y>=Y: |R| <= |A1_J1/A3_J3| t^{(J1^2-J3^2)/64} (1+e1)/(1-e3)
    t=(-2*PI*Y).exp()
    def eps(A,J0):
        s=arb(0)
        for J,a in enumerate(A):
            if a is not None and J>J0: s+=abs(a)/abs(A[J0])*t**(arb(J*J-J0*J0)/64)
        return s+tailJ(A,Y)/abs(A[J0])
    e1,e3=eps(A1,J1),eps(A3,J3)
    tailbound=abs(A1[J1])/abs(A3[J3])*t**(arb(J1*J1-J3*J3)/64)*(1+e1)/(1-e3)
    # ---- box [0,L] x [1,Y]: branch and bound
    def box_up(x0,y0,hx,hy):
        tau=acb(arb(x0+hx/2,hx/2),arb(y0+hy/2,hy/2))
        tb=tailJ(A1,arb(y0)); tb3=tailJ(A3,arb(y0))
        num=Fval(A1,tau); den=Fval(A3,tau)
        ad=abs(den)-tb3
        if not (ad>0): return float('inf')
        return float(((abs(num)+tb)/ad).upper())
    def pt(x,y):
        tau=acb(x,y); return float(abs(Fval(A1,tau)/Fval(A3,tau)).lower())
    best=0.0
    for i in range(33):
        for j in range(11):
            best=max(best,pt(L*i/32,1+5*j/10))
    heap=[]
    nx,ny=8*L,20
    for i in range(nx):
        for j in range(ny):
            hx=L/nx; hy=5/ny
            u=box_up(i*hx,1+j*hy,hx,hy)
            heapq.heappush(heap,(-u,i*hx,1+j*hy,hx,hy))
    target=best*1.02
    while -heap[0][0]>target and len(heap)<200000:
        u,x0,y0,hx,hy=heapq.heappop(heap)
        best=max(best,pt(x0+hx/2,y0+hy/2))
        target=best*1.02
        for dx in (0,1):
            for dy in (0,1):
                uu=box_up(x0+dx*hx/2,y0+dy*hy/2,hx/2,hy/2)
                heapq.heappush(heap,(-uu,x0+dx*hx/2,y0+dy*hy/2,hx/2,hy/2))
    boxsup=-heap[0][0]
    K=max(boxsup,float(tailbound.upper()))
    results.append(((J1,J3),L,best,boxsup,float(tailbound.upper()),K,len(heap)))
    print('type (%d,%d) L=%d  lower=%.6f  box upper=%.6f  tail(y>=6) upper=%.6f  => K=%.6f  boxes=%d  %.0fs'%(J1,J3,L,best,boxsup,float(tailbound.upper()),K,len(heap),time.time()-t0),flush=True)
KNG=max(r[5] for r in results)
print('K_NG (certified upper bound) =',KNG)
pickle.dump(results,open('supK.pkl','wb'))
