# Verify the Weil transformation formulas for theta_c(tau)=sum_{J=c mod 32} exp(2 pi i tau J^2/64)
import mpmath as mp
mp.mp.dps=30
def th(c,tau,M=60):
    return mp.fsum(mp.exp(2j*mp.pi*tau*(32*m+c)**2/64) for m in range(-M,M+1))
tau=mp.mpc('0.37','0.81')
err=0
for c in range(32):
    lhs=th(c,-1/tau)
    rhs=(1j/(32*tau))**(-0.5)/32*mp.fsum(mp.exp(2j*mp.pi*c*d/32)*th(d,tau) for d in range(32))
    # principal branch check: (i/(32 tau))^(-1/2)
    err=max(err,abs(lhs-rhs)/abs(lhs))
    err=max(err,abs(th(c,tau+1)-mp.exp(2j*mp.pi*c*c/64)*th(c,tau))/abs(th(c,tau)))
print('max relative error S,T formulas:',mp.nstr(err,5))
# S1 = th6 - th10, S3 = th2 - th14 against definition S_a = sum psi_a(j) q^{j^2/16}
def psi(a,j):
    r=a-4
    if (j-r)%8: return 0
    return -1 if ((j-r)//8)%2 else 1
def S(a,tau,M=200):
    return mp.fsum(psi(a,j)*mp.exp(2j*mp.pi*tau*j*j/16) for j in range(-M,M+1) if psi(a,j))
print('S1 check', mp.nstr(abs(S(1,tau)-(th(6,tau)-th(10,tau))),5), ' S3 check', mp.nstr(abs(S(3,tau)-(th(2,tau)-th(14,tau))),5))
