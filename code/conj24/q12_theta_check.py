# Check Q12 = q^{-1} S1/S5 with S1 = th10 - th14, S5 = th2 - th22, th_c = sum_{J = c mod 48} q^{J^2/96}
import mpmath as mp
mp.mp.dps=30
tau=mp.mpc('0.31','0.27'); q=mp.exp(2j*mp.pi*tau)
def th(c,M=40): return mp.fsum(mp.exp(2j*mp.pi*tau*(48*m+c)**2/96) for m in range(-M,M+1))
Q=mp.fprod((1-q**(12*j+1))*(1-q**(12*j+11))/((1-q**(12*j+5))*(1-q**(12*j+7))) for j in range(200))
print(abs(Q - (th(10)-th(14))/(th(2)-th(22))/q))
# S-transform: th_c(-1/tau) = (i/(48 tau))^{-1/2}/48 * sum_d e(cd/48) th_d(tau)
c=10
lhs=mp.fsum(mp.exp(2j*mp.pi*(-1/tau)*(48*m+c)**2/96) for m in range(-40,41))
rhs=(1j/(48*tau))**(-0.5)/48*mp.fsum(mp.exp(2j*mp.pi*c*d/48)*th(d) for d in range(48))
print(abs(lhs-rhs))
