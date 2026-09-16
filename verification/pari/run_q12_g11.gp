read("finite.gp");
e12(d)=my(r=d%12); if(r==1||r==11,1,if(r==5||r==7,-1,0));
c=coeffpolys(e12,780);
print("Q12 [2,3], n<=400: ", check(c,0,400,2,3,Vec("+-+--+-+-++-"),0));
bad=check(c,0,780,-1,-6441134808/10^10,Vec("+++++------+"),0);
print("Q12 [-1,delta1+], n<=780 (expect only n=22): ", bad);
\\ c_22: exactly one root in (-1,delta1+) and sign at -1 correct
p=c[23]; print("c22 roots in [-1,delta1+]: ", polsturm(p,[-1,-6441134808/10^10]), "  sign at -1: ", sign(subst(p,x,-1+10^-30)));
e11(d)=if(d%11,1,0);
c=coeffpolys(e11,310);
P=c[22]; F=factor(P)[,1]; g=0; for(i=1,#F, if(poldegree(F[i])==18, g=F[i]));
print("G11: degree-18 factor found: ", g!=0, "  its roots in [1.7584535519419,2]: ", polsturm(g,[17584535519419/10^13,2]));
sg=sign(subst(g,x,19/10));
print("G11 [gamma_lo,2], n<=310: ", check(c,0,310,17584535519419/10^13,2,Vec("+--+++-+--+"),[21,g,sg]));
quit
