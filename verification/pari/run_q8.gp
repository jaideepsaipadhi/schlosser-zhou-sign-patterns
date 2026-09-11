read("finite.gp");
e8(d)=my(r=d%8); if(r==3||r==5,-1,if(r==1||r==7,1,0));   \\ Q8: (1-q^d)^{+1} for d=1,7; ^{-1} for 3,5
c=coeffpolys(e8,700);
pat=Vec("+-++-+--");
print("Q8 [8/3,4], n<=700 problems: ", check(c,0,700,8/3,4,pat,0));
pn=Vec("+++----+");
q=x^2-7*x-6; sg=sign(subst(q,x,-9/10));
print("Q8 [-1,-772/1000] (c_4 / (x^2-7x-6)), n<=400 problems: ", check(c,0,400,-1,-772/1000,pn,[4,q,sg]));
print("Q8 [beta+,8/3], n<=700 problems: ", check(c,0,700,2664479110226973/10^15,8/3,pat,0));
quit
