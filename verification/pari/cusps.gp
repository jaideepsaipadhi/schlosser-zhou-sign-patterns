chi(a,j,m,off)={my(r=a-off); if((j-r)%m, 0, if(((j-r)/m)%2, -1, 1))};
Gs(a,h,k,l,M,m,off)={my(P=polcyclo(M,'y), s=0); for(t=0,M-1, my(v=chi(a,t,m,off)); if(v, s+=v*('y^((h*t^2+l*t)%M)+'y^((h*t^2-l*t)%M)))); lift(Mod(s,P))};
Jl(a,h,k,M,m,off)={my(l=0); while(Gs(a,h,k,l,M,m,off)==0, l++); l};
typeQ8(K)={my(bad=0,cnt=0); for(k=1,K, for(h=0,k-1, if(gcd(h,k)==1, cnt++;
   my(M=16*k, j1=Jl(1,h,k,M,8,4), j3=Jl(3,h,k,M,8,4), g1=(j1<j3), p1=(k%8==0 && (h%8==3||h%8==5)));
   if(g1!=p1 || (g1 && [j1,j3]!=[2,6]), bad++)))); [cnt,bad]};
typeQ12(K)={my(bad=0,cnt=0); for(k=1,K, for(h=0,k-1, if(gcd(h,k)==1, cnt++;
   my(M=24*k, j1=Jl(1,h,k,M,12,6), j5=Jl(5,h,k,M,12,6), A1=(j1<j5), A2=(j1>j5), B1=(k%12==0 && (h%12==5||h%12==7)), B2=(k%12==0 && (h%12==1||h%12==11)));
   if(A1!=B1 || A2!=B2, bad++)))); [cnt,bad]};
print("Q8 cusps k<=40 [checked, mismatches]: ", typeQ8(40));
print("Q12 cusps k<=24 [checked, mismatches]: ", typeQ12(24));
quit
