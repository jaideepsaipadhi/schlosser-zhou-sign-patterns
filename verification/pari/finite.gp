\\ Independent re-implementation (PARI/GP). d_n = n! c_n in Z[x]:  d_n = x * sum_k b_k * (n-1)!/(n-k)! * d_{n-k}.
coeffpolys(epsf, N) = {
  my(b=vector(N), d=vector(N+1));
  for(k=1,N, b[k] = -sumdiv(k,t, epsf(t)*t));
  d[1]=1;
  for(n=1,N, my(s=0, fac=1); for(k=1,n, s += b[k]*fac*d[n-k+1]; fac*= (n-k)); d[n+1]=x*s);
  d
};
\\ sign of c_n = sign of d_n (n!>0).  Checks as described.
check(c, N1, N2, a, bb, pat, special) = {
  my(L=#pat, bad=List(), p, w, f, s, nr);
  for(n=N1,N2,
    p=c[n+1]; w=if(pat[n%L+1]=="+",1,-1); f=1;
    if(special && n==special[1], p=p/special[2]; f*=special[3]);
    while(subst(p,x,a)==0 && p!=0, p=p/(x-a));
    while(subst(p,x,bb)==0 && p!=0, p=p/(x-bb); f=-f);
    nr=polsturm(p,[a,bb]);
    s=sign(subst(p,x,a))*f;
    if(nr!=0 || s!=w, listput(bad,[n,nr,s,w])));
  bad
};
