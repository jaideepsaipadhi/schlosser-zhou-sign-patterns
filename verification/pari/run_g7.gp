read("finite.gp");
e7(d)=if(d%7,1,0);
t=getabstime(); c=coeffpolys(e7,897); print("built in s ",(getabstime()-t)/1000.);
lo=487350758538673426330055236816/10^29; hi=487350758538673426343362426758/10^29;
t=getabstime(); print("G7 [2,3], n<=897: ", check(c,0,897,2,3,Vec("+--+++-"),0), "  s ",(getabstime()-t)/1000.);
t=getabstime(); print("G7 [3,delta_lo], n<=897: ", check(c,0,897,3,lo,Vec("+-++---"),0), "  s ",(getabstime()-t)/1000.);
t=getabstime(); print("G7 [delta_lo,delta_hi], n<=897 (expect only 897): ", check(c,0,897,lo,hi,Vec("+-++---"),0), "  s ",(getabstime()-t)/1000.);
p=c[898]; print("c897: roots in [lo,hi]=",polsturm(p,[lo,hi]),"  roots in [hi,4.87350758538675509]=",polsturm(p,[hi,487350758538675509/10^17]), " signs lo,hi: ",sign(subst(p,x,lo)),sign(subst(p,x,hi)));
quit
