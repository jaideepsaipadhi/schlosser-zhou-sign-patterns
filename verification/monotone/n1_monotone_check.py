# Verifies the hypothesis behind "one check at n1 covers all n>=n1" (Lemma on envelopes):
# every competitor/dominant ratio is bounded by C * y^a * exp(-g*y) with y = sqrt(C_n) ~ sqrt(2n),
# g = alpha_dom - alpha_comp > 0, and this bound is decreasing for y >= a/g.  We use a = 5 for all terms
# (the largest power arising: N^2 counting factor ~ y^2, derivative factors D_k ~ y, error-derivative ~ N^3 log N).
import math
A_POW=5
def alpha(k,B): return 2*math.pi/k*math.sqrt(B)
cases=[]
def g7(name,a,b,n1,dom_k=7,dom_B=None,corr_dom=False):
    # competitors: k=7 correction (if b>4, unless it is the dominant term), k=14,21,28 main and correction
    Bd = (a/2-2) if corr_dom else a/2
    ad = alpha(dom_k,Bd)
    comps=[]
    if b>4 and not corr_dom: comps.append(('k7 corr',alpha(7,b/2-2)))
    for k in (14,21,28):
        comps.append(('k%d main'%k,alpha(k,b/2)))
        if b>4: comps.append(('k%d corr'%k,alpha(k,b/2-2)))
    comps.append(('error (const or poly)',0.0))
    y1=math.sqrt(2*n1-b/2)
    g=min(ad-ac for _,ac in comps); worst=min(comps,key=lambda c:ad-c[1])
    cases.append((name,n1,y1,g,A_POW/g if g>0 else float('inf'),worst[0]))
g7('G7 dominance [2,3)',2,3,898); g7('G7 dominance (3,4.87]',3,4.87,898); g7('G7 dominance [4.87,4.87350758538]',4.87,4.87350758538,898)
g7('G7 window, residues !=1',4.87350758538,4.87350758539,898); g7('G7 window residue 1 (n>=2500)',4.87350758538,4.87350758539,2500)
g7('G7 derivative (window)',4.87350758538,4.87350758539,898); g7('G7 eps near 3',2.999999999,3.000000001,898)
g7('G7 delta=1',1,1,400); g7('G7 delta=3',3,3,400); g7('G7 delta=5 (main-dominated residues)',5,5,400)
g7('G7 delta=5 residues 1,2,6 (correction dominant; k7,14,21 main vanish)',5,5,400,corr_dom=True)
# G11: dominant k=11 main, B=5d/6; competitors k=22,33,44 main (+corr if d>2.4)
def g11(name,a,b,n1):
    ad=alpha(11,5*a/6); comps=[(k,alpha(k,5*b/6)) for k in (22,33,44)]
    if b>2.4: comps+= [('k11 corr',alpha(11,5*b/6-2))]+[('k%d corr'%k,alpha(k,5*b/6-2)) for k in (22,33,44)]
    y1=math.sqrt(2*n1-5*b/6); g=min(ad-c[1] for c in comps)
    cases.append((name,n1,y1,g,A_POW/g,'-'))
g11('G11 dominance [gamma,2]',1.7584,2,300); g11('G11 delta=1',1,1,600); g11('G11 delta=3',3,3,600)
# special: G7 delta=5 corr-dominant: k=14,21 main amplitudes vanish exactly, so exclude them
ad=alpha(7,0.5); comps=[alpha(14,0.5),alpha(21,0.5),alpha(28,2.5)]
cases.append(('G7 delta=5 corr-dominant (exact zeros used)',400,math.sqrt(800-2.5),min(ad-c for c in comps),A_POW/min(ad-c for c in comps),'k28 main'))
ok=True
print('%-62s %6s %7s %8s %8s  %s'%('check','n1','y(n1)','gap g','a/g','status'))
for name,n1,y1,g,ag,w in cases:
    if 'residues 1,2,6 (correction dominant' in name: continue
    st='OK' if (g>0 and y1>=ag) else 'FAIL'; ok&=(st=='OK')
    print('%-62s %6d %7.2f %8.4f %8.2f  %s'%(name,n1,y1,g,ag,st))
print('ALL HYPOTHESES SATISFIED:',ok)
