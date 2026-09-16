import sys, time
from g7_dom import check
import mpmath as mp
mp.mp.dps=40
def tile(a,b,n1,pat,minw,out,depth=0):
    ok=all(check(mp.nstr(a,40),mp.nstr(b,40),n1,r,pat)[0] for r in RES) 
    # float() loses precision for tiny widths: use string-precise version
    if ok: out['ok']+=1; out.setdefault('leaves',[]).append((a,b)); return
    if b-a<minw: out['fail'].append((a,b)); return
    m=(a+b)/2; tile(a,m,n1,pat,minw,out,depth+1); tile(m,b,n1,pat,minw,out,depth+1)
a,b=mp.mpf(sys.argv[1]),mp.mpf(sys.argv[2]); n1=int(sys.argv[3]); pat=sys.argv[4]; minw=mp.mpf(sys.argv[5]); RES=[int(c) for c in sys.argv[6]] if len(sys.argv)>6 else range(7)
out={'ok':0,'fail':[]}; t0=time.time()
tile(a,b,n1,pat,minw,out)
print('leaves:',[(mp.nstr(x,12),mp.nstr(y,12)) for x,y in out.get('leaves',[])][:8]); print('[%s,%s] n>=%d: %d intervals certified, %d unresolved (<%s): %s  %.0fs'%(mp.nstr(a,12),mp.nstr(b,12),n1,out['ok'],len(out['fail']),mp.nstr(minw,3),[(mp.nstr(x,15),mp.nstr(y,15)) for x,y in out['fail'][:4]],time.time()-t0))
