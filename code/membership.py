"""Random-point membership test of fibre completeness.
Draw a random point Y' of V (random complex state, random branches of hb,hc), track the single path
from F6(Y') to c0 (via a random intermediate target) and test whether the endpoint is in the known set S.
If S misses a fraction f of the sheets, misses occur with frequency ~f.  Misses are new sheets."""
import numpy as np, sys, json, time
from move_fibre import *
from track_robust import track_many
from scipy.spatial import cKDTree
inp=sys.argv[1]; T=int(sys.argv[2]); seed=int(sys.argv[3]); outp=sys.argv[4]
d=np.load(inp); S=d['S']; c0=d['c0']
tree=cKDTree(np.c_[S.real,S.imag]/(1+np.abs(S).max(1))[:,None])
rng=np.random.default_rng(seed); rc=lambda *n:(rng.normal(size=n)+1j*rng.normal(size=n))/np.sqrt(2)
Ys=[];starts=[];mids=[]
for t in range(T):
    sc=rng.choice([0.3,1.0,3.0]); q=rc(7)*sc+np.array([1.2,0.4,0.9,0,0,0,0])*rng.integers(0,2)
    l,x,y=q[:3]; hb=rng.choice([1,-1])/np.sqrt(x*x+y*y); hc=rng.choice([1,-1])/np.sqrt((x-l)**2+y*y)
    Y=np.r_[q,1/l,hb,hc]; Ys.append(Y); c=F6(Y); starts.append(c); mids.append(0.5*(c+c0)*(1+rc(7)))
Ys=np.array(Ys); res=[]; t0=time.time()
for i in range(T):
    y=Ys[i:i+1].copy(); ok=True
    for a,b in ((starts[i],mids[i]),(mids[i],c0)):
        y,o=track_many(y,a,b,E,C,off,md); ok&=bool(o[0])
    if ok:
        dist=tree.query(y[0].real.tolist()+y[0].imag.tolist() if False else np.r_[y[0].real,y[0].imag]/(1+np.abs(y[0]).max()))[0]
        res.append(dict(ok=True,hit=bool(dist<1e-7),dist=float(dist),end=[str(v) for v in y[0]]))
    else: res.append(dict(ok=False))
hits=sum(r.get('hit',False) for r in res); oks=sum(r['ok'] for r in res)
print(f'tested {T} tracked {oks} hits {hits} misses {oks-hits} time {time.time()-t0:.0f}s')
json.dump(dict(set_size=len(S),tested=T,tracked=oks,hits=hits,misses=oks-hits,results=res),open(outp,'w'))
