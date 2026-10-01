"""PHASE 9: transport the whole fibre to real physical targets in diverse dynamical regimes.
For each target: endpoints, failures, real physical sheets, membership of the true source state,
min relative separation of A7 over ALL endpoints, and real-physical A7 separation."""
import numpy as np, sys, json
from move_fibre import *
from track_robust import track_many
from scipy.spatial import cKDTree
d=np.load(sys.argv[1]); S=d['S']; c0=d['c0']; m=np.array([1.,2.,3.])
def chart(r,v):
    a=r[1]-r[0]; b=r[2]-r[0]; va=v[1]-v[0]; vb=v[2]-v[0]; l=np.hypot(*a); c,s=a/l; R=np.array([[c,s],[-s,c]])
    a,b,va,vb=R@a,R@b,R@va,R@vb
    return np.array([l,b[0],b[1],va[0],va[1],vb[0],vb[1],1/l,1/np.hypot(*b),1/np.hypot(*(b-a))],complex)
P=lambda *z:np.array(z,float).reshape(3,2)
regimes={
    'hierarchical_binary_plus_far':(P(0,0,0.25,0,2.5,2.0),P(0,-1.0,0,2.0,-0.2,0.3)),
 'near_collision_12':(P(0,0,0.2,0,0.5,1.0),P(0,-1.2,0,2.4,0.1,0.0)),
}
out={}; rng=np.random.default_rng(99)
for name,(r,v) in regimes.items():
    X=chart(r,v); cT=F6(X); mid=0.5*(c0+cT)*(1+0.5*(rng.normal(size=7)+1j*rng.normal(size=7)))
    Ys=S.copy(); ok=np.ones(len(S),bool)
    for a,b in ((c0,mid),(mid,cT)):
        Ys,o=track_many(Ys,a,b,E,C,off,md); ok&=o
    G=Ys[ok]
    Xn=np.c_[G.real,G.imag]/(1+np.abs(G).max(1))[:,None]; prs=cKDTree(Xn).query_pairs(1e-7)
    drop={j for i,j in prs}; ndup=len(prs); G=G[[i for i in range(len(G)) if i not in drop]]
    np.save(f'/home/user/3BP-Research/data_regime_{name}.npy',G)
    a7=np.array([A7(g) for g in G])
    t=cKDTree(np.c_[a7.real,a7.imag]); dd,_=t.query(np.c_[a7.real,a7.imag],k=2); rel=dd[:,1]/(1+np.abs(a7))
    real=np.abs(G.imag).max(1)<1e-7*(1+np.abs(G).max(1)); phys=real&(G[:,0].real>0)&(G[:,8].real>0)&(G[:,9].real>0)
    src=np.abs(G-X).max(1)/(1+np.abs(X).max())
    js=int(np.argmin(src))
    pa7=np.sort(a7[phys].real)
    out[name]=dict(tracked=int(ok.sum()),failed=int((~ok).sum()),coincident_endpoint_pairs_merged=ndup,unique_endpoints=len(G),
        real=int(real.sum()),real_physical=int(phys.sum()),source_found=bool(src[js]<1e-7),
        min_rel_A7_gap_all=float(rel.min()),A7_source=float(a7[js].real),
        nearest_other_A7_rel_to_source=float((np.abs(np.delete(a7,js)-a7[js])/(1+abs(a7[js]))).min()),
        physical_A7=pa7.tolist())
    json.dump(out,open(sys.argv[2],'w'),indent=1)
    print(name,{k:v for k,v in out[name].items() if k!='physical_A7'},flush=True)
json.dump(out,open(sys.argv[2],'w'),indent=1)
