"""Completeness test: apply FRESH random monodromy loops (new seed, different radii) to every certified
sheet; count endpoints not already in the set. Also re-track the (sheet,loop) pairs that failed in the run."""
import numpy as np, sys, json, time
from move_fibre import *
d=np.load('/home/user/3BP-Research/data_orbit_123_s7.npz'); S=d['S']; c0=d['c0']; loops=d['loops']; failed=d['failed']
import pickle; R=pickle.load(open('/home/user/3BP-Research/certificates/arb_fibre_123.pkl','rb'))['res']
S=np.array([r[6] for r in R])  # arb-polished centres
seed=int(sys.argv[1]); nl=int(sys.argv[2]); SC=float(sys.argv[3])
rng=np.random.default_rng(seed); rc=lambda n:(rng.normal(size=n)+1j*rng.normal(size=n))/np.sqrt(2)
out={}; t=time.time()
for li in range(nl):
    p,q=c0*(1+SC*rc(7)),c0*(1+SC*rc(7))
    Ys=S.copy(); ok=np.ones(len(S),bool)
    for a,b in ((c0,p),(p,q),(q,c0)):
        Ys,o=track_many(Ys,a,b,E,C,off,md); ok&=o
    dist,idx=match(Ys[ok],S)
    new=int((dist>1e-7).sum()); perm_ok=len(set(idx[dist<=1e-7]))==int((dist<=1e-7).sum())
    out[f'loop{li}']=dict(scale=SC,tracked=int(ok.sum()),failed=int((~ok).sum()),new_points=new,injective=bool(perm_ok))
    print(li,out[f'loop{li}'],f'{time.time()-t:.0f}s',flush=True)
json.dump(out,open(f'/home/user/3BP-Research/certificates/closure_test_seed{seed}_sc{SC}.json','w'),indent=1)
