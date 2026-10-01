import numpy as np
from numeval import *; import track as T, numba as nb
polys=load('/home/user/3BP-Research/data_jets_123_K8.pkl')
E,C,off=pack(polys[:7]+constraint_polys()); md=int(E.max())
d=np.load('/home/user/3BP-Research/data_fibre_123_s1.npz'); S=d['S']; c0=d['c0']
rng=np.random.default_rng(5)
p=c0*(1+(rng.normal(size=7)+1j*rng.normal(size=7))/np.sqrt(2))
Ys,ok=T.track_many(S[:200].copy(),c0,p,E,C,off,md)
print('fail',(~ok).sum())
bad=np.where(~ok)[0]
for i in bad[:8]:
    print(i, np.linalg.norm(S[i]), np.linalg.norm(Ys[i]), np.round(np.abs(Ys[i]),2))
