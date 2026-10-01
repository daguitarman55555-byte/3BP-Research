import numpy as np,sys
from numeval import *; from track import *
polys=load('/home/user/3BP-Research/data_jets_123_K8.pkl')
E,C,off=pack(polys[:7]+constraint_polys()); md=int(E.max())
rng=np.random.default_rng(1)
Y0=np.array([1.2,0.4,0.9,0.3,-0.2,0.1,0.4,0,0,0],dtype=complex); Y0[7]=1/Y0[0];Y0[8]=1/np.sqrt(Y0[1]**2+Y0[2]**2);Y0[9]=1/np.sqrt((Y0[1]-Y0[0])**2+Y0[2]**2)
F,J=evalsys(Y0,E,C,off,md); print(F); print(np.linalg.cond(J))
c0=F[:7].copy()
for sc in [1e-3,1e-2,1e-1,0.3,1]:
    c1=c0*(1+sc*(rng.normal(size=7)+1j*rng.normal(size=7)))
    Y,ok=track(Y0.copy(),c0,c1,E,C,off,md); print(sc,ok,np.round(Y[:3],4))
