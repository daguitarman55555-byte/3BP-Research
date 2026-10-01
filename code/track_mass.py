"""Parameter homotopy in (masses, target):  P(Y; m) = c.  Used to transport the fibre between mass triples."""
import numpy as np, numba as nb, pickle
from numeval import evalsys, constraint_polys, pack
def load_massvar(path):
    J=pickle.load(open(path,'rb'))[:7]; polys=[]
    for Jk in J:
        it=list(Jk.items()); polys.append((np.array([e for e,_ in it],dtype=np.int64).reshape(-1,13),np.array([p/q for _,(p,q) in it])))
    for E,C in constraint_polys(): polys.append((np.c_[E,np.zeros((len(E),3),dtype=np.int64)],C))
    return pack(polys)
@nb.njit(cache=True)
def Gm(Y,p,E,C,off,md):
    Z=np.concatenate((Y,p[:3])); F,J=evalsys(Z,E,C,off,md)
    for k in range(7): F[k]-=p[3+k]
    return F,J
@nb.njit(cache=True)
def tang(Y,p,dp,E,C,off,md):
    F,J=Gm(Y,p,E,C,off,md); rhs=-(J[:,10:13]@dp[:3]).astype(np.complex128)
    for k in range(7): rhs[k]+=dp[3+k]
    return np.linalg.solve(J[:,:10].copy(),rhs)
@nb.njit(cache=True)
def trackm(Y,p0,p1,E,C,off,md,hmax):
    s=0.0; h=min(0.02,hmax); dp=p1-p0; n=0
    while s<1.0:
        if h<1e-13 or n>200000: return Y,False
        if s+h>1.0: h=1.0-s
        k1=tang(Y,p0+s*dp,dp,E,C,off,md); k2=tang(Y+0.5*h*k1,p0+(s+0.5*h)*dp,dp,E,C,off,md)
        k3=tang(Y+0.5*h*k2,p0+(s+0.5*h)*dp,dp,E,C,off,md); k4=tang(Y+h*k3,p0+(s+h)*dp,dp,E,C,off,md)
        Yp=Y+h/6*(k1+2*k2+2*k3+k4); Yc=Yp.copy(); ok=False; prev=1e300
        for it in range(4):
            F,J=Gm(Yc,p0+(s+h)*dp,E,C,off,md); d=np.linalg.solve(J[:,:10].copy(),F); nd=np.linalg.norm(d); Yc=Yc-d
            if it>0 and nd>0.5*prev: break
            prev=nd
            if nd<1e-9*(1+np.linalg.norm(Yc)): ok=True; break
        if ok and np.linalg.norm(Yc-Yp)<0.01*(1+np.linalg.norm(Y)):
            Y=Yc; s+=h; h=min(h*1.6,hmax); n+=1
            if np.linalg.norm(Y)>1e12: return Y,False
        else: h*=0.5
    for it in range(10):
        F,J=Gm(Y,p1,E,C,off,md); d=np.linalg.solve(J[:,:10].copy(),F); Y=Y-d
        if np.linalg.norm(d)<1e-10*(1+np.linalg.norm(Y)): return Y,True
    return Y,False
@nb.njit(cache=True,parallel=True)
def trackm_many(Ys,p0,p1,E,C,off,md):
    out=np.empty_like(Ys); ok=np.zeros(Ys.shape[0],dtype=np.bool_)
    for i in nb.prange(Ys.shape[0]):
        y,o=trackm(Ys[i].copy(),p0,p1,E,C,off,md,0.05)
        if not o: y,o=trackm(Ys[i].copy(),p0,p1,E,C,off,md,0.005)
        out[i]=y; ok[i]=o
    return out,ok
