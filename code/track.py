import numpy as np, numba as nb
from numeval import evalsys
@nb.njit(cache=True)
def G(Y,c,E,C,off,md):
    F,J=evalsys(Y,E,C,off,md)
    for k in range(c.shape[0]): F[k]-=c[k]
    return F,J
@nb.njit(cache=True)
def newton(Y,c,E,C,off,md,it,tol):
    for _ in range(it):
        F,J=G(Y,c,E,C,off,md)
        dY=np.linalg.solve(J,F)
        Y=Y-dY
        if np.linalg.norm(dY)<=tol*(1+np.linalg.norm(Y)): return Y,True,np.linalg.norm(dY)
    return Y,False,np.linalg.norm(dY)
@nb.njit(cache=True)
def tangent(Y,c,dc,E,C,off,md):
    F,J=G(Y,c,E,C,off,md)
    rhs=np.zeros(Y.shape[0],dtype=np.complex128)
    rhs[:dc.shape[0]]=dc
    return np.linalg.solve(J,rhs)
@nb.njit(cache=True)
def track(Y,c0,c1,E,C,off,md):
    s=0.0; h=0.02; dc=c1-c0; nsteps=0
    while s<1.0:
        if h<1e-13 or nsteps>200000: return Y,False
        if s+h>1.0: h=1.0-s
        k1=tangent(Y,c0+s*dc,dc,E,C,off,md)
        k2=tangent(Y+0.5*h*k1,c0+(s+0.5*h)*dc,dc,E,C,off,md)
        k3=tangent(Y+0.5*h*k2,c0+(s+0.5*h)*dc,dc,E,C,off,md)
        k4=tangent(Y+h*k3,c0+(s+h)*dc,dc,E,C,off,md)
        Yp=Y+h/6*(k1+2*k2+2*k3+k4)
        # corrector: require quick contraction
        Yc=Yp.copy(); ok=False; prev=1e300
        for it in range(4):
            F,J=G(Yc,c0+(s+h)*dc,E,C,off,md)
            d=np.linalg.solve(J,F); nd=np.linalg.norm(d)
            Yc=Yc-d
            if it>0 and nd>0.5*prev: break
            prev=nd
            if nd<1e-9*(1+np.linalg.norm(Yc)): ok=True; break
        if ok and np.linalg.norm(Yc-Yp)<0.01*(1+np.linalg.norm(Y)):
            Y=Yc; s+=h; h=min(h*1.6,0.1); nsteps+=1
            if np.linalg.norm(Y)>1e12: return Y,False
        else:
            h*=0.5
    Y,ok,_=newton(Y,c1,E,C,off,md,10,1e-10)
    return Y,ok
@nb.njit(cache=True,parallel=True)
def track_many(Ys,c0,c1,E,C,off,md):
    n=Ys.shape[0]; out=np.empty_like(Ys); ok=np.zeros(n,dtype=np.bool_)
    for i in nb.prange(n):
        y,o=track(Ys[i].copy(),c0,c1,E,C,off,md); out[i]=y; ok[i]=o
    return out,ok
