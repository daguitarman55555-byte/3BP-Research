"""Compile exact jet polynomials into numba evaluators for the square system
G(Y;c) = [A_k(Y)-c_k (k=0..6), ha*l-1, hb^2(x^2+y^2)-1, hc^2((x-l)^2+y^2)-1],
Y=(l,x,y,u,v,w,z,ha,hb,hc)."""
import numpy as np, pickle, numba as nb
def load(path, K=6):
    J=pickle.load(open(path,'rb'))
    polys=[]
    for k in range(len(J)):
        items=list(J[k].items())
        E=np.array([e for e,_ in items],dtype=np.int64).reshape(-1,10)
        C=np.array([p/q for _,(p,q) in items],dtype=np.float64)
        polys.append((E,C))
    return polys
def constraint_polys():
    # l*ha-1 ; hb^2 x^2 + hb^2 y^2 -1 ; hc^2 (x^2 -2xl + l^2 + y^2) -1
    def mono(**kw):
        e=[0]*10; idx=dict(l=0,x=1,y=2,u=3,v=4,w=5,z=6,ha=7,hb=8,hc=9)
        for k,v in kw.items(): e[idx[k]]=v
        return e
    P1=([mono(l=1,ha=1),mono()],[1.,-1.])
    P2=([mono(x=2,hb=2),mono(y=2,hb=2),mono()],[1.,1.,-1.])
    P3=([mono(x=2,hc=2),mono(x=1,l=1,hc=2),mono(l=2,hc=2),mono(y=2,hc=2),mono()],[1.,-2.,1.,1.,-1.])
    return [(np.array(E,dtype=np.int64),np.array(C)) for E,C in (P1,P2,P3)]
def pack(polylist):
    off=np.zeros(len(polylist)+1,dtype=np.int64)
    for i,(E,C) in enumerate(polylist): off[i+1]=off[i]+len(C)
    E=np.vstack([E for E,_ in polylist]); C=np.concatenate([C for _,C in polylist]).astype(np.complex128)
    return E,C,off
@nb.njit(cache=True)
def evalsys(Y,E,C,off,maxdeg):
    n=Y.shape[0]; m=off.shape[0]-1
    pw=np.ones((n,maxdeg+1),dtype=np.complex128)
    for i in range(n):
        for d in range(1,maxdeg+1): pw[i,d]=pw[i,d-1]*Y[i]
    F=np.zeros(m,dtype=np.complex128); Jm=np.zeros((m,n),dtype=np.complex128)
    for p in range(m):
        for t in range(off[p],off[p+1]):
            val=C[t]
            for i in range(n): val*=pw[i,E[t,i]]
            F[p]+=val
            for i in range(n):
                e=E[t,i]
                if e>0:
                    g=C[t]*e*pw[i,e-1]
                    for j in range(n):
                        if j!=i: g*=pw[j,E[t,j]]
                    Jm[p,i]+=g
    return F,Jm
