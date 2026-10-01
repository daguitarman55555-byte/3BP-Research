"""Parameter homotopy of the complete (numerically saturated) fibre from c0 to a new target c1.
Also used as an independent completeness check: new endpoints must be distinct, none must fail,
and (for a closure test) applying fresh loops must produce no new points."""
import numpy as np, sys, json, pickle
from numeval import load, constraint_polys, pack, evalsys
from track import track_many
polys=load('/home/user/3BP-Research/data_jets_123_K8.pkl')
E,C,off=pack(polys[:7]+constraint_polys()); md=int(E.max())
E7,C7,off7=pack([polys[7]])
def F6(Y): return evalsys(np.asarray(Y,complex),E,C,off,md)[0][:7]
def A7(Y): return evalsys(np.asarray(Y,complex),E7,C7,off7,int(E7.max()))[0][0]
def move(S,c0,c1,mid=None):
    path=[c0]+([mid] if mid is not None else [])+[c1]; Ys=S.copy(); okall=np.ones(len(S),bool)
    for a,b in zip(path[:-1],path[1:]):
        Ys,ok=track_many(Ys,a,b,E,C,off,md); okall&=ok
    return Ys,okall
def distinct(Ys,tol=1e-7):
    from scipy.spatial import cKDTree
    X=np.c_[Ys.real,Ys.imag]/(1+np.abs(Ys).max(1))[:,None]
    return len(cKDTree(X).query_pairs(tol))
def match(A,B,tol=1e-7):
    from scipy.spatial import cKDTree
    t=cKDTree(np.c_[B.real,B.imag]); d,i=t.query(np.c_[A.real,A.imag]); return d/(1+np.abs(A).max(1)),i
