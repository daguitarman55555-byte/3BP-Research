"""PHASE 5: direct substitution of candidate hidden symmetries into the exact jets A^(0..8)."""
import numpy as np, sys
from numeval import load, pack, evalsys
m=np.array([1.,2.,3.]); polys=load('/home/user/3BP-Research/data_jets_123_K8.pkl')
packs=[pack([p]) for p in polys]
def jets(r,v):
    # r,v: 3x2 inertial; rotate so a along x-axis (A rotation invariant; complex-safe via SO(2,C))
    a=r[1]-r[0]; b=r[2]-r[0]; va=v[1]-v[0]; vb=v[2]-v[0]
    l=np.sqrt(a@a+0j); c,s=a[0]/l,a[1]/l; R=np.array([[c,s],[-s,c]])
    a,b,va,vb=R@a,R@b,R@va,R@vb
    Y=np.array([a[0],b[0],b[1],va[0],va[1],vb[0],vb[1],1/a[0],1/np.sqrt(b@b),1/np.sqrt((b-a)@(b-a))])
    return np.array([evalsys(Y,*P,int(P[0].max()))[0][0] for P in packs])
def com(r,v): M=m.sum(); return r-(m@r)/M, v-(m@v)/M
cr=lambda p,q:p[0]*q[1]-p[1]*q[0]
rng=np.random.default_rng(3)
def first_fail(r,v,r2,v2):
    j1,j2=jets(r,v),jets(r2,v2); rel=np.abs(j1-j2)/(1+np.abs(j1))
    k=np.where(rel>1e-9)[0]; return (int(k[0]) if len(k) else None), rel
cands={}
for trial in range(3):
    r=rng.normal(size=(3,2)); v=rng.normal(size=(3,2)); r,v=com(r,v)
    I=(m*(r**2).sum(1)).sum(); J=sum(m[i]*cr(r[i],v[i]) for i in range(3)); w=J/I
    perp=np.c_[-r[:,1],r[:,0]]
    tests={
     'spin_flip (shape velocities kept, J->-J)': (r, v-2*w*perp),
     'reflection y->-y (A->-A expected)': (r*[1,-1], v*[1,-1]),
     'time reversal v->-v': (r,-v),
     'reflection+time reversal': (r*[1,-1], -v*[1,-1]),
     'shape-velocity flip (rigid part kept)': (r, -(v-w*perp)+w*perp),
     'radial-velocity flip (dilation reversed)': (r, v-2*((r*v).sum()*1.0/ (r*r).sum())*r),
    }
    for name,(r2,v2) in tests.items():
        k,rel=first_fail(r,v,r2,v2); cands.setdefault(name,[]).append(k)
for name,ks in cands.items(): print(f'{name:45s} first differing derivative order over 3 random states: {ks}')
