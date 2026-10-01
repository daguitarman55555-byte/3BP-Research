import numpy as np, sys, json
from numeval import load, constraint_polys, pack, evalsys
from certify import certify
masses=sys.argv[1]; f=sys.argv[2]; outp=sys.argv[3]
polys=load(f'/home/user/3BP-Research/data_jets_{masses}_K8.pkl')
E,C,off=pack(polys[:7]+constraint_polys()); md=int(E.max())
E7,C7,off7=pack([polys[7]]); md7=int(E7.max())
d=np.load(f); S=d['S']; c=d['c0']
res=certify(S,c,E,C,off,md,E7,C7,off7,md7,evalsys)
ok=np.array([r[0] for r in res]); rr=np.array([r[1] for r in res])
a7=np.array([r[2] for r in res]); rad=np.array([r[3] for r in res])
# pairwise disjointness of balls (sup-norm) and of A7 discs
Sc=S[ok]; r=rr[ok]
from scipy.spatial import cKDTree
X=np.c_[Sc.real,Sc.imag]; tree=cKDTree(X)
pairs=tree.query_pairs(2*r.max()*1.0001+1e-300,p=np.inf)
ball_overlaps=sum(1 for i,j in pairs if np.abs(Sc[i]-Sc[j]).max()<=r[i]+r[j])
A=a7[ok]; R=rad[ok]; Z=np.c_[A.real,A.imag]; t2=cKDTree(Z)
p7=t2.query_pairs(2*R.max()*1.0001+1e-300)
a7_overlaps=sum(1 for i,j in p7 if abs(A[i]-A[j])<=R[i]+R[j])
gaps=np.sort(np.abs(A[:,None]-A[None,:])[np.triu_indices(len(A),1)])[:3] if len(A)<6000 else None
summary=dict(n_input=len(S),n_certified=int(ok.sum()),max_radius=float(np.nanmax(rr)),ball_overlaps=ball_overlaps,
             max_A7_radius=float(np.nanmax(rad)),A7_disc_overlaps=a7_overlaps,
             real_physical=int(((np.abs(Sc.imag).max(1)<=r)&(Sc[:,0].real>0)&(Sc[:,8].real>0)&(Sc[:,9].real>0)).sum()))
print(summary); json.dump(summary,open(outp,'w'),indent=1)
np.savez(outp.replace('.json','.npz'),centers=Sc,radii=r,A7=A,A7rad=R)
