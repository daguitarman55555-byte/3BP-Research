"""Exact consequences of the arb certificates: distinct roots, pairwise-disjoint A7 enclosures."""
import pickle, numpy as np, json, sys
from scipy.spatial import cKDTree
D=pickle.load(open(sys.argv[1],'rb')); res=D['res']
ok=np.array([r[1] for r in res]); assert ok.all()
rad=np.array([r[2] for r in res]); a7=np.array([r[3] for r in res]); a7r=np.array([r[4] for r in res])
X=np.array([r[6] for r in res]); prec=np.array([r[5] for r in res])
# balls: each coordinate disc radius <= rad (sup norm). distinct iff some coordinate separates by > rad_i+rad_j
T=cKDTree(np.c_[X.real,X.imag]); pairs=T.query_pairs(2*rad.max()*2+1e-300,p=np.inf)
ov=[(i,j) for i,j in pairs if np.abs(X[i]-X[j]).max()<=rad[i]+rad[j]]
T7=cKDTree(np.c_[a7.real,a7.imag]); p7=T7.query_pairs(2*a7r.max()*1.5+1e-300)
ov7=[(i,j) for i,j in p7 if abs(a7[i]-a7[j])<=a7r[i]+a7r[j]]
# minimal A7 separation (relative)
dd,ii=T7.query(np.c_[a7.real,a7.imag],k=2); gap=dd[:,1]
relgap=gap/(1+np.abs(a7))
s=dict(certified=int(ok.sum()),precisions_used={int(p):int((prec==p).sum()) for p in set(prec)},
 max_ball_radius=float(rad.max()),ball_overlaps=len(ov),max_A7_radius=float(a7r.max()),A7_overlaps=len(ov7),
 min_abs_A7_gap=float(gap.min()),min_rel_A7_gap=float(relgap.min()),min_gap_over_radius=float((gap/(a7r+a7r[ii[:,1]])).min()),
 base_point_index=0, real_points=int((np.abs(X.imag).max(1)<=rad).sum()))
print(json.dumps(s,indent=1)); json.dump(s,open(sys.argv[2],'w'),indent=1)
