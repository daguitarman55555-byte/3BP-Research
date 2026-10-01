"""Transport the fibre from masses (1,2,3) to masses (10,17,29) (path through complex masses),
then run fresh-loop and random-membership completeness tests at the new masses and check A7 separation."""
import numpy as np, sys, json, time
from track_mass import load_massvar, trackm_many
from numeval import load, constraint_polys, pack, evalsys
import track_robust as TR
from scipy.spatial import cKDTree
Em,Cm,offm=load_massvar('/home/user/3BP-Research/data_jets_massvar_K7.pkl'); mdm=int(Em.max())
d=np.load(sys.argv[1]); S=d['S']; c0=d['c0']
P2=load('/home/user/3BP-Research/data_jets_101729_K8.pkl'); E,C,off=pack(P2[:7]+constraint_polys()); md=int(E.max())
E7,C7,off7=pack([P2[7]])
F6=lambda Y:evalsys(Y,E,C,off,md)[0][:7]; A7=lambda Y:evalsys(Y,E7,C7,off7,int(E7.max()))[0][0]
rng=np.random.default_rng(17); rc=lambda n:(rng.normal(size=n)+1j*rng.normal(size=n))/np.sqrt(2)
m0=np.array([1,2,3],complex); m1=np.array([10,17,29],complex)
Y1=np.array([1.2,0.4,0.9,0.3,-0.2,0.1,0.4],complex)+0.1*rc(7); Y1=np.r_[Y1,1/Y1[0],1/np.sqrt(Y1[1]**2+Y1[2]**2),1/np.sqrt((Y1[1]-Y1[0])**2+Y1[2]**2)]
c1=F6(Y1)
pm=np.r_[0.5*(m0+m1)*(1+0.3*rc(3)),0.5*(c0+c1)*(1+0.5*rc(7))]
t=time.time(); Ys=S.copy(); ok=np.ones(len(S),bool)
for a,b in ((np.r_[m0,c0],pm),(pm,np.r_[m1,c1])):
    Ys,o=trackm_many(Ys,a,b,Em,Cm,offm,mdm); ok&=o
G=Ys[ok]; X=np.c_[G.real,G.imag]/(1+np.abs(G).max(1))[:,None]; tr=cKDTree(X)
dup=len(tr.query_pairs(1e-7)); src=np.abs(G-Y1).max(1).min()/(1+np.abs(Y1).max())
a7=np.array([A7(g) for g in G]); dd,_=cKDTree(np.c_[a7.real,a7.imag]).query(np.c_[a7.real,a7.imag],k=2)
res=dict(tracked=int(ok.sum()),failed=int((~ok).sum()),duplicate_endpoints=dup,source_state_found=bool(src<1e-7),
         min_rel_A7_gap=float((dd[:,1]/(1+np.abs(a7))).min()),time=time.time()-t)
print(res,flush=True)
# random membership tests at new masses
hits=0; tested=0; misses=[]
for k in range(800):
    q=rc(7)*rng.choice([0.3,1.0,3.0]); l,x,y=q[:3]
    Y=np.r_[q,1/l,rng.choice([1,-1])/np.sqrt(x*x+y*y),rng.choice([1,-1])/np.sqrt((x-l)**2+y*y)]
    cs=F6(Y); mid=0.5*(cs+c1)*(1+rc(7)); y=Y[None].copy(); okk=True
    for a,b in ((cs,mid),(mid,c1)):
        y,o=TR.track_many(y,a,b,E,C,off,md); okk&=bool(o[0])
    if okk:
        tested+=1; dist=tr.query(np.r_[y[0].real,y[0].imag]/(1+np.abs(y[0]).max()))[0]
        if dist<1e-7: hits+=1
        else: misses.append([str(v) for v in y[0]])
res.update(membership_tested=tested,membership_hits=hits)
print(res); json.dump(res,open('/home/user/3BP-Research/certificates/mass_101729_transport.json','w'),indent=1)
np.savez('/home/user/3BP-Research/data_fibre_101729.npz',S=G,c0=c1,a7=a7)
