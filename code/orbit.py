"""Orbit closure of the sixth-jet fibre under a FIXED finite set of monodromy loops.
Loops are triangles c0->p_i->q_i->c0.  Output: set closed under all loops (if no failures)."""
import numpy as np, sys, time
sys.path.insert(0,'/home/user/3BP-Research/code')
from numeval import load, constraint_polys, pack, evalsys
from track import track_many, newton
masses=sys.argv[1]; seed=int(sys.argv[2]); nloops=int(sys.argv[3]); SC=float(sys.argv[4])
polys=load(f'/home/user/3BP-Research/data_jets_{masses}_K8.pkl')
E,C,off=pack(polys[:7]+constraint_polys()); md=int(E.max())
E7,C7,off7=pack([polys[7]]); E8,C8,off8=pack([polys[8]])
rng=np.random.default_rng(seed)
def rc(n): return (rng.normal(size=n)+1j*rng.normal(size=n))/np.sqrt(2)
l=1.2+0.1*rc(1)[0]; x,y=0.4+0.1*rc(1)[0],0.9+0.1*rc(1)[0]; u,v,w,z=np.array([0.3,-0.2,0.1,0.4])+0.1*rc(4)
Y0=np.array([l,x,y,u,v,w,z,1/l,1/np.sqrt(x*x+y*y),1/np.sqrt((x-l)**2+y*y)])
c0=evalsys(Y0,E,C,off,md)[0][:7].copy()
loops=[(c0*(1+SC*rc(7)),c0*(1+SC*rc(7))) for _ in range(nloops)]
sols=[Y0]; queue=[0]; log=open(f'/home/user/3BP-Research/logs/orbit_{masses}_s{seed}.log','w')
def find(Y,A):
    d=np.linalg.norm(A-Y,axis=1)/(1+np.linalg.norm(Y)); i=int(d.argmin()); return i if d[i]<1e-7 else -1
t0=time.time(); totfail=0; failed=[]
while queue:
    batch=queue[:400]; queue=queue[400:]
    for li,(p,q) in enumerate(loops):
        Ys=np.array([sols[i] for i in batch]); okall=np.ones(len(batch),bool)
        for (a,b) in ((c0,p),(p,q),(q,c0)):
            Ys,ok=track_many(Ys,a,b,E,C,off,md); okall&=ok
        totfail+=int((~okall).sum()); failed+= [(batch[j],li) for j in np.where(~okall)[0]]
        A=np.array(sols)
        for Y in Ys[okall]:
            if find(Y,A)<0:
                sols.append(Y); queue.append(len(sols)-1); A=np.array(sols)
    msg=f'total {len(sols)} queue {len(queue)} fails {totfail} t={time.time()-t0:.0f}s'; print(msg,flush=True); log.write(msg+'\n'); log.flush()
S=np.array(sols)
res=np.array([np.linalg.norm(evalsys(Y,E,C,off,md)[0]-np.r_[c0,0,0,0])/(1+np.abs(c0).max()) for Y in S])
cond=np.array([np.linalg.cond(evalsys(Y,E,C,off,md)[1]) for Y in S])
a7=np.array([evalsys(Y,E7,C7,off7,int(E7.max()))[0][0] for Y in S])
a8=np.array([evalsys(Y,E8,C8,off8,int(E8.max()))[0][0] for Y in S])
np.savez(f'/home/user/3BP-Research/data_orbit_{masses}_s{seed}.npz',S=S,c0=c0,a7=a7,a8=a8,Y0=Y0,loops=np.array(loops),failed=np.array(failed))
gap=np.abs(a7[:,None]-a7[None,:])/(1+np.abs(a7[:,None])); np.fill_diagonal(gap,np.inf)
msg=(f'FINAL N_found={len(S)} fails={totfail} max_res={res.max():.2e} max_cond={cond.max():.2e} '
     f'min_rel_gap_A7={gap.min():.3e}; partners of Y0 (A7 rel<1e-8): {int((np.abs(a7-a7[0])<1e-8*abs(a7[0])).sum())-1}; '
     f'real-physical count: {int(((np.abs(S.imag).max(1)<1e-9)&(S[:,0].real>0)&(S[:,8].real>0)&(S[:,9].real>0)).sum())}')
print(msg); log.write(msg+'\n')
