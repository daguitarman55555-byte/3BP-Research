"""Monodromy saturation of the complex sixth-jet fibre F_6^{-1}(c0) on the lifted variety V."""
import numpy as np, sys, time, json
sys.path.insert(0,'/home/user/3BP-Research/code')
from numeval import load, constraint_polys, pack, evalsys
from track import track_many, newton
masses=sys.argv[1]; seed=int(sys.argv[2]); maxloops=int(sys.argv[3]); stopafter=int(sys.argv[4]); SC=float(sys.argv[5])
polys=load(f'/home/user/3BP-Research/data_jets_{masses}_K8.pkl')
P6=polys[:7]+constraint_polys(); E,C,off=pack(P6); md=int(E.max())
E7,C7,off7=pack([polys[7]])
rng=np.random.default_rng(seed)
def rc(n,s=1.0): return s*(rng.normal(size=n)+1j*rng.normal(size=n))/np.sqrt(2)
# random complex base point on V
l=1.2+0.1*rc(1)[0]; x,y=0.4+0.1*rc(1)[0],0.9+0.1*rc(1)[0]; u,v,w,z=np.array([0.3,-0.2,0.1,0.4])+0.1*rc(4)
hb=1/np.sqrt(x*x+y*y); hc=1/np.sqrt((x-l)**2+y*y)
Y0=np.array([l,x,y,u,v,w,z,1/l,hb,hc],dtype=np.complex128)
F,_=evalsys(Y0,E,C,off,md); c0=F[:7].copy()
sols=[Y0]
def key_match(Y,S):
    if not S: return False
    A=np.array(S); d=np.linalg.norm(A-Y,axis=1)/(1+np.linalg.norm(Y))
    return d.min()<1e-7
log=open(f'/home/user/3BP-Research/logs/monodromy_{masses}_s{seed}.log','w')
noinc=0; t0=time.time()
for loop in range(maxloops):
    c1=c0*(1+SC*rc(7)); c2=c0*(1+SC*rc(7))
    Ys=np.array(sols); nfail=0
    for (a,b) in ((c0,c1),(c1,c2),(c2,c0)):
        Ys,ok=track_many(Ys,a,b,E,C,off,md); nfail+=int((~ok).sum()); Ys=Ys[ok]
    new=0
    for Y in Ys:
        if not key_match(Y,sols): sols.append(Y); new+=1
    noinc = 0 if new else noinc+1
    msg=f'loop {loop} tracked {len(Ys)+nfail} fail {nfail} new {new} total {len(sols)} t={time.time()-t0:.0f}s'
    print(msg,flush=True); log.write(msg+'\n'); log.flush()
    if noinc>=stopafter: break
S=np.array(sols)
# A7 values
a7=np.array([evalsys(Y,E7,C7,off7,int(E7.max()))[0][0] for Y in S])
np.savez(f'/home/user/3BP-Research/data_fibre_{masses}_s{seed}.npz',S=S,c0=c0,a7=a7,Y0=Y0)
d=np.abs(a7[:,None]-a7[None,:])+np.eye(len(a7))*1e300
msg=f'FINAL N_found={len(S)}  min |A7_i-A7_j| = {d.min():.3e}  scale {np.abs(a7).max():.3e}; #sols with |A7-A7(Y0)|<1e-6 rel: {(np.abs(a7-a7[0])<1e-6*np.abs(a7[0])).sum()}'
print(msg); log.write(msg+'\n')
