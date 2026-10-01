"""Audit Phase-1 item 4: move the complete fibre to the third-stage (1,2,3) target state and locate
its 37 certified real roots; count all real physical preimages at that target."""
import numpy as np, json, fractions
from move_fibre import *
B=json.load(open('/tmp/claude-0/-home-user-3BP-Research/cf728f3e-d8fd-52fa-9da4-da8bb4e48b1d/scratchpad/in/a98694b9-Fourth_Stage_Area_History_Evidence/Fourth_Stage_Area_History/inputs/THIRD_STAGE_sixth_jet_boxes.json'))
l,x,y,u,v,w,z=B['target_state']; T=np.array([l,x,y,u,v,w,z,1/l,1/np.hypot(x,y),1/np.hypot(x-l,y)],complex)
cT=F6(T)
roots=np.array([[float(fractions.Fraction(s)) for s in r['center']] for r in B['roots']])
print('their roots shape',roots.shape)
d=np.load('/home/user/3BP-Research/data_orbit_123_ext4.npz'); S=d['S']; c0=d['c0']
rng=np.random.default_rng(11); mid=0.5*(c0+cT)*(1+0.3*(rng.normal(size=7)+1j*rng.normal(size=7)))
Ys,ok=move(S,c0,cT,mid)
print('tracked',ok.sum(),'fail',(~ok).sum(),'coincident endpoint pairs',distinct(Ys[ok]))
G=Ys[ok]; real=np.abs(G.imag).max(1)<1e-8*(1+np.abs(G).max(1))
phys=real&(G[:,0].real>0)&(G[:,8].real>0)&(G[:,9].real>0)
print('real endpoints',real.sum(),'real physical',phys.sum())
R=G[phys].real
# compare with their 37 (their coordinates: 7 state coords)
dd=np.array([np.abs(R[:,:7]-r[:7]).max(1).min() for r in roots])
print('their roots found (maxdiff<1e-6):',(dd<1e-6).sum(),'of',len(roots))
json.dump(dict(tracked=int(ok.sum()),failed=int((~ok).sum()),real=int(real.sum()),real_physical=int(phys.sum()),
   third_stage_roots_recovered=int((dd<1e-6).sum()),third_stage_roots=len(roots),
   physical_A7=[float(A7(g).real) for g in G[phys]]),open('/home/user/3BP-Research/certificates/audit_third_stage_37.json','w'),indent=1)
np.save('/home/user/3BP-Research/data_fibre_at_third_stage_target.npy',G)
