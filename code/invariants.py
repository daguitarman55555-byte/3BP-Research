"""PHASE 6: number of distinct values of classical invariants over the (near-complete) sixth-jet fibre.
#distinct values of f over the fibre is a LOWER bound for [K(f):K] (exact count if the fibre is complete)."""
import numpy as np, sys, json
d=np.load(sys.argv[1]); S=d['S']; m=np.array([1.,2.,3.]); M=m.sum()
l,x,y,u,v,w,z,ha,hb,hc=S.T
r=np.stack([np.zeros_like(l),np.zeros_like(l),l,0*l,x,y],1).reshape(-1,3,2)
vv=np.stack([0*l,0*l,u,v,w,z],1).reshape(-1,3,2)
rc=r-(m[None,:,None]*r).sum(1,keepdims=True)/M; vc=vv-(m[None,:,None]*vv).sum(1,keepdims=True)/M
T=0.5*(m[None,:]*(vc**2).sum(2)).sum(1); U=m[0]*m[1]*ha+m[0]*m[2]*hb+m[1]*m[2]*hc
inv=dict(H=T-U, J=(m[None,:]*(rc[:,:,0]*vc[:,:,1]-rc[:,:,1]*vc[:,:,0])).sum(1),
 I=(m[None,:]*(rc**2).sum(2)).sum(1), Idot=2*(m[None,:]*(rc*vc).sum(2)).sum(1), r12sq=l*l, r13sq=x*x+y*y, r23sq=(x-l)**2+y*y,
 r12=1/ha, r13=1/hb, r23=1/hc)
out={}
for k,f in inv.items():
    z=np.round(f/(1+abs(f)).max()*1e7)  # cluster at relative 1e-7
    vals=np.unique(np.round(f.real,6)+1j*np.round(f.imag,6)); out[k]=int(len(vals))
print(len(S),out); json.dump(dict(fibre_size=len(S),distinct_values=out),open(sys.argv[2],'w'),indent=1)
