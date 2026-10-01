"""Robust re-transport to the third-stage target (exact rational sixth jet) + arb certification of all
real physical endpoints at the EXACT rational target; comparison with the 37 third-stage roots."""
import numpy as np, json, fractions, flint
from move_fibre import F6, A7, E, C, off, md
from track_robust import track_many
from scipy.spatial import cKDTree
from cert_arb import System, ev
import cert_arb
B=json.load(open('/tmp/claude-0/-home-user-3BP-Research/cf728f3e-d8fd-52fa-9da4-da8bb4e48b1d/scratchpad/in/a98694b9-Fourth_Stage_Area_History_Evidence/Fourth_Stage_Area_History/inputs/THIRD_STAGE_sixth_jet_boxes.json'))
l,x,y,u,v,w,z=B['target_state']; T=np.array([l,x,y,u,v,w,z,1/l,1/np.hypot(x,y),1/np.hypot(x-l,y)],complex)
Sys=System('/home/user/3BP-Research/data_jets_123_K8.pkl')
Q=flint.fmpq; pt=[Q(l),Q(x),Q(y),Q(u),Q(v),Q(w),Q(z),Q(1,l),Q(1,4),Q(1,5)]
flint.ctx.prec=300
cexact=[ev(Sys.sys[k],[flint.acb(p) for p in pt],Sys.md) for k in range(7)]   # exact rationals as balls of radius 0
cT=np.array([complex(float(c.real.mid())) for c in cexact])
roots=np.array([[float(fractions.Fraction(s)) for s in r['center']] for r in B['roots']])
d=np.load('/home/user/3BP-Research/data_orbit_123_ext4.npz'); S=d['S']; c0=d['c0']
rng=np.random.default_rng(11); mid=0.5*(c0+cT)*(1+0.3*(rng.normal(size=7)+1j*rng.normal(size=7)))
Ys=S.copy(); ok=np.ones(len(S),bool)
for a,b in ((c0,mid),(mid,cT)):
    Ys,o=track_many(Ys,a,b,E,C,off,md); ok&=o
G=Ys[ok]; real=np.abs(G.imag).max(1)<1e-7*(1+np.abs(G).max(1)); phys=real&(G[:,0].real>0)&(G[:,8].real>0)&(G[:,9].real>0)
R=G[phys].real
# also add the third-stage roots themselves as candidates (they are certified there), dedupe
cand=np.vstack([R,np.c_[roots,1/roots[:,0],1/np.hypot(roots[:,1],roots[:,2]),1/np.hypot(roots[:,1]-roots[:,0],roots[:,2])]])
# arb certify at exact target
def kraw(xv):
    flint.ctx.prec=212
    xt=[flint.acb(float(t)) for t in xv]
    for _ in range(8):
        F=[ev(Sys.sys[k],xt,Sys.md)-(cexact[k] if k<7 else 0) for k in range(10)]
        Jm=flint.acb_mat(Sys.J(xt)); dd=Jm.solve(flint.acb_mat([[f] for f in F])); xt=[flint.acb((xt[i]-dd[i,0]).mid()) for i in range(10)]
    F=[ev(Sys.sys[k],xt,Sys.md)-(cexact[k] if k<7 else 0) for k in range(10)]
    Jm=flint.acb_mat(Sys.J(xt)); Yi=flint.acb_mat([[flint.acb(Jm[i,j].mid()) for j in range(10)] for i in range(10)]).inv()
    Yi=flint.acb_mat([[flint.acb(Yi[i,j].mid()) for j in range(10)] for i in range(10)])
    YF=Yi*flint.acb_mat([[f] for f in F]); base=max(float(abs(YF[i,0]).upper()) for i in range(10))+1e-300
    for rf in (4,16,64,256):
        r=flint.arb(base*rf); X=[flint.acb(xt[i].real+flint.arb(0,r),xt[i].imag+flint.arb(0,r)) for i in range(10)]
        M=flint.acb_mat([[1 if i==j else 0 for j in range(10)] for i in range(10)])-Yi*flint.acb_mat(Sys.J(X))
        K=flint.acb_mat([[xt[i]] for i in range(10)])-YF+M*flint.acb_mat([[flint.acb(flint.arb(0,r),flint.arb(0,r))] for _ in range(10)])
        if all(abs(K[i,0].real-xt[i].real).upper()<r and abs(K[i,0].imag-xt[i].imag).upper()<r for i in range(10)):
            # real root: the unique root in X is real since the system is real and X is conjugation-symmetric
            a7=ev(Sys.a7,X,Sys.md)
            return True,[float(t.real.mid()) for t in xt],float(r),float(a7.real.mid()),float(a7.rad())+float(abs(a7.real.rad()))
    return False,None,None,None,None
certs=[kraw(cv) for cv in cand]
good=[c for c in certs if c[0]]
P=np.array([c[1] for c in good]); rads=np.array([c[2] for c in good])
# dedupe certified boxes
keep=[]
for i in range(len(P)):
    if all(np.abs(P[i]-P[j]).max()>rads[i]+rads[j] for j in keep): keep.append(i)
P=P[keep]; A7v=np.array([good[i][3] for i in keep]); A7r=np.array([good[i][4] for i in keep])
phys_ok=(P[:,0]>0)&(P[:,8]>0)&(P[:,9]>0)
o=np.argsort(A7v); disj=all(A7v[o[i+1]]-A7v[o[i]]>A7r[o[i]]+A7r[o[i+1]] for i in range(len(o)-1))
inR=[bool(np.abs(P[:,:7]-r).max(1).min()<1e-6) for r in roots]
res=dict(tracked=int(ok.sum()),failed=int((~ok).sum()),real_endpoints=int(real.sum()),real_physical_endpoints=int(phys.sum()),
   certified_distinct_real_physical=int(phys_ok.sum()),A7_pairwise_disjoint=bool(disj),third_stage_roots_among_certified=int(sum(inR)),
   third_stage_roots_found_by_transport=int(sum(bool(np.abs(R[:,:7]-r).max(1).min()<1e-6) for r in roots)) if len(R) else 0)
print(res); json.dump(res,open('/home/user/3BP-Research/certificates/audit_third_stage_target.json','w'),indent=1)
np.save('/home/user/3BP-Research/certificates/third_stage_target_real_physical_certified.npy',P)
