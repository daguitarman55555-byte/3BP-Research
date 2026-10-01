"""PHASE 1 audit: exact rank of DF_6 (and DF_7 shape) at a rational physical point.
Heron triangle 14-13-15: a=(14,0), b=(5,12): |a|=14,|b|=13,|b-a|=15 -> ha,hb,hc rational."""
import flint, sys, json
sys.path.insert(0,'.')
from jets import build
Q=flint.fmpq
def run(m,vel):
    ctx,P=build(*m,7)
    l,x,y=Q(14),Q(5),Q(12); u,v,w,z=[Q(t) for t in vel]
    ha,hb,hc=Q(1,14),Q(1,13),Q(1,15)
    pt=[l,Q(0),x,y,u,v,w,z,ha,hb,hc]          # a1,a2,b1,b2,u1,u2,w1,w2,ha,hb,hc
    # reduced coordinates q=(l,x,y,u,v,w,z) -> var indices; chain rule through h's
    idx={'l':0,'x':2,'y':3,'u':4,'v':5,'w':6,'z':7}
    dh={ # partial derivatives of (ha,hb,hc) wrt reduced coords
     'l':{8:-ha**2,10:hc**3*(x-l)},'x':{9:-hb**3*x,10:-hc**3*(x-l)},'y':{9:-hb**3*y,10:-hc**3*y}}
    rows=[]
    for k,p in enumerate(P):
        row=[]
        for name,i in idx.items():
            val=p.derivative(i)(*pt)
            for j,c in dh.get(name,{}).items(): val+=p.derivative(j)(*pt)*c
            row.append(val)
        rows.append(row)
    M=flint.fmpq_mat(rows[:7]); M8=flint.fmpq_mat(rows)
    return M.det(), M8.rank(), [p(*pt) for p in P]
out={}
for m,vel in [((1,2,3),(1,-2,3,1)),((3,5,7),(Q(1,3),Q(2,5),Q(-1,2),Q(3,7)))]:
    det,r8,vals=run(m,vel)
    out[str(m)]={'det_DF6':str(det),'rank_DF7':r8,'jet':[str(v) for v in vals]}
    print(m,'det DF6 =',float(det),' exact nonzero:',det!=0,' rank DF7:',r8)
json.dump(out,open('/home/user/3BP-Research/certificates/exact_rank_DF6.json','w'),indent=1)
