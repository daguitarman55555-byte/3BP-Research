"""Positive control analysis (masses 1,2,2): count A7-coincident pairs and verify they are the
swap(2,3)+reflection partners (exact symmetry) by mapping each sheet through the symmetry."""
import numpy as np, json
from scipy.spatial import cKDTree
d=np.load('/home/user/3BP-Research/data_fibre_122.npz'); G=d['S']; a7=d['a7']
z=a7/(1+np.abs(a7)); Z=np.c_[np.log(np.abs(a7)+1e-300),np.angle(a7)]
t=cKDTree(Z); prs=t.query_pairs(1e-9)
co=[(i,j) for i,j in prs if abs(a7[i]-a7[j])<=1e-9*(1+abs(a7[i]))]
def sigma(Y):
    l,x,y,u,v,w,zz,ha,hb,hc=Y
    # bodies: r1=0,r2=(l,0),r3=(x,y); swap 2<->3 then reflect (X,Y)->(X,-Y): new a=(x,-y), new b=(l,0)
    a=np.array([x,-y]); b=np.array([l,0]); va=np.array([w,-zz]); vb=np.array([u,-v])
    L=np.sqrt(a@a); c,s=a/L; R=np.array([[c,s],[-s,c]]); a,b,va,vb=R@a,R@b,R@va,R@vb
    return np.array([a[0],b[0],b[1],va[0],va[1],vb[0],vb[1],1/a[0],hc*0+1/np.sqrt(b@b) if False else ha* (1 if True else 1),0,0])
def sigma2(Y):
    l,x,y,u,v,w,zz,ha,hb,hc=Y
    a=np.array([x,-y]); b=np.array([l,0]); va=np.array([w,-zz]); vb=np.array([u,-v])
    L=np.sqrt(a@a); c,s=a/L; R=np.array([[c,s],[-s,c]]); a,b,va,vb=R@a,R@b,R@va,R@vb
    # inverse distances transform as labels: new |a|=old r13 (hb), new |b|=old r12 (ha), new |b-a|=old r23 (hc)
    out=np.array([a[0],b[0],b[1],va[0],va[1],vb[0],vb[1],0,ha,hc])
    out[7]=1/out[0]; return out, hb*out[0]   # second: consistency of branch hb*l' (=+-1)
T=cKDTree(np.c_[G.real,G.imag]/(1+np.abs(G).max(1))[:,None])
img=0; branch_ok=0
for i,j in co[:2000]:
    s,chk=sigma2(G[i])
    if abs(chk+1)<1e-6: s=np.r_[-s[:7],-s[7],s[8],s[9]]   # rotate by pi back into the h_a=1/l gauge
    dist=min(np.abs(s-G[j]).max(),np.abs(sigma2(G[j])[0]*0+s-G[j]).max())/(1+np.abs(G[j]).max())
    if dist<1e-6: img+=1
res=dict(endpoints=len(G),A7_coincident_pairs=len(co),pairs_checked=min(len(co),2000),pairs_equal_to_symmetry_image=img)
print(res); json.dump(res,open('/home/user/3BP-Research/certificates/mass_122_positive_control.json','w'),indent=1)
