"""PHASE 1 item 5 audit (independent engine): pair-Jacobian ranks at the opposite-spin equilateral
pair (total mass 9, unit side, lambda=0 (turning), omega=+-1), masses (1,2,3)*1.5.
Taylor recurrence in mpmath (60 digits) + central finite differences (step 1e-20) in the 7 chart coords."""
import mpmath as mp, json
mp.mp.dps=60
def jets(q,m,K):
    l,x,y,u,v,w,z=q; m1,m2,m3=m
    a=[[l],[mp.mpf(0)]]; b=[[x],[y]]; va=[[u],[v]]; vb=[[w],[z]]
    def series_inv_cube_norm(p,n):  # coefficients of |p|^{-3} up to order n via h' = -h^3 p.p'
        pass
    # generic Taylor via Picard iteration on truncated power series (simple, O(K^3))
    T=lambda c:[mp.mpf(0)]*(K+1) if c is None else c
    import itertools
    def mul(A,B): return [mp.fsum(A[j]*B[k-j] for j in range(k+1)) for k in range(K+1)]
    def integ(A,c0): return [c0]+[A[k-1]/k for k in range(1,K+1)]
    # state series
    S={'a1':[l]+[0]*K,'a2':[mp.mpf(0)]+[0]*K,'b1':[x]+[0]*K,'b2':[y]+[0]*K,
       'u1':[u]+[0]*K,'u2':[v]+[0]*K,'w1':[w]+[0]*K,'w2':[z]+[0]*K}
    h0={'ha':1/mp.sqrt(l*l),'hb':1/mp.sqrt(x*x+y*y),'hc':1/mp.sqrt((x-l)**2+y*y)}
    for key,val in h0.items(): S[key]=[val]+[0]*K
    for it in range(K+1):
        a1,a2,b1,b2=S['a1'],S['a2'],S['b1'],S['b2']; c1=[b1[k]-a1[k] for k in range(K+1)]; c2=[b2[k]-a2[k] for k in range(K+1)]
        u1,u2,w1,w2=S['u1'],S['u2'],S['w1'],S['w2']
        ha3=mul(mul(S['ha'],S['ha']),S['ha']); hb3=mul(mul(S['hb'],S['hb']),S['hb']); hc3=mul(mul(S['hc'],S['hc']),S['hc'])
        A1=[-(m1+m2)*p+m3*(q1-r) for p,q1,r in zip(mul(a1,ha3),mul(c1,hc3),mul(b1,hb3))]
        A2=[-(m1+m2)*p+m3*(q1-r) for p,q1,r in zip(mul(a2,ha3),mul(c2,hc3),mul(b2,hb3))]
        B1=[-(m1+m3)*p-m2*q1-m2*r for p,q1,r in zip(mul(b1,hb3),mul(c1,hc3),mul(a1,ha3))]
        B2=[-(m1+m3)*p-m2*q1-m2*r for p,q1,r in zip(mul(b2,hb3),mul(c2,hc3),mul(a2,ha3))]
        cu1=[w1[k]-u1[k] for k in range(K+1)]; cu2=[w2[k]-u2[k] for k in range(K+1)]
        dha=[-t for t in mul(ha3,[p+q1 for p,q1 in zip(mul(a1,u1),mul(a2,u2))])]
        dhb=[-t for t in mul(hb3,[p+q1 for p,q1 in zip(mul(b1,w1),mul(b2,w2))])]
        dhc=[-t for t in mul(hc3,[p+q1 for p,q1 in zip(mul(c1,cu1),mul(c2,cu2))])]
        S={'a1':integ(u1,l),'a2':integ(u2,mp.mpf(0)),'b1':integ(w1,x),'b2':integ(w2,y),
           'u1':integ(A1,u),'u2':integ(A2,v),'w1':integ(B1,w),'w2':integ(B2,z),
           'ha':integ(dha,h0['ha']),'hb':integ(dhb,h0['hb']),'hc':integ(dhc,h0['hc'])}
    Aser=[(p-q1)/2 for p,q1 in zip(mul(S['a1'],S['b2']),mul(S['a2'],S['b1']))]
    return [Aser[k]*mp.factorial(k) for k in range(K+1)]
def state(lam,om):
    s3=mp.sqrt(3); bx,by=mp.mpf(1)/2,s3/2
    return [mp.mpf(1),bx,by,lam,om,lam*bx-om*by,lam*by+om*bx]
m=[mp.mpf(3)/2,mp.mpf(3),mp.mpf(9)/2]; K=11
X,Y=state(0,1),state(0,-1); eps=mp.mpf(10)**-20
def jac(q):
    cols=[]
    for i in range(7):
        qp=list(q); qm=list(q); qp[i]+=eps; qm[i]-=eps
        cols.append([(p-r)/(2*eps) for p,r in zip(jets(qp,m,K),jets(qm,m,K))])
    return mp.matrix(cols).T
JX,JY=jac(X),jac(Y)
FX,FY=jets(X,m,K),jets(Y,m,K)
out={'jet_difference_X_minus_Y':[mp.nstr(a-b,5) for a,b in zip(FX,FY)]}
for k in (6,9,10,11):
    M=mp.matrix(k+1,14)
    for r in range(k+1):
        for c in range(7): M[r,c]=JX[r,c]; M[r,7+c]=-JY[r,c]
    sv=mp.svd_r(M,compute_uv=False)
    out[f'pair_rows_0..{k}_singular_values']=[mp.nstr(s,4) for s in sv]
print(json.dumps(out,indent=1)); json.dump(out,open('/home/user/3BP-Research/certificates/audit_equilateral_pair.json','w'),indent=1)
