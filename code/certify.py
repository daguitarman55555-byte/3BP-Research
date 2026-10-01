"""Rigorous a-posteriori certification of sixth-jet fibre points (norm-form Krawczyk / Kantorovich test)
and enclosure of A^(7) on each certified ball.
For p = sum c_a y^a with majorant Q = sum |c_a| y^a:
  |p(x+d) - p(x)| <= Q(|x|+r) - Q(|x|)   for |d_i|<=r  (componentwise, complex).
Floating-point rounding is bounded a priori by gamma = 8*(nterms+deg+4)*u relative to the majorant
evaluated at |x| (+r), which dominates every partial sum (standard recursive-summation bound,
doubled for complex arithmetic, doubled again for safety)."""
import numpy as np, numba as nb, sys
U=2.0**-53
@nb.njit(cache=True)
def maj(E,Cabs,off,ax,md):
    n=ax.shape[0]; m=off.shape[0]-1
    pw=np.ones((n,md+1))
    for i in range(n):
        for d in range(1,md+1): pw[i,d]=pw[i,d-1]*ax[i]
    Q=np.zeros(m); DQ=np.zeros((m,n))
    for p in range(m):
        for t in range(off[p],off[p+1]):
            val=Cabs[t]
            for i in range(n): val*=pw[i,E[t,i]]
            Q[p]+=val
            for i in range(n):
                e=E[t,i]
                if e>0:
                    g=Cabs[t]*e*pw[i,e-1]
                    for j in range(n):
                        if j!=i: g*=pw[j,E[t,j]]
                    DQ[p,i]+=g
    return Q,DQ
def certify(S,c,E,C,off,md,E7,C7,off7,md7,evalsys):
    nterm=np.diff(off); gam=8*(nterm+md+4)*U; gam7=8*(len(C7)+md7+4)*U
    Cabs=np.abs(C)*(1+2*U); C7abs=np.abs(C7)*(1+2*U)
    tgt=np.r_[c,0,0,0]
    out=[]
    for x in S:
        F,J=evalsys(x,E,C,off,md); F=F-tgt
        ax=np.abs(x)
        Q0,DQ0=maj(E,Cabs,off,ax,md)
        Fbound=np.abs(F)+gam*Q0+np.abs(tgt)*U*2           # enclosure radius of true F(x)
        Y=np.linalg.inv(J)
        R=np.eye(len(x))-Y@J
        Rn=np.abs(R).sum(1).max()+ (np.abs(Y)@(gam[:,None]*DQ0)).sum(1).max()+1e-15
        Yn=np.abs(Y).sum(1).max()
        YF=(np.abs(Y)@Fbound).max()*(1+1e-12)
        ok=False; rr=np.nan
        for r in YF*np.array([1.5,2,4,8,16,64,256,1e3,1e4]):
            Q1,DQ1=maj(E,Cabs,off,ax+r,md)
            Delta=((DQ1-DQ0)*(1+gam[:,None])+gam[:,None]*DQ1)      # bound on |J(y)-J~(x)| entrywise
            Dn=(np.abs(Y)@Delta).sum(1).max()
            if YF+(Rn+Dn)*r < r*(1-1e-9) and Rn+Dn<1:
                ok=True; rr=r; break
        if ok:
            a7,_=evalsys(x,E7,C7,off7,md7); a7=a7[0]
            Q7a,_=maj(E7,C7abs,off7,ax,md7); Q7b,_=maj(E7,C7abs,off7,ax+rr,md7)
            rad7=(Q7b[0]-Q7a[0])*(1+gam7)+gam7*Q7b[0]
            out.append((True,rr,a7,rad7))
        else: out.append((False,np.nan,np.nan,np.nan))
    return out
