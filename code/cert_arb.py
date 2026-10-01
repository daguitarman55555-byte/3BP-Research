"""Arb (ball arithmetic, python-flint) Krawczyk certification for fibre points the double-precision
test could not settle, plus rigorous A^(7) enclosures. Fully rigorous interval arithmetic."""
import flint, pickle, numpy as np
from flint import acb, arb, ctx
from numeval import constraint_polys
def to_terms(path):
    J=pickle.load(open(path,'rb'))
    polys=[[(tuple(int(t) for t in e),flint.fmpq(p,q)) for e,(p,q) in Jk.items()] for Jk in J]
    for E,Cc in constraint_polys():
        polys_c=[(tuple(int(t) for t in e),flint.fmpq(int(c))) for e,c in zip(E,Cc)]
        polys.append(polys_c)
    return polys   # indices 0..8 jets, 9..11 constraints
def deriv(terms,i):
    out=[]
    for e,c in terms:
        if e[i]>0:
            f=list(e); f[i]-=1; out.append((tuple(f),c*e[i]))
    return out
def ev(terms,X,md):
    pw=[[acb(1)] for _ in X]
    for i,x in enumerate(X):
        for d in range(md): pw[i].append(pw[i][-1]*x)
    s=acb(0)
    for e,c in terms:
        t=acb(c)
        for i,k in enumerate(e):
            if k: t*=pw[i][k]
        s+=t
    return s
class System:
    def __init__(self,path):
        P=to_terms(path); self.sys=P[:7]+P[9:12]; self.a7=P[7]
        self.D=[[deriv(p,i) for i in range(10)] for p in self.sys]
        self.md=max(max(e) for p in P for e,_ in p)
    def F(self,X,c):
        return [ev(p,X,self.md)-(acb(c[k]) if k<7 else 0) for k,p in enumerate(self.sys)]
    def J(self,X): return [[ev(d,X,self.md) for d in row] for row in self.D]
def krawczyk(Sys,x,c,prec=212,rfac=(4,16,64,256)):
    ctx.prec=prec
    xt=[acb(complex(v).real,complex(v).imag) for v in x]
    cc=[acb(complex(v).real,complex(v).imag) for v in c]   # target taken as exact binary values
    # Newton polish in arb midpoints
    for _ in range(6):
        F=Sys.F(xt,cc); Jm=flint.acb_mat(Sys.J(xt)); d=Jm.solve(flint.acb_mat([[f] for f in F]))
        xt=[acb((xt[i]-d[i,0]).mid()) for i in range(10)]
    F=Sys.F(xt,cc); Jm=flint.acb_mat(Sys.J(xt))
    Jmid=flint.acb_mat([[acb(Jm[i,j].mid()) for j in range(10)] for i in range(10)])
    Y=Jmid.inv(); Y=flint.acb_mat([[acb(Y[i,j].mid()) for j in range(10)] for i in range(10)])
    YF=Y*flint.acb_mat([[f] for f in F])
    base=max(float(abs(YF[i,0]).upper()) for i in range(10))+1e-300
    for rf in rfac:
        r=arb(base*rf)
        X=[acb(xt[i].real+arb(0,r),xt[i].imag+arb(0,r)) for i in range(10)]  # box containing disc of radius r
        JX=flint.acb_mat(Sys.J(X))
        M=flint.acb_mat([[1 if i==j else 0 for j in range(10)] for i in range(10)])-Y*JX
        Dx=flint.acb_mat([[acb(arb(0,r),arb(0,r))] for _ in range(10)])
        K=flint.acb_mat([[xt[i]] for i in range(10)])-YF+M*Dx
        inside=all(abs(K[i,0].real-xt[i].real).upper()<r and abs(K[i,0].imag-xt[i].imag).upper()<r for i in range(10))
        if inside:
            a7=ev(Sys.a7,X,Sys.md)
            return True,float(r.upper())*2**0.5,complex(float(a7.real.mid()),float(a7.imag.mid())),float(a7.rad())*2+float(abs(a7.real.rad()))+float(abs(a7.imag.rad())), [complex(float(v.real.mid()),float(v.imag.mid())) for v in xt]
    return False,None,None,None,None
