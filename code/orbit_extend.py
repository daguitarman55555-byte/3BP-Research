"""Extend a fibre set by BFS under NEW random loops (orbit closure), starting from all known sheets."""
import numpy as np, sys, time
from move_fibre import *
from scipy.spatial import cKDTree
inp=sys.argv[1]; seed=int(sys.argv[2]); nl=int(sys.argv[3]); SC=float(sys.argv[4]); outp=sys.argv[5]
d=np.load(inp); S=list(d['S']); c0=d['c0']
rng=np.random.default_rng(seed); rc=lambda n:(rng.normal(size=n)+1j*rng.normal(size=n))/np.sqrt(2)
loops=[(c0*(1+SC*rc(7)),c0*(1+SC*rc(7))) for _ in range(nl)]
def norm(Y): return np.c_[Y.real,Y.imag]/(1+np.abs(Y).max(-1))[:,None]
tree=cKDTree(norm(np.array(S))); extra=[]
def known(Y):
    if tree.query(norm(Y[None])[0])[0]<1e-7: return True
    if extra and np.min(np.abs(np.array(extra)-Y).max(1)/(1+np.abs(Y).max()))<1e-7: return True
    return False
queue=list(range(len(S))); t=time.time(); fails=0
while queue:
    batch=queue[:2000]; queue=queue[2000:]
    for (p,q) in loops:
        Ys=np.array([S[i] for i in batch]); ok=np.ones(len(batch),bool)
        for a,b in ((c0,p),(p,q),(q,c0)):
            Ys,o=track_many(Ys,a,b,E,C,off,md); ok&=o
        fails+=int((~ok).sum())
        for Y in Ys[ok]:
            if not known(Y):
                S.append(Y); extra.append(Y); queue.append(len(S)-1)
    print(f'total {len(S)} new {len(extra)} queue {len(queue)} fails {fails} t={time.time()-t:.0f}s',flush=True)
    np.savez(outp,S=np.array(S),c0=c0,loops=np.array(loops))
