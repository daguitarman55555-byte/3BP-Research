"""BFS from NEW sheets only (e.g. membership misses) under a set of loops; adds every new image."""
import numpy as np, sys, time, json
from move_fibre import *
from track_robust import track_many
from scipy.spatial import cKDTree
inp=sys.argv[1]; seeds_json=sys.argv[2]; seed=int(sys.argv[3]); nl=int(sys.argv[4]); outp=sys.argv[5]
d=np.load(inp); S=list(d['S']); c0=d['c0']
rng=np.random.default_rng(seed); rc=lambda n:(rng.normal(size=n)+1j*rng.normal(size=n))/np.sqrt(2)
loops=[(c0*(1+1.0*rc(7)),c0*(1+1.0*rc(7))) for _ in range(nl)]
nrm=lambda Y:np.r_[Y.real,Y.imag]/(1+np.abs(Y).max())
tree=cKDTree(np.array([nrm(Y) for Y in S])); extra=[]
def known(Y):
    if tree.query(nrm(Y))[0]<1e-7: return True
    return bool(extra) and min(np.linalg.norm(nrm(Y)-nrm(e)) for e in extra)<1e-7
queue=[]
for f in seeds_json.split(','):
    for r in json.load(open(f))['results']:
        if r.get('ok') and not r['hit']:
            Y=np.array([complex(s) for s in r['end']])
            if not known(Y): S.append(Y); extra.append(Y); queue.append(len(S)-1)
print('seeds',len(queue)); t=time.time()
while queue:
    batch=queue[:2000]; queue=queue[2000:]
    for (p,q) in loops:
        Ys=np.array([S[i] for i in batch])
        ok=np.ones(len(batch),bool)
        for a,b in ((c0,p),(p,q),(q,c0)):
            Ys,o=track_many(Ys,a,b,E,C,off,md); ok&=o
        for Y in Ys[ok]:
            if not known(Y): S.append(Y); extra.append(Y); queue.append(len(S)-1)
    print(f'total {len(S)} new {len(extra)} queue {len(queue)} t={time.time()-t:.0f}s',flush=True)
np.savez(outp,S=np.array(S),c0=c0)
