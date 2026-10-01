import numpy as np, sys, pickle, time
from multiprocessing import Pool
from cert_arb import System, krawczyk
f=sys.argv[1]; out=sys.argv[2]
d=np.load(f); S=d['S']; c=d['c0']
Sys=None
def work(i):
    global Sys
    if Sys is None: Sys=System('/home/user/3BP-Research/data_jets_123_K8.pkl')
    for prec in (212,424):
        try:
            r=krawczyk(Sys,S[i],c,prec=prec)
            if r[0]: return (i,)+r[:4]+(prec,r[4])
        except Exception as e: pass
    return (i,False,None,None,None,None,None)
if __name__=='__main__':
    t=time.time()
    with Pool(4) as p: res=p.map(work,range(len(S)),chunksize=50)
    pickle.dump(dict(res=res,c0=c),open(out,'wb'))
    print('certified',sum(1 for r in res if r[1]),'of',len(S),'time',time.time()-t)
