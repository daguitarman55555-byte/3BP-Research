"""Arb-certify only the sheets not already certified, then merge with the earlier certificate."""
import numpy as np, sys, pickle, time
from multiprocessing import Pool
from cert_arb import System, krawczyk
d=np.load(sys.argv[1]); S=d['S']; c=d['c0']; old=pickle.load(open(sys.argv[2],'rb')); n0=len(old['res'])
assert np.allclose(old['c0'],c)
Sys=None
def work(i):
    global Sys
    if Sys is None: Sys=System('/home/user/3BP-Research/data_jets_123_K8.pkl')
    for prec in (212,424,848):
        try:
            r=krawczyk(Sys,S[i],c,prec=prec)
            if r[0]: return (i,)+r[:4]+(prec,r[4])
        except Exception: pass
    return (i,False,None,None,None,None,None)
if __name__=='__main__':
    t=time.time()
    with Pool(4) as p: res=p.map(work,range(n0,len(S)),chunksize=20)
    allres=list(old['res'])+res
    pickle.dump(dict(res=allres,c0=c),open(sys.argv[3],'wb'))
    print('new certified',sum(1 for r in res if r[1]),'of',len(res),'total',len(allres),'time',time.time()-t)
