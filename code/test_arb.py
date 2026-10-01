import numpy as np, time
from cert_arb import System, krawczyk
S=System('/home/user/3BP-Research/data_jets_123_K8.pkl')
d=np.load('/home/user/3BP-Research/data_fibre_123_s1.npz')
t=time.time()
for i in [0,1,2]:
    print(krawczyk(S,d['S'][i],d['c0'])[:4], time.time()-t)
