"""Write the square sixth-jet fibre system over F_p at a uniformly random F_p target, for msolve."""
import pickle, random, sys
p=int(sys.argv[1]); seed=int(sys.argv[2]); K=int(sys.argv[3]) if len(sys.argv)>3 else 6
random.seed(seed)
J=pickle.load(open('/home/user/3BP-Research/data_jets_123_K8.pkl','rb'))
V=['l','x','y','u','v','w','z','ha','hb','hc']
def term(e,c):
    mon='*'.join(f'{V[i]}^{int(k)}' if k>1 else V[i] for i,k in enumerate(e) if k)
    return f'{c}*{mon}' if mon else f'{c}'
polys=[]
for k in range(K+1):
    s=[]
    for e,(a,b) in J[k].items():
        c=a*pow(b,-1,p)%p; s.append(term(e,c))
    s.append(str((-random.randrange(1,p))%p))   # random target
    polys.append('+'.join(s))
polys+=['l*ha-1','hb^2*x^2+hb^2*y^2-1','hc^2*x^2-2*hc^2*x*l+hc^2*l^2+hc^2*y^2-1']
print(','.join(V)); print(p); print(',\n'.join(polys))
import sys as _s
print('max total degree',max(sum(int(t) for t in e) for k in range(K+1) for e in J[k]),file=_s.stderr)
