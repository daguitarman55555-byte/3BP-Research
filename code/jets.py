"""Exact polynomial area jets A^(k) on the lifted state variety.
Variables (inertial initial data, rotation fixed by a2=0 after differentiation):
 a=(a1,a2), b=(b1,b2), va=(u1,u2), vb=(w1,w2), inverse distances ha,hb,hc.
Constraints: ha^2|a|^2=1, hb^2|b|^2=1, hc^2|b-a|^2=1.  G=1."""
import flint, functools, pickle, os, sys
NAMES=['a1','a2','b1','b2','u1','u2','w1','w2','ha','hb','hc']
def build(m1,m2,m3,K):
    ctx=flint.fmpq_mpoly_ctx.get(NAMES,'deglex')
    a1,a2,b1,b2,u1,u2,w1,w2,ha,hb,hc=ctx.gens()
    m1,m2,m3=[flint.fmpq(m) for m in (m1,m2,m3)]
    c1,c2=b1-a1,b2-a2; cu1,cu2=w1-u1,w2-u2
    acc_a1=-(m1+m2)*a1*ha**3+m3*(c1*hc**3-b1*hb**3)
    acc_a2=-(m1+m2)*a2*ha**3+m3*(c2*hc**3-b2*hb**3)
    acc_b1=-(m1+m3)*b1*hb**3-m2*c1*hc**3-m2*a1*ha**3
    acc_b2=-(m1+m3)*b2*hb**3-m2*c2*hc**3-m2*a2*ha**3
    rhs=[u1,u2,w1,w2,acc_a1,acc_a2,acc_b1,acc_b2,
         -ha**3*(a1*u1+a2*u2),-hb**3*(b1*w1+b2*w2),-hc**3*(c1*cu1+c2*cu2)]
    gens=ctx.gens()
    def D(p): return sum((p.derivative(i)*rhs[i] for i in range(11)),ctx.from_dict({}))
    A=(a1*b2-a2*b1)/2
    out=[A]
    for k in range(K): out.append(D(out[-1])); print('order',k+1,'terms',len(out[-1]),file=sys.stderr,flush=True)
    return ctx,out
if __name__=='__main__':
    m=tuple(int(x) for x in sys.argv[1:4]); K=int(sys.argv[4])
    ctx,out=build(*m,K)
    # set a2=0 and store as dict strings
    res=[]
    for p in out:
        d={k:v for k,v in p.to_dict().items() if k[1]==0}
        res.append({tuple(k[:1]+k[2:]):(int(v.p),int(v.q)) for k,v in d.items()})
    pickle.dump(res,open(f'/home/user/3BP-Research/data_jets_{m[0]}{m[1]}{m[2]}_K{K}.pkl','wb'))
    print([len(r) for r in res])
