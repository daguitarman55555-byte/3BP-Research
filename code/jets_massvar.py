"""Exact area jets with the masses as polynomial variables (for mass homotopies)."""
import flint, pickle, sys
NAMES=['a1','a2','b1','b2','u1','u2','w1','w2','ha','hb','hc','m1','m2','m3']
K=int(sys.argv[1])
ctx=flint.fmpq_mpoly_ctx.get(NAMES,'deglex')
a1,a2,b1,b2,u1,u2,w1,w2,ha,hb,hc,m1,m2,m3=ctx.gens()
c1,c2=b1-a1,b2-a2; cu1,cu2=w1-u1,w2-u2
rhs=[u1,u2,w1,w2,-(m1+m2)*a1*ha**3+m3*(c1*hc**3-b1*hb**3),-(m1+m2)*a2*ha**3+m3*(c2*hc**3-b2*hb**3),
     -(m1+m3)*b1*hb**3-m2*c1*hc**3-m2*a1*ha**3,-(m1+m3)*b2*hb**3-m2*c2*hc**3-m2*a2*ha**3,
     -ha**3*(a1*u1+a2*u2),-hb**3*(b1*w1+b2*w2),-hc**3*(c1*cu1+c2*cu2)]
def D(p): return sum((p.derivative(i)*rhs[i] for i in range(11)),ctx.from_dict({}))
out=[(a1*b2-a2*b1)/2]
for k in range(K): out.append(D(out[-1])); print(k+1,len(out[-1]),file=sys.stderr)
res=[{tuple(k[:1]+k[2:]):(int(v.p),int(v.q)) for k,v in p.to_dict().items() if k[1]==0} for p in out]
pickle.dump(res,open(f'/home/user/3BP-Research/data_jets_massvar_K{K}.pkl','wb'))
