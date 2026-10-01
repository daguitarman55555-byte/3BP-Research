"""Cross-check exact jet polynomials against direct high-accuracy integration of the inertial 3-body ODE."""
import numpy as np, pickle, math
from scipy.integrate import solve_ivp
J=pickle.load(open('/home/user/3BP-Research/data_jets_123_K8.pkl','rb'))
m=np.array([1.,2.,3.])
l,x,y,u,v,w,z=1.2,0.4,0.9,0.3,-0.2,0.1,0.4
Y=[l,x,y,u,v,w,z,1/l,1/math.hypot(x,y),1/math.hypot(x-l,y)]
jet=[sum(p/q*np.prod([Y[i]**int(e[i]) for i in range(10)]) for e,(p,q) in Jk.items()) for Jk in J]
# inertial bodies: r1=0, r2=a, r3=b, velocities 0,(u,v),(w,z)
r=np.array([[0,0],[l,0],[x,y]],float); vv=np.array([[0,0],[u,v],[w,z]],float)
def f(t,s):
    r=s[:6].reshape(3,2); acc=np.zeros((3,2))
    for i in range(3):
        for j in range(3):
            if i!=j: d=r[j]-r[i]; acc[i]+=m[j]*d/np.linalg.norm(d)**3
    return np.r_[s[6:],acc.ravel()]
def A(s): r=s[:6].reshape(3,2); a=r[1]-r[0]; b=r[2]-r[0]; return 0.5*(a[0]*b[1]-a[1]*b[0])
h=0.005; ts=np.arange(-4,5)*h
sol=solve_ivp(f,(0,4*h),np.r_[r.ravel(),vv.ravel()],rtol=1e-13,atol=1e-15,dense_output=True)
solb=solve_ivp(f,(0,-4*h),np.r_[r.ravel(),vv.ravel()],rtol=1e-13,atol=1e-15,dense_output=True)
Avals=np.array([A(sol.sol(t)) if t>=0 else A(solb.sol(t)) for t in ts])
taylor=np.array([sum(jet[k]*t**k/math.factorial(k) for k in range(9)) for t in ts])
print('max |A_ode - Taylor_8| =',np.abs(Avals-taylor).max(), ' (truncation O(h^9) ~',abs(jet[8])*(4*h)**9/math.factorial(9),')')
