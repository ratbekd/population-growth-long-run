import numpy as np
from scipy.optimize import brentq
gam,eta,phi,th,sig=0.26,0.07,0.138,0.83,1.0; G=gam*(eta+th); a0,xi=0.39,0.146
gbar,omega,nu=0.0975,100.0,11.0
ell=lambda a:(1+gam*th)/(1+a)
def simulate(zeta,tau=0.0,rho=2.0,Ebar=2.0,T=60,A0=1.0,h0=1.0,L0=1.0,Pbar=None):
    A,h,L=A0,h0,L0
    if Pbar is None:  # normalise so that E0=1 with tau=0 and alpha=a0
        Pbar=A0**sig*h0*ell(a0)*L0/A0**zeta
    out=[]
    for t in range(T):
        Einf=Ebar*(1+rho*tau)
        def E_of(a): P=A**sig*h*ell(a)*L/A**zeta; return Einf/(1+P/Pbar)
        F=lambda a: a-a0*(E_of(a)/L)**xi
        try: a=brentq(F,G+1e-10,50.0)
        except ValueError: a=np.nan
        E=E_of(a) if np.isfinite(a) else E_of(G)
        eps=E/L
        if not np.isfinite(a) or a<=G+1e-9:
            out.append((t,L,0.0,eps,E,np.nan,np.nan,A,h)); break
        n=(a-G)/(phi*(1+a)); s=eta*gam*phi/(a-G); qs=gam*th*phi/(a-G)
        H=h*ell(a)*L; gA=gbar/(1+omega*np.exp(-H/A))
        y=A**sig*h*ell(a)*(1-tau)
        out.append((t,L,n,eps,E,gA,y,A,h))
        h=nu*((1-tau)*qs*A**sig*h)**th*s**eta; A=A*(1+gA); L=L*n
        if L<1e-6: break
    o=np.array(out,dtype=float)
    return o
if __name__=="__main__":
    thr=sig/(1-th); print(f"decoupling threshold zeta* = sigma/(1-theta) = {thr:.3f}")
    ann=lambda x: x**(1/30)-1
    for lab,z,tau in [("decoupling zeta=7, tau=0",7,0.0),("decoupling zeta=7, tau=0.1",7,0.1),
                      ("threshold zeta=5.88",thr,0.0),("weak zeta=4, tau=0",4,0.0),("weak zeta=4, tau=0.1",4,0.1),("no decoupling zeta=0",0,0.0)]:
        o=simulate(z,tau)
        t,L,n,eps,E,gA,y=o[:,0],o[:,1],o[:,2],o[:,3],o[:,4],o[:,5],o[:,6]
        gy=np.r_[np.nan,y[1:]/y[:-1]]
        idx=[i for i in (0,5,10,20,30,40,59) if i<len(o)]
        print(f"\n{lab}  (ran {len(o)} gens)")
        for i in idx: print(f"   gen {int(t[i]):2d}: L={L[i]:9.4f} n={n[i]:.4f} eps={eps[i]:.3f} E={E[i]:.3f} gA_ann={ann(1+gA[i]) if np.isfinite(gA[i]) else np.nan:.4%} gy_ann={ann(gy[i]) if np.isfinite(gy[i]) else np.nan:.4%}")
