import numpy as np
gam,eta,phi0,th=0.26,0.07,0.138,0.83; G=gam*(eta+th); a0=0.39
K=lambda a:a*(1+G)/((a-G)*(1+a)); L0=-(a0-G)/(1+G)
def pois(y,x,it=30):
    X=np.column_stack([np.ones_like(x),x]); b=np.array([np.log(y.mean()),0.0])
    for _ in range(it):
        mu=np.exp(X@b); g=X.T@(y-mu); H=(X*mu[:,None]).T@X; st=np.linalg.solve(H,g); b+=st
        if np.abs(st).max()<1e-10: break
    return b[1]
def sim(C,m,xi,chi,rng,sde=0.327,within=0.235,sdq=0.9):
    le=np.repeat(rng.normal(0,sde*np.sqrt(1-within),C),m)+rng.normal(0,sde*np.sqrt(within),C*m); e=np.exp(le)
    a=a0*e**xi; ph=phi0*e**(-chi); nbar=(a-G)/(ph*(1+a)); Q=gam*th/(1+a)
    clu=np.repeat(rng.normal(0,0.2,C),m)
    n=rng.poisson(nbar*np.exp(clu)).astype(float); lq=np.log(Q)+clu+rng.normal(0,sdq,e.size)
    pn=pois(n,le); x=le-le.mean(); pq=(x@(lq-lq.mean()))/(x@x); return pq/pn
rng=np.random.default_rng(5); chi0=0.146*K(a0)
print(f"sd(ln eps)=0.327 (GSS ln size sd 2.18 x (1-beta)=0.15), 6 obs/cluster, share noise sd 0.9")
for C in (6300,25000,100000):
    null=np.array([sim(C,6,0,chi0,rng) for _ in range(200)]); sd0=null.std(ddof=1)
    out=f"N={C*6:7d}  null sd {sd0:.3f} size {(null<-1.645*sd0).mean():.2f} |"
    for lab,xi,chi in [("taste",0.146,0.0),("mixed chi=.1",0.146,0.1)]:
        lam=np.array([sim(C,6,xi,chi,rng) for _ in range(200)])
        true=L0*K(a0)*xi/(K(a0)*xi+chi)
        out+=f" {lab}: true {true:+.3f} mean {lam.mean():+.3f} sd {lam.std(ddof=1):.3f} power {(lam<-1.645*sd0).mean():.2f} |"
    print(out,flush=True)
