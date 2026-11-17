import numpy as np, statsmodels.api as sm, warnings; warnings.filterwarnings("ignore")
gam,eta,phi0,th=0.26,0.07,0.138,0.83; G=gam*(eta+th); a0=0.39
K=lambda a:a*(1+G)/((a-G)*(1+a)); L0=-(a0-G)/(1+G)
def sim(C,m,xi,chi,rng):
    le=np.repeat(rng.normal(0,1.5,C),m); e=np.exp(le)
    a=a0*e**xi; ph=phi0*e**(-chi)
    nbar=(a-G)/(ph*(1+a)); Q=gam*th/(1+a)
    clu=np.repeat(rng.normal(0,0.2,C),m)
    n=rng.poisson(np.clip(nbar*np.exp(clu),1e-6,50))
    lq=np.log(Q)+clu+rng.normal(0,0.45,e.size)
    X=sm.add_constant(le)
    pn=sm.GLM(n,X,family=sm.families.Poisson()).fit().params[1]
    pq=np.linalg.lstsq(X,lq,rcond=None)[0][1]
    return pn,pq
rng=np.random.default_rng(2026)
cases=[("taste only",0.146,0.0),("taste + time cost",0.146,0.1),("taste + time cost",0.146,0.3),
       ("time cost only",0.0,0.146*K(a0))]
print(f"K(a0)={K(a0):.4f}  Lambda0={L0:.4f}")
for C in (500,1000):
    null=np.array([ (lambda r:r[1]/r[0])(sim(C,20,0.0,0.146*K(a0),rng)) for _ in range(200)])
    sd0=null.std(ddof=1)
    for lab,xi,chi in cases:
        R=np.array([sim(C,20,xi,chi,rng) for _ in range(200)])
        lam=R[:,1]/R[:,0]
        xih=-R[:,1]*(1+a0)/a0; chih=R[:,0]+R[:,1]*(1+G)/(a0-G)
        true=L0*K(a0)*xi/(K(a0)*xi+chi) if (xi+chi)>0 else 0
        share=lam/L0
        print(f"C={C:5d} N={C*20:6d} {lab:18s} xi={xi:.3f} chi={chi:.3f} | "
              f"Lam' true {true:+.4f} mean {lam.mean():+.4f} sd {lam.std(ddof=1):.4f} "
              f"rej0 {(lam<-1.96*sd0).mean():.2f} | share {share.mean():.2f} | "
              f"xi^ {xih.mean():.3f} ({xih.std(ddof=1):.3f})  chi^ {chih.mean():.3f} ({chih.std(ddof=1):.3f})",flush=True)
