import pandas as pd, numpy as np
exec(open('cex_est.py').read().split("print(f\"taste-only band")[0].replace("B=250","B=150"))
def ols(y,X,w):
    Xw=X*w[:,None]; return np.linalg.solve(Xw.T@X+1e-8*np.eye(X.shape[1]),Xw.T@y)
def run2(g,specs,lab,lev=0,B=200,seed=3):
    """specs: list of (name, column, 'ppml'|'ols'); first is fertility. Ratios of hab coefs."""
    cols=[s[1] for s in specs]; g=g.dropna(subset=cols+['hab','AGE_REF','SEX_REF','FINLWT21']); g=g[np.isfinite(g[cols]).all(1)]
    X=X_of(g,lev); j=list(X.columns).index('hab'); Xv=X.values; w=g.FINLWT21.values; Y=[g[cn].values for cn in cols]
    def fit(ii,b0=None):
        out=[];bb=[]
        for k,(nm,cn,kind) in enumerate(specs):
            b=ppml(Y[k][ii],Xv[ii],w[ii],it=(80 if b0 is None else 10),b0=(None if b0 is None else b0[k])) if kind=='ppml' else ols(Y[k][ii],Xv[ii],w[ii])
            out.append(b[j]); bb.append(b)
        return np.array(out),bb
    est,b0=fit(np.arange(len(g)))
    codes,u=pd.factorize(g.clu); idx=[np.where(codes==k)[0] for k in range(len(u))]; rng=np.random.default_rng(seed); bs=[]
    for _ in range(B):
        ii=np.concatenate([idx[k] for k in rng.integers(0,len(u),len(u))])
        try: bs.append(fit(ii,b0)[0])
        except Exception: pass
    bs=np.array(bs); bs=bs[(np.abs(bs-est)<0.5).all(1)]; se=bs.std(0)
    print(f"{lab}  N={len(g)} clusters={len(u)} B={len(bs)}")
    for k,(nm,cn,kind) in enumerate(specs):
        line=f"    {nm:34s} {est[k]:+.4f} ({se[k]:.4f})"
        if k>0:
            r=bs[:,k]/bs[:,0]; lo,hi=np.percentile(r,[2.5,97.5]); line+=f"   ratio {est[k]/est[0]:+.3f} [{lo:+.3f},{hi:+.3f}]"
        print(line,flush=True)
w=c[c.wage.between(3,200)&(c.NO_EARNR<=2)&(c['T']>c.E)&(c['T']>0)].copy()
w['lnC']=np.log(w['T']-w.E); w['lnw']=np.log(w.wage); w['lnCw']=w.lnC-w.lnw; w['lnsal']=np.log(w.sal); w['lnhrs']=np.log(w.hrs)
w['lnTw']=np.log(w['T'])-w.lnw; w['anyE']=(w.E>0)*1.0
cap=w.Qw.quantile(.99); w['Qw99']=w.Qw.clip(upper=cap); w['lnSHw']=np.log(w.SH.clip(lower=1))-w.lnw
w['lnNSw']=np.log((w['T']-w.E-w.SH).clip(lower=1))-w.lnw
S=[('fertility: children <18 (PPML)','PERSLT18','ppml'),('ln(consumption/wage)','lnCw','ols'),('  ln consumption','lnC','ols'),
   ('  ln wage','lnw','ols'),('  ln annual hours','lnhrs','ols'),('ln(non-shelter consumption/wage)','lnNSw','ols'),
   ('education/wage, winsorised (PPML)','Qw99','ppml'),('any education spending (PPML)','anyE','ppml')]
run2(w,S,"HEADLINE: ref aged 25-44, 1-2 adults, earners")
S2=[S[0],S[1],S[6]]
run2(w,S2,"+ education of reference person",lev=1)
run2(w[w.MARITAL1==1],S2,"married"); run2(w[w.PERSLT18>0],[('fertility','PERSLT18','ppml'),S[1],S[6]],"households with children")
run2(w[w.yr<=2012],S2,"2000-2012"); run2(w[w.yr>2012],S2,"2013-2026"); run2(w[w.BLS_URBN==1],S2,"urban only")
c['hrs_']=c.hrs; run2(c,[S[0],('annual hours incl. zeros (PPML)','hrs_','ppml')],"LABOUR SUPPLY: all 1-2 adult households")
run2(c[c.MARITAL1==1],[S[0],('annual hours incl. zeros (PPML)','hrs_','ppml')],"LABOUR SUPPLY: married")
