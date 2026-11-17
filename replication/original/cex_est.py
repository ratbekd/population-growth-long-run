import pandas as pd, numpy as np, sys
gam,eta,phi,th=0.26,0.07,0.138,0.83; G=gam*(eta+th); a0=0.39; astar=(phi+G)/(1-phi)
c=pd.read_pickle('cex_cu.pkl')
c=c[c.AGE_REF.between(25,44)&(c.adults.round()<=2)&(c.adults.round()>=1)].copy()
c['Qw']=c.E/c.wage; c['Cw']=(c['T']-c.E)/c.wage; c['Tw']=c['T']/c.wage; c['ET']=c.E/c['T'].where(c['T']>0)
def X_of(g,lev):
    X=pd.DataFrame({'c':1.0,'hab':g.hab,'age':g.AGE_REF,'age2':g.AGE_REF**2/100,'fem':(g.SEX_REF==2)*1.0},index=g.index)
    for v in ['MARITAL1','REF_RACE','REGION','yr']:
        X=X.join(pd.get_dummies(g[v].fillna(-1).astype(int),prefix=v,drop_first=True).astype(float))
    if lev>=1: X=X.join(pd.get_dummies(g.EDUC_REF.fillna(-1).astype(int),prefix='ed',drop_first=True).astype(float))
    keep=['c']+[k for k in X.columns if k!='c' and X[k].std()>0]; return X[keep]
def ppml(y,X,w,it=80,b0=None):
    if b0 is None: b=np.zeros(X.shape[1]); b[0]=np.log(max(np.average(y,weights=w),1e-9))
    else: b=b0.copy()
    for _ in range(it):
        mu=np.exp(np.clip(X@b,-30,30)); W=w*mu; z=X@b+(y-mu)/mu
        bn=np.linalg.solve((X*W[:,None]).T@X+1e-6*np.eye(X.shape[1]),(X*W[:,None]).T@z)
        if np.max(np.abs(bn-b))<1e-8: b=bn;break
        b=bn
    return b
def run(g,ys,lab,lev=0,B=250,seed=3):
    g=g.dropna(subset=ys+['hab','AGE_REF','SEX_REF','FINLWT21']); X=X_of(g,lev); j=list(X.columns).index('hab'); Xv=X.values; w=g.FINLWT21.values
    Y=[g[y].values for y in ys]; b0=[ppml(y,Xv,w) for y in Y]; est=np.array([b[j] for b in b0])
    codes,u=pd.factorize(g.clu); idx=[np.where(codes==k)[0] for k in range(len(u))]; rng=np.random.default_rng(seed); bs=[]
    for _ in range(B):
        ii=np.concatenate([idx[k] for k in rng.integers(0,len(u),len(u))])
        try: bs.append([ppml(y[ii],Xv[ii],w[ii],it=10,b0=b)[j] for y,b in zip(Y,b0)])
        except Exception: pass
    bs=np.array(bs); ok=(np.abs(bs-est)<0.5).all(1); bs=bs[ok]; se=bs.std(0)
    out=f"{lab:34s} N={len(g):6d} cl={len(u):3d} B={len(bs)} | "+" ".join(f"{y}={e:+.4f}({s:.4f})" for y,e,s in zip(ys,est,se))
    for k in range(1,len(ys)):
        r=bs[:,k]/bs[:,0]; lo,hi=np.percentile(r,[2.5,97.5]); out+=f" | L[{ys[k]}]={est[k]/est[0]:+.3f} [{lo:+.3f},{hi:+.3f}]"
    print(out,flush=True); return est,se
print(f"taste-only band [{-(astar-G)/(1+G):.3f},{-(a0-G)/(1+G):.3f}]; any cost channel: 0 for Qw and Cw\n")
w=c[c.wage.between(3,200)&(c.NO_EARNR<=2)]
run(w,['PERSLT18','Qw','Cw','Tw'],"headline (workers, wage 3-200)")
run(w,['PERSLT18','Qw','Cw','Tw'],"+ education of ref (bad control)",lev=1)
run(w[w.AGE_REF.between(30,44)],['PERSLT18','Qw','Cw'],"ref age 30-44")
run(w[w.MARITAL1==1],['PERSLT18','Qw','Cw'],"married")
run(w[w.yr<=2012],['PERSLT18','Qw','Cw'],"2000-2012"); run(w[w.yr>2012],['PERSLT18','Qw','Cw'],"2013-2026")
run(w[w.BLS_URBN==1],['PERSLT18','Qw','Cw'],"urban only (POPSIZE variation)")
run(w.dropna(subset=['ET']),['PERSLT18','ET'],"spec check: educ share of spending")
run(c,['PERSLT18','hrs'],"labour supply (all, hours incl. 0)")
run(c[c.MARITAL1==1],['PERSLT18','hrs'],"labour supply, married")
