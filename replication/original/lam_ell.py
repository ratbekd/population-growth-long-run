import pandas as pd, numpy as np, sys
gam,eta,phi,th=0.26,0.07,0.138,0.83; G=gam*(eta+th); a0=0.39; astar=(phi+G)/(1-phi)
L0=lambda a:-(a-G)/(1+G)
d=pd.read_pickle("gss_sub.pkl")
d=d[(d.year>=1975)&d['size'].notna()&d.vpsu.notna()].copy()
print("size==0 share:",(d['size']==0).mean().round(4))
d['sizek']=d['size'].clip(lower=0.5)
d['hab']=-np.log(d.sizek)
d['hours']=np.where(d.wrkstat.isin([1,2]),d.hrs1,np.where(d.wrkstat.isin([3,4,5,6,7,8]),0.0,np.nan))
d['sph']=np.where(d.spwrksta.isin([1,2]),d.sphrs1,np.where(d.spwrksta.isin([3,4,5,6,7,8]),0.0,np.nan))
d['hh_hours']=np.where(d.marital==1,d.hours+d.sph,np.nan)
d['psu']=d.year.astype(int).astype(str)+"_"+d.vstrat.astype(int).astype(str)+"_"+d.vpsu.astype(int).astype(str)
d['wt']=d.wtssps
d['linc']=np.log(d.realinc.where(d.realinc>0))
d['ideal']=d.chldidel.where(d.chldidel<=7)

def X_of(g,level):
    X=pd.DataFrame({'c':1.0,'hab':g.hab,'age':g.age,'age2':g.age**2/100,'fem':(g.sex==2)*1.0,'att':g.attend},index=g.index)
    for v in ['marital','race','region','year']:
        X=X.join(pd.get_dummies(g[v].astype(int),prefix=v,drop_first=True).astype(float))
    if level>=1: X['educ']=g.educ
    if level>=2: X['linc']=g.linc
    return X
def ppml(y,X,w,it=60,b0=None):
    if b0 is None:
        b=np.zeros(X.shape[1]); b[0]=np.log(max(np.average(y,weights=w),1e-6))
    else: b=b0.copy()
    for _ in range(it):
        mu=np.exp(np.clip(X@b,-20,20)); z=X@b+(y-mu)/mu; W=w*mu
        bn=np.linalg.solve((X*W[:,None]).T@X+1e-9*np.eye(X.shape[1]),(X*W[:,None]).T@z)
        if np.max(np.abs(bn-b))<1e-9: b=bn;break
        b=bn
    return b
def run(g,yn,yl,level,lab,B=200,seed=1):
    need=[yn,yl,'hab','age','sex','attend','marital','race','region','wt']+(['educ'] if level>=1 else [])+(['linc'] if level>=2 else [])
    g=g.dropna(subset=need); X=X_of(g,level); keep=['c']+[c for c in X.columns if c!='c' and X[c].std()>0]; X=X[keep]
    cols=list(X.columns); j=cols.index('hab'); Xv=X.values; w=g.wt.values; yN=g[yn].values; yL=g[yl].values
    bN=ppml(yN,Xv,w); bL=ppml(yL,Xv,w); pn=bN[j]; pl=bL[j]
    codes,u=pd.factorize(g.psu); idx=[np.where(codes==k)[0] for k in range(len(u))]
    rng=np.random.default_rng(seed); bs=[]
    for _ in range(B):
        pick=rng.integers(0,len(u),len(u)); ii=np.concatenate([idx[k] for k in pick])
        try:
            a=ppml(yN[ii],Xv[ii],w[ii],it=8,b0=bN)[j]; b=ppml(yL[ii],Xv[ii],w[ii],it=8,b0=bL)[j]; bs.append((a,b,b/a))
        except Exception: pass
    bs=np.array(bs); ci=np.percentile(bs[:,2],[2.5,97.5]); se=bs.std(0)
    print(f"{lab:36s} N={len(g):6d} cl={len(u):5d} pi_n={pn:+.4f}({se[0]:.4f}) pi_l={pl:+.4f}({se[1]:.4f}) "
          f"Lam={pl/pn:+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]",flush=True)
    return pn,pl,pl/pn,ci
print(f"targets: taste-only band [{L0(astar):.3f},{L0(a0):.3f}]; pure time cost 0; pure goods cost (-0.128,0)\n")
p=d[d.age.between(25,54)]
run(p,'childs','hours',0,"headline: age 25-54, no bad controls")
run(p,'childs','hours',1," + education")
run(p,'childs','hours',2," + education + income")
run(p[p.sex==2],'childs','hours',0,"women")
run(p[p.sex==1],'childs','hours',0,"men")
run(p[p.wrkstat.isin([1,2])],'childs','hours',0,"workers only (conditional hours)")
run(p[p.marital==1],'childs','hh_hours',0,"married: household hours (own+spouse)")
run(d[d.age.between(40,54)],'childs','hours',0,"age 40-54")
run(p[p.year<2000],'childs','hours',0,"1975-1998"); run(p[p.year>=2000],'childs','hours',0,"2000-2022")
run(p.dropna(subset=['ideal']),'ideal','hours',0,"ideal children as fertility measure")
