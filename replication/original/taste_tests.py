import pandas as pd, numpy as np, warnings; warnings.filterwarnings("ignore")
def ppml(y,X,w,it=80,b0=None):
    b=np.zeros(X.shape[1]) if b0 is None else b0.copy()
    if b0 is None: b[0]=np.log(max(np.average(y,weights=w),1e-9))
    for _ in range(it):
        mu=np.exp(np.clip(X@b,-30,30)); W=w*mu; z=X@b+(y-mu)/mu
        bn=np.linalg.solve((X*W[:,None]).T@X+1e-8*np.eye(X.shape[1]),(X*W[:,None]).T@z)
        if np.max(np.abs(bn-b))<1e-8: return bn
        b=bn
    return b
def ols(y,X,w):
    Xw=X*w[:,None]; return np.linalg.solve(Xw.T@X+1e-8*np.eye(X.shape[1]),Xw.T@y)
def boot(fn,g,clu,B=200,seed=1):
    codes,u=pd.factorize(g[clu]); idx=[np.where(codes==k)[0] for k in range(len(u))]; rng=np.random.default_rng(seed); est=fn(g.reset_index(drop=True)); bs=[]
    for _ in range(B):
        ii=np.concatenate([idx[k] for k in rng.integers(0,len(u),len(u))]); bs.append(fn(g.iloc[ii].reset_index(drop=True)))
    return np.array(est),np.array(bs)

# ================= TEST 1: GSS childhood habitat vs current habitat =================
print("="*90); print("TEST 1  GSS: does habitat at age 16 predict fertility, holding current size of place fixed?")
a=pd.read_pickle("gss_sub.pkl"); e=pd.read_pickle("gss_extra.pkl")[['year','id','sibs','paeduc','maeduc','mobile16']]
d=a.merge(e,on=['year','id'],how='left')
d=d[(d.year>=1975)&d['size'].notna()&d.vpsu.notna()&d.res16.between(1,6)].copy()
d['hab']=-np.log(d['size'].clip(lower=0.5))
# RES16 ordered from roomiest to most crowded; habitat at 16 = -log of representative population
res_pop={2:0.5,1:2.0,3:20.0,4:120.0,5:250.0,6:1000.0}   # farm, open country, town<50k, city 50-250k, suburb of big city, big city>250k (thousands)
d['hab16']=-np.log(d.res16.map(res_pop))
d['psu']=d.year.astype(int).astype(str)+"_"+d.vstrat.astype(int).astype(str)+"_"+d.vpsu.astype(int).astype(str)
d['wt']=d.wtssps; d['ideal']=d.chldidel.where(d.chldidel<=7)
d['moved']=(d.mobile16.isin([2,3])).astype(float)   # different city or state from age 16
def X_of(g,extra):
    X=pd.DataFrame({'c':1.0,'age':g.age,'age2':g.age**2/100,'fem':(g.sex==2)*1.0,'att':g.attend},index=g.index)
    for k in extra: X[k]=g[k]
    for v in ['marital','race','region','year']: X=X.join(pd.get_dummies(g[v].astype(int),prefix=v,drop_first=True).astype(float))
    return X[['c']+[k for k in X.columns if k!='c' and X[k].std()>0]]
def spec(g,y,extra,est='ppml'):
    X=X_of(g,extra); b=(ppml if est=='ppml' else ols)(g[y].values,X.values,g.wt.values); cols=list(X.columns)
    return [b[cols.index(k)] for k in extra]
p=d[d.age.between(40,64)].dropna(subset=['childs','hab','hab16','age','sex','attend','marital','race','region'])
for lab,g,extra in [("current habitat only",p,['hab']),("habitat at 16 only",p,['hab16']),("both",p,['hab','hab16']),
                    ("both + siblings + parents' educ",p.dropna(subset=['sibs','paeduc','maeduc']).assign(lsibs=lambda x:np.log1p(x.sibs)),['hab','hab16','lsibs','paeduc','maeduc']),
                    ("movers only (left city of age 16)",p[p.moved==1],['hab','hab16']),
                    ("stayers only (same city as age 16)",p[p.mobile16==1],['hab','hab16'])]:
    est,bs=boot(lambda gg: spec(gg,'childs',extra),g,'psu',B=150)
    se=bs.std(0); print(f"  completed fertility, age 40-64 | {lab:36s} N={len(g):6d} | "+"  ".join(f"{k}={v:+.4f}({s:.4f})" for k,v,s in zip(extra[:2],est[:2],se[:2])),flush=True)
q=d[d.age.between(25,64)].dropna(subset=['ideal','hab','hab16','age','sex','attend','marital','race','region'])
est,bs=boot(lambda gg: spec(gg,'ideal',['hab','hab16']),q,'psu',B=150); se=bs.std(0)
print(f"  ideal number of children, 25-64 | both                            N={len(q):6d} | hab={est[0]:+.4f}({se[0]:.4f})  hab16={est[1]:+.4f}({se[1]:.4f})",flush=True)

# ================= TEST 2: WVS ideal vs realised gap =================
print("="*90); print("TEST 2  WVS waves 2-4: does the gap between ideal and realised fertility widen with density?")
w=pd.read_pickle("/home/claude/work/ivs_clean.pkl")
for c in ['D017','X011','X049','F028','X007']: w[c]=w[c].where(w[c]>=0)
w['age']=w.X003.where(w.X003.between(15,99)); w['sex']=w.X001.where(w.X001.isin([1,2]))
w['cw']=w.S003.astype(str)+"_"+w.S002.astype(str); w['wt']=w.S017.where(w.S017>0,1.0)
mid={1:1000,2:3162,3:7071,4:14142,5:31623,6:70711,7:223607,8:1000000}; w['hab']=-np.log(w.X049.map(mid)/1e3)
s=w[w.S002.isin([2,3,4])&w.age.between(45,70)].dropna(subset=['D017','X011','hab','sex','X007','F028']).copy()
s=s[s.D017<=8]; s['X011c']=s.X011.clip(upper=5); s['D017c']=s.D017.clip(upper=5)   # same top-code on both
s['gap']=s.D017c-s.X011c
def Xw(g):
    X=pd.DataFrame({'c':1.0,'hab':g.hab,'age':g.age,'age2':g.age**2/100,'fem':(g.sex==2)*1.0,'rel':g.F028},index=g.index)
    X=X.join(pd.get_dummies(g.X007.astype(int),prefix='m',drop_first=True).astype(float)).join(pd.get_dummies(g.cw,prefix='cw',drop_first=True).astype(float))
    return X
def wfn(g):
    X=Xw(g); j=list(X.columns).index('hab'); Xv=X.values; ww=g.wt.values
    bi=ppml(g.D017c.values,Xv,ww)[j]; br=ppml(g.X011c.values,Xv,ww)[j]; bg=ols(g.gap.values,Xv,ww)[j]
    return [bi,br,bi/br,bg]
est,bs=boot(wfn,s,'cw',B=150); se=bs.std(0); ci=np.percentile(bs[:,2],[2.5,97.5])
print(f"  N={len(s)}  country-waves={s.cw.nunique()}  mean ideal={np.average(s.D017c,weights=s.wt):.2f}  mean realised={np.average(s.X011c,weights=s.wt):.2f}")
print(f"  elasticity of IDEAL fertility    = {est[0]:+.4f} ({se[0]:.4f})")
print(f"  elasticity of REALISED fertility = {est[1]:+.4f} ({se[1]:.4f})")
print(f"  ratio ideal/realised             = {est[2]:+.3f}  [{ci[0]:+.3f},{ci[1]:+.3f}]   (tastes only -> ~1; costs only -> 0)")
print(f"  slope of gap (ideal - realised)  = {est[3]:+.4f} ({se[3]:.4f}) children per log point (costs only -> negative: gap widens with density)")
