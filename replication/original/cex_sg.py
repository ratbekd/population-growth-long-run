import pandas as pd, numpy as np
src=open('cex_est2.py').read().split("w=c[c.wage.between(3,200)")[0]; exec(src)
w=c[c.wage.between(3,200)&(c.NO_EARNR<=2)&(c['T']>c.E)].copy()
w['C']=4*(w['T']-w.E)/w.adults.round(); w['Ed']=4*w.E.clip(upper=w.E.quantile(.99))/w.adults.round(); w['W']=w.wage*2000   # per adult, annual
w['habc']=w.hab-np.average(w.hab,weights=w.FINLWT21)
def fit(g,y,ii=None):
    X=X_of(g,0).drop(columns=['hab']); X['hab']=g.habc; X['W']=g.W/1e4; X['Wh']=g.W/1e4*g.habc
    Xv=X.values; wt=g.FINLWT21.values; yv=g[y].values
    if ii is not None: Xv,wt,yv=Xv[ii],wt[ii],yv[ii]
    b=ols(yv,Xv,wt); k=list(X.columns); return b[k.index('W')],b[k.index('Wh')]
def run(g,lab,B=200):
    g=g.dropna(subset=['C','W','habc','AGE_REF','SEX_REF','FINLWT21'])
    Xn=X_of(g,0); jn=list(Xn.columns).index('hab'); bn0=ppml(g.PERSLT18.values,Xn.values,g.FINLWT21.values); pn=bn0[jn]
    est={y:fit(g,y) for y in('C','Ed')}
    codes,u=pd.factorize(g.clu); idx=[np.where(codes==k)[0] for k in range(len(u))]; rng=np.random.default_rng(4); bs=[]
    for _ in range(B):
        ii=np.concatenate([idx[k] for k in rng.integers(0,len(u),len(u))])
        a=ppml(g.PERSLT18.values[ii],Xn.values[ii],g.FINLWT21.values[ii],it=10,b0=bn0)[jn]
        r=[a]
        for y in('C','Ed'):
            b,dl=fit(g,y,ii); r+=[b,dl/b]
        bs.append(r)
    bs=np.array(bs)
    print(f"{lab}: N={len(g)}  pi_n={pn:+.4f} ({bs[:,0].std():.4f})")
    for k,y in enumerate(('C','Ed')):
        b,dl=est[y]; e=dl/b; r=bs[:,2+2*k]/bs[:,0]; 
        print(f"   {y:3s} slope in W/1e4 = {b:8.1f} ({bs[:,1+2*k].std():.1f}); d ln slope/d hab = {e:+.4f} ({bs[:,2+2*k].std():.4f}); ratio to pi_n = {e/pn:+.3f} [{np.percentile(r,2.5):+.3f},{np.percentile(r,97.5):+.3f}]",flush=True)
print("Stone-Geary slope test: taste-only predicts ratio in [-0.160,-0.126]; cost channels predict 0")
run(w,"headline"); run(w[w.MARITAL1==1],"married"); run(w[w.yr<=2012],"2000-2012"); run(w[w.yr>2012],"2013-2026"); run(w[w.wage.between(8,80)],"wage 8-80")
