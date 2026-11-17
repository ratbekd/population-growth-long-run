import pandas as pd, numpy as np
src=open('cex_est2.py').read().split("w=c[c.wage.between(3,200)")[0]
src=src.replace("keep=['c']+[k for k in X.columns","X=X.join(g[['lnw']]) if lev>=2 else X\n    keep=['c']+[k for k in X.columns")
exec(src)
_X0=X_of
def X_of(g,lev):
    X=_X0(g,1 if lev==3 else (0 if lev==2 else lev))
    return X.join(g[['lnw']]) if lev>=2 else X
w=c[c.wage.between(3,200)&(c.NO_EARNR<=2)&(c['T']>c.E)&(c['T']>0)].copy()
w['lnC']=np.log(w['T']-w.E); w['lnw']=np.log(w.wage); w['lnCw']=w.lnC-w.lnw; w['lnhrs']=np.log(w.hrs)
w['Qw99']=w.Qw.clip(upper=w.Qw.quantile(.99))
S=[('fertility: children <18 (PPML)','PERSLT18','ppml'),('ln(consumption/wage)','lnCw','ols'),('ln annual hours','lnhrs','ols'),('education/wage, winsorised (PPML)','Qw99','ppml')]
run2(w,S,"CONTROLLING FOR LN WAGE",lev=2)
run2(w,S,"CONTROLLING FOR LN WAGE AND EDUCATION",lev=3)
