"""Reproduce Table 6; audit legacy bootstrap and refit without coefficient trimming.

The sample and outcome construction are deliberately held fixed in this solver audit.
They are documented separately and are not certified as a structural identification design.
"""
import os, sys
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
LOCAL=Path(r'C:\Users\ratbekd\.codex\visualizations\2026\11\11\01a1c8fa-6b21-7f31-a0ee-8ee897cb4798')
sys.path.insert(0,str(LOCAL/'cex_deps'))
import numpy as np, pandas as pd, scipy, scipy.sparse as sp, scipy.linalg as la
import statsmodels.api as sm
import argparse, json, time, platform, warnings
HERE=Path(__file__).resolve().parent; OUT=HERE/'results'; CACHE=LOCAL/'cex_replication'
OUT.mkdir(exist_ok=True)

class FitFailure(RuntimeError): pass

def solve_system(H,rhs):
 with warnings.catch_warnings():
  warnings.simplefilter('error',la.LinAlgWarning)
  try: return la.solve(H,rhs,assume_a='pos')
  except (la.LinAlgError,la.LinAlgWarning):
   return la.lstsq(H,rhs,cond=1e-12,lapack_driver='gelsd')[0]

class Model:
 def __init__(self,X,y,w,kind):
  self.kind=kind; self.scale=np.sqrt(np.mean(X*X,axis=0)); self.scale[0]=1
  self.X=sp.csr_matrix(X/self.scale); self.y=y.astype(float)
  self.yscale=max(np.average(y,weights=w),1e-12) if kind=='ppml' else 1.
  self.Y=self.y/self.yscale; self.wbase=w/np.mean(w)
 def fit(self,mult=None,start=None,maxiter=200):
  w=self.wbase if mult is None else self.wbase*mult
  X=self.X; Y=self.Y; denom=w.sum()*max(1,np.average(np.abs(Y),weights=w))
  b=np.zeros(X.shape[1]) if start is None else start.copy()
  gram=(X.T @ X.multiply(w[:,None])).toarray()
  ev,vec=la.eigh(gram)
  null=vec[:,ev<ev[-1]*1e-12]
  if null.shape[1] and np.linalg.norm(null[1,:])>1e-6:
   raise FitFailure('habitat coefficient not identified in resampled design')
  if self.kind=='ols':
   H=(X.T @ X.multiply(w[:,None])).toarray(); rhs=np.asarray(X.T@(w*Y)).ravel()
   try: b=solve_system(H,rhs)
   except la.LinAlgError: raise FitFailure('rank-deficient weighted design')
   score=float(np.max(np.abs(X.T@(w*(Y-X@b))))/denom)
   return b,{'converged':score<1e-8,'iterations':1,'score':score,'backtracks':0}
  backtracks=0
  for it in range(maxiter):
   eta=np.asarray(X@b); mu=np.exp(eta)
   if not np.isfinite(mu).all(): raise FitFailure('nonfinite fitted mean')
   grad=np.asarray(X.T@(w*(Y-mu))).ravel(); score=float(np.max(np.abs(grad))/denom)
   if score<1e-10:
    return b,{'converged':True,'iterations':it,'score':score,'backtracks':backtracks}
   H=(X.T @ X.multiply((w*mu)[:,None])).toarray()
   try: step=solve_system(H,grad)
   except la.LinAlgError: raise FitFailure('singular Hessian / separation candidate')
   objective=float(np.dot(w,mu-Y*eta)); rate=1.
   for bt in range(45):
    bn=b+rate*step; en=np.asarray(X@bn)
    if en.max()<700:
     objn=float(np.dot(w,np.exp(en)-Y*en))
     if np.isfinite(objn) and objn<=objective+1e-12*max(1,abs(objective)): break
    rate*=.5
   else: raise FitFailure('line search failed')
   backtracks+=bt; b=bn
  raise FitFailure(f'iteration limit, last score={score}')
 def coefficients(self,b):
  c=b/self.scale
  if self.kind=='ppml': c[0]+=np.log(self.yscale)
  return c

def legacy_fit(X,y,w,kind,start=None):
 X=sp.csr_matrix(X)
 if kind=='ols': return np.linalg.solve((X.T@X.multiply(w[:,None])).toarray()+1e-8*np.eye(X.shape[1]),X.T@(w*y))
 b=np.zeros(X.shape[1]) if start is None else start.copy()
 if start is None: b[0]=np.log(max(np.average(y,weights=w),1e-9))
 for _ in range(80 if start is None else 10):
  eta=X@b; mu=np.exp(np.clip(eta,-30,30)); W=w*mu; z=eta+(y-mu)/mu
  bn=np.linalg.solve((X.T@X.multiply(W[:,None])).toarray()+1e-6*np.eye(X.shape[1]),X.T@(W*z))
  if np.max(np.abs(bn-b))<1e-8: b=bn; break
  b=bn
 return b

def design(g,wage=False):
 X=pd.DataFrame({'c':1.,'hab':g.hab,'age':g.AGE_REF,'age2':g.AGE_REF**2/100,'fem':(g.SEX_REF==2)*1.},index=g.index)
 for v in ['MARITAL1','REF_RACE','REGION','yr']:
  X=X.join(pd.get_dummies(g[v].fillna(-1).astype(int),prefix=v,drop_first=True).astype(float))
 X=X[['c']+[k for k in X.columns if k!='c' and X[k].std()>0]]
 if wage: X=X.join(g[['lnw']])
 return X

def samples():
 allc=pd.read_pickle(CACHE/'cex_cu.pkl'); flow=[]
 def record(label,g): flow.append({'step':label,'N':len(g),'clusters':g.clu.nunique()})
 record('CUs after original location requirements',allc)
 c=allc[allc.AGE_REF.between(25,44)].copy(); record('Reference age 25-44',c)
 c=c[c.adults.round().between(1,2)].copy(); record('Rounded mean adult count 1 or 2 (panel C base)',c)
 c['Qw']=c.E/c.wage; c['Cw']=(c['T']-c.E)/c.wage
 w=c[c.wage.between(3,200)].copy(); record('Hourly wage in [3,200]',w)
 w=w[w.NO_EARNR<=2].copy(); record('Mean earners <=2',w)
 w=w[(w['T']>w.E)&(w['T']>0)].copy(); record('Positive noneducation and total expenditure',w)
 w['lnC']=np.log(w['T']-w.E); w['lnw']=np.log(w.wage); w['lnCw']=w.lnC-w.lnw
 w['lnhrs']=np.log(w.hrs); w['lnNSw']=np.log((w['T']-w.E-w.SH).clip(lower=1))-w.lnw
 w['Qw99']=w.Qw.clip(upper=w.Qw.quantile(.99)); w['anyE']=(w.E>0)*1.
 c['hrs_']=c.hrs
 pd.DataFrame(flow).to_csv(OUT/'sample_flow.csv',index=False)
 (OUT/'construction_details.json').write_text(json.dumps({'winsorized_variable':'education / hourly wage','Qw_cap':float(w.Qw.quantile(.99)),'cached_rows':len(allc)},indent=2))
 A=[('children','PERSLT18','ppml'),('consumption_share','lnCw','ols'),('consumption','lnC','ols'),('wage','lnw','ols'),('hours','lnhrs','ols'),('nonshelter_share','lnNSw','ols'),('education_share','Qw99','ppml'),('any_education','anyE','ppml')]
 B=[A[0],A[1],A[4],A[6]]
 C=[A[0],('hours','hrs_','ppml')]
 return {'A':(w,A,False),'B':(w,B,True),'C':(c,C,False)}

def run(panel,g,specs,wage,B,legacy_B):
 cols=[s[1] for s in specs]
 g=g.dropna(subset=cols+['hab','AGE_REF','SEX_REF','FINLWT21'])
 g=g[np.isfinite(g[cols]).all(axis=1)].copy()
 assert (g.FINLWT21>0).all()
 Xdf=design(g,wage); X=Xdf.to_numpy(); j=list(Xdf).index('hab'); w=g.FINLWT21.to_numpy()
 codes,clusters=pd.factorize(g.clu); C=len(clusters)
 models=[]; bases=[]; oldbases=[]; rows=[]; diag=[]
 for name,col,kind in specs:
  y=g[col].to_numpy(); model=Model(X,y,w,kind); b,d=model.fit(); coef=model.coefficients(b)
  old=legacy_fit(X,y,w,kind)
  if kind=='ppml':
   check=sm.GLM(y,X,family=sm.families.Poisson(),freq_weights=w/w.mean()).fit(maxiter=200,tol=1e-10)
   assert check.converged, name
  else: check=sm.WLS(y,X,weights=w).fit()
  err=float(np.max(np.abs(coef-check.params)))
  assert err<1e-5,(name,err)
  rows.append({'panel':panel,'outcome':name,'N':len(g),'clusters':C,'legacy_coefficient':float(old[j]),'repaired_coefficient':float(coef[j]),'statsmodels_max_abs_difference':err})
  models.append(model); bases.append(b); oldbases.append(old); diag.append({'outcome':name,**d})
  print('POINT',panel,name,len(g),C,round(coef[j],7),'legacy',round(old[j],7),flush=True)
 pd.DataFrame(rows).to_csv(OUT/f'point_{panel}.csv',index=False)
 draws=np.full((B,len(specs)),np.nan); olddraws=np.full((min(B,legacy_B),len(specs)),np.nan)
 logs=[]; rng=np.random.default_rng(3); t=time.monotonic()
 for r in range(B):
  pick=rng.integers(0,C,C); counts=np.bincount(pick,minlength=C); mult=counts[codes]
  log={'replication':r,'panel':panel,'failures':[],'max_iterations':0,'max_score':0.,'retry_count':0}
  for k,model in enumerate(models):
   try:
    try: b,d=model.fit(mult,bases[k])
    except FitFailure:
     log['retry_count']+=1; b,d=model.fit(mult,None,maxiter=500)
    if not d['converged']: raise FitFailure('score tolerance failed')
    draws[r,k]=model.coefficients(b)[j]
    log['max_iterations']=max(log['max_iterations'],d['iterations']); log['max_score']=max(log['max_score'],d['score'])
   except Exception as e: log['failures'].append({'outcome':specs[k][0],'reason':str(e)})
   if r<legacy_B:
    try: olddraws[r,k]=legacy_fit(X,model.y,w*mult,model.kind,oldbases[k])[j]
    except Exception: pass
  logs.append(log)
  if (r+1)%25==0:
   print('BOOT',panel,r+1,'/',B,'elapsed',round(time.monotonic()-t,1),'failed',sum(bool(q['failures']) for q in logs),flush=True)
  if (r+1)%100==0:
   np.save(OUT/f'draws_{panel}.npy',draws)
 np.save(OUT/f'draws_{panel}.npy',draws); np.save(OUT/f'legacy_draws_{panel}.npy',olddraws)
 (OUT/f'diagnostics_{panel}.json').write_text(json.dumps({'point':diag,'bootstrap':logs,'columns':list(Xdf)},indent=2))
 estimates=np.array([r['repaired_coefficient'] for r in rows]); oldest=np.array([r['legacy_coefficient'] for r in rows])
 oldkeep=np.isfinite(olddraws).all(axis=1)&(np.abs(olddraws-oldest)<.5).all(axis=1)
 for k,row in enumerate(rows):
  valid=np.isfinite(draws[:,k]); row['bootstrap_attempted']=B; row['bootstrap_valid']=int(valid.sum())
  row['bootstrap_se']=float(np.std(draws[valid,k],ddof=1)) if valid.any() else None
  row['legacy_valid_before_trim']=int(np.isfinite(olddraws).all(axis=1).sum()); row['legacy_retained']=int(oldkeep.sum())
  row['legacy_bootstrap_se']=float(np.std(olddraws[oldkeep,k],ddof=0)) if oldkeep.any() else None
  if k:
   valid=valid&np.isfinite(draws[:,0])&(draws[:,0]!=0)
   ratios=draws[valid,k]/draws[valid,0]
   row['ratio']=float(estimates[k]/estimates[0]); row['ratio_lo'],row['ratio_hi']=np.quantile(ratios,[.025,.975]).tolist()
   if oldkeep.any(): row['legacy_ratio_lo'],row['legacy_ratio_hi']=np.quantile(olddraws[oldkeep,k]/olddraws[oldkeep,0],[.025,.975]).tolist()
  row['denominator_nonpositive_draws']=int((draws[:,0]<=0).sum())
 pd.DataFrame(rows).to_csv(OUT/f'results_{panel}.csv',index=False)
 print('COMPLETE',panel,'legacy retained',oldkeep.sum(),'/',len(olddraws),'failed new',sum(bool(q['failures']) for q in logs),flush=True)

if __name__=='__main__':
 parser=argparse.ArgumentParser(); parser.add_argument('--bootstrap',type=int,default=999); parser.add_argument('--legacy-bootstrap',type=int,default=200); parser.add_argument('--panels',default='ABC'); args=parser.parse_args()
 (OUT/'environment.json').write_text(json.dumps({'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'statsmodels':sm.__version__},indent=2))
 data=samples()
 for panel in args.panels: run(panel,*data[panel],args.bootstrap,args.legacy_bootstrap)
