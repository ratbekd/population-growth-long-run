import pandas as pd, numpy as np, glob, os
keep=['NEWID','QINTRVMO','QINTRVYR','FINLWT21','EDUCAPQ','EDUCACQ','TOTEXPPQ','TOTEXPCQ','FSALARYX','FSALARYM','FINCBTAX','FINCBTXM',
 'FAM_SIZE','PERSLT18','AGE_REF','AGE2','SEX_REF','MARITAL1','EDUC_REF','REF_RACE','REGION','BLS_URBN','POPSIZE','SMSASTAT','STATE','PSU',
 'INC_HRS1','INC_HRS2','INCWEEK1','INCWEEK2','FAM_TYPE','NO_EARNR','CUTENURE','SHELTPQ','SHELTCQ','ROOMSQ']
fs=sorted(glob.glob('/home/claude/data/cex/fmli*.csv')); out=[]
for f in fs:
    b=os.path.basename(f)[4:-4]          # e.g. 001x, 002, 191
    yy,q,x=b[:2],b[2],b[3:]
    if q=='1' and x=='' and os.path.exists(f'/home/claude/data/cex/fmli{yy}1x.csv'): continue   # prefer the 'x' first-quarter files
    d=pd.read_csv(f,encoding='latin-1',low_memory=False,dtype=str)
    d.columns=[c.strip().upper() for c in d.columns]
    d=d[[c for c in keep if c in d.columns]].copy(); d['file']=b
    out.append(d)
d=pd.concat(out,ignore_index=True)
for c in d.columns:
    if c not in('file','STATE','PSU'): d[c]=pd.to_numeric(d[c],errors='coerce')
print(d.shape, d.QINTRVYR.min(), d.QINTRVYR.max()); d.to_pickle('cex_raw.pkl')
print(d.groupby('QINTRVYR').size().to_dict())
print("POPSIZE x period:"); d['per']=np.where(d.QINTRVYR<=2005,'00-05',np.where(d.QINTRVYR<=2014,'06-14','15-19'))
print(pd.crosstab(d.POPSIZE,d.per,normalize='columns').round(3))
print("BLS_URBN:",d.BLS_URBN.value_counts(dropna=False).to_dict()); print("SMSASTAT:",d.SMSASTAT.value_counts(dropna=False).to_dict())
print(pd.crosstab(d.POPSIZE,d.BLS_URBN))
print("educ>0 share:",((d.EDUCAPQ.fillna(0)+d.EDUCACQ.fillna(0))>0).mean().round(3))
print("STATE missing:",d.STATE.isna().mean().round(3),"n states",d.STATE.nunique(),"PSU nonmissing:",d.PSU.notna().mean().round(3) if 'PSU' in d else None)
