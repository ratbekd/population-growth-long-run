import pandas as pd, numpy as np
d=pd.read_pickle('cex_raw.pkl')
d['CU']=d.file.str[:2]+"_"+(d.NEWID//10).astype('int64').astype(str)   # NEWID = CU id + interview number
intno=(d.NEWID%10).astype(int); yq=d.QINTRVYR*4+((d.QINTRVMO-1)//3)
d['CU']=(d.NEWID//10).astype('int64').astype(str)+"_"+(yq-(intno-1)).astype(int).astype(str)
d['E']=d.EDUCAPQ.fillna(0)+d.EDUCACQ.fillna(0); d['T']=d.TOTEXPPQ.fillna(0)+d.TOTEXPCQ.fillna(0)
d['SH']=d.SHELTPQ.fillna(0)+d.SHELTCQ.fillna(0)
d['sal']=d.FSALARYM.where(d.FSALARYM.notna(),d.FSALARYX)
h1=(d.INC_HRS1.fillna(0)*d.INCWEEK1.fillna(0)); h2=(d.INC_HRS2.fillna(0)*d.INCWEEK2.fillna(0))
d['hrs']=h1+h2; d['nwork']=(h1>0).astype(int)+(h2>0).astype(int)
d['adults']=d.FAM_SIZE-d.PERSLT18
d['yr']=d.QINTRVYR
num=['E','T','SH','sal','hrs','nwork','adults','PERSLT18','FAM_SIZE','FINLWT21','AGE_REF','NO_EARNR','yr']
cat=['SEX_REF','MARITAL1','EDUC_REF','REF_RACE','REGION','BLS_URBN','POPSIZE','STATE','FAM_TYPE','CUTENURE']
g=d.sort_values(['CU','NEWID']).groupby('CU')
c=g[num].mean().join(g[cat].first()).join(g.size().rename('nint'))
c['yr']=c.yr.round().astype(int)
old={1:8e6,2:2.19e6,3:6.3e5,4:2.0e5,5:4.0e4}; new={1:9e6,2:2.24e6,3:7.1e5,4:2.2e5,5:3.5e4}
c['pop']=np.where(c.yr>=2015,c.POPSIZE.map(new),c.POPSIZE.map(old)); c.loc[c.BLS_URBN==2,'pop']=5e3
c['hab']=-np.log(c['pop']/1e3)
c['wage']=np.where((c.hrs>0)&(c.sal>0),c.sal/c.hrs,np.nan)
c=c[c.POPSIZE.notna()&c.BLS_URBN.notna()].copy()
c['clu']=c.STATE.fillna('NA').astype(str)+"_"+c.POPSIZE.astype(int).astype(str)+"_"+c.BLS_URBN.astype(int).astype(str)
c.to_pickle('cex_cu.pkl')
print("CUs:",len(c)," mean interviews:",c.nint.mean().round(2)," clusters:",c.clu.nunique())
s=c[c.AGE_REF.between(25,44)]
print("age 25-44:",len(s)," with wage:",s.wage.notna().sum()," E>0 share:",(s.E>0).mean().round(3))
print(s[['E','T','sal','hrs','wage','PERSLT18','adults']].describe(percentiles=[.01,.5,.99]).round(1).T)
print("sd hab:",s.hab.std().round(3)," by pop:",s.groupby('pop').size().to_dict())
