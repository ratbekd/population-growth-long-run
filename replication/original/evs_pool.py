import pandas as pd, numpy as np, pycountry, pycountry_convert as pc, warnings; warnings.filterwarnings("ignore")
exec(open('taste_tests.py').read().split("# ================= TEST 1")[0])
# ---- harmonise WVS and EVS ----
w=pd.read_pickle("/home/claude/work/ivs_clean.pkl"); w['study']='W'; w['wv']=w.S002
e=pd.read_pickle("evs_sub.pkl"); e['study']='E'; e['wv']=e.evsw
cols=['study','wv','S003','S017','S020','D017','F028','X001','X003','X007','X011','X049']
d=pd.concat([w[cols],e[cols]],ignore_index=True)
for c in ['D017','X011','X049','F028','X007']: d[c]=d[c].where(d[c]>=0)
d['year']=d.S020.astype(int)
dup=set(zip(e.S003,e.S020))&set(zip(w.S003,w.S020)); d=d[~((d.study=='E')&pd.Series(list(zip(d.S003,d.S020))).isin(dup).values)]
d['age']=d.X003.where(d.X003.between(15,99)); d['sex']=d.X001.where(d.X001.isin([1,2]))
d['wt']=d.S017.where(d.S017>0,1.0); d['cw']=d.study+"_"+d.S003.astype(int).astype(str)+"_"+d.wv.astype(int).astype(str)
d['period']=pd.cut(d.year,[1980,1986,1994,2001,2005,2010,2017,2030],labels=False)
def iso3(n):
    n=int(n)
    if n==158: return 'TWN'
    try: return pycountry.countries.get(numeric=f"{n:03d}").alpha_3
    except Exception: return {197:'CYP',909:None}.get(n)
d['iso3']=d.S003.map({k:iso3(k) for k in d.S003.dropna().unique()})
def cont(a3):
    try: return pc.country_alpha2_to_continent_code(pc.country_alpha3_to_country_alpha2(a3))
    except Exception: return {'TWN':'AS','XKX':'EU','TLS':'AS'}.get(a3,'OT')
d['cont']=d.iso3.map({k:cont(k) for k in d.iso3.dropna().unique()})
p=pd.read_csv("/home/claude/work2/wpp.csv")
for k in ['year','pop_k','density','tfr']: p[k]=pd.to_numeric(p[k],errors='coerce')
gd=pd.read_csv("/home/claude/data/gdp/gdp.csv").rename(columns={'Country Code':'iso3','Year':'year','Value':'gdp'})
d=d.merge(p[['iso3','year','density','pop_k']],on=['iso3','year'],how='left').merge(gd[['iso3','year','gdp']],on=['iso3','year'],how='left')
d['nhab']=-np.log(d.density); d['lgdppc']=np.log(d.gdp/(d.pop_k*1e3))
mid={1:1000,2:3162,3:7071,4:14142,5:31623,6:70711,7:223607,8:1000000}; d['lhab']=-np.log(d.X049.map(mid)/1e3)
print("unmatched iso:",sorted(d.loc[d.iso3.isna(),'S003'].unique()),"  missing GDP country-years:",d.loc[d.D017.notna()&d.gdp.isna()].drop_duplicates(['iso3','year'])[['iso3','year']].values.tolist()[:12])
def Xn(g,extra,fe):
    X=pd.DataFrame({'c':1.0,'age':g.age,'age2':g.age**2/100,'fem':(g.sex==2)*1.0,'rel':g.F028},index=g.index)
    for k in extra: X[k]=g[k]
    X=X.join(pd.get_dummies(g.X007.astype(int),prefix='m',drop_first=True).astype(float))
    for f in fe: X=X.join(pd.get_dummies(g[f].astype(str),prefix=f,drop_first=True).astype(float))
    return X[['c']+[k for k in X.columns if k!='c' and X[k].std()>0]]
def nspec(g,y,extra,fe):
    X=Xn(g,extra,fe); b=ppml(g[y].values,X.values,g.wt.values); cols=list(X.columns); return [b[cols.index(k)] for k in extra]
def show(lab,g,y,extra,fe,cl='iso3',B=150):
    g=g.dropna(subset=[y]+extra+['age','sex','F028','X007']).reset_index(drop=True)
    est,bs=boot(lambda gg: nspec(gg,y,extra,fe),g,cl,B=B); se=bs.std(0)
    print(f"  {lab:46s} N={len(g):6d} countries={g.iso3.nunique():3d} cw={g.cw.nunique():3d} | "+"  ".join(f"{k}={v:+.4f}({s:.4f}) p={2*min((bs[:,i]<=0).mean(),(bs[:,i]>=0).mean()):.3f}" for i,(k,v,s) in enumerate(zip(extra,est,se))),flush=True)
    return est,se
s=d[d.D017.notna()&(d.D017<=8)]
print("\nPANEL C pooled WVS + EVS: IDEAL fertility on national habitat")
show("WVS only, wave FE (replication)",s[s.study=='W'].assign(wvf=lambda x:x.study+x.wv.astype(str)),'D017',['nhab','lgdppc'],['wvf'])
s=s.assign(wvf=s.study+s.wv.astype(str))
show("pooled, national + GDP, study-wave FE",s,'D017',['nhab','lgdppc'],['wvf'])
show("pooled, national + GDP, wave + continent FE",s,'D017',['nhab','lgdppc'],['wvf','cont'])
show("pooled, national + local + GDP, study-wave FE",s[s.lhab.notna()],'D017',['nhab','lhab','lgdppc'],['wvf'])
show("pooled, national + local + GDP, + continent FE",s[s.lhab.notna()],'D017',['nhab','lhab','lgdppc'],['wvf','cont'])
eu=s[s.cont=='EU']
show("Europe only, national + GDP, study-wave FE",eu,'D017',['nhab','lgdppc'],['wvf'])
show("Europe only, national + local + GDP",eu[eu.lhab.notna()],'D017',['nhab','lhab','lgdppc'],['wvf'])
show("pooled, national + GDP, country + wave FE (within)",s,'D017',['nhab','lgdppc'],['wvf','iso3'])
print("\nPANEL A pooled: ideal vs realised, ages 45-70, country-wave FE")
a=s[s.age.between(45,70)&s.X011.notna()&s.lhab.notna()].copy(); a['X011c']=a.X011.clip(upper=5); a['D017c']=a.D017.clip(upper=5)
def afn(g):
    X=pd.DataFrame({'c':1.0,'hab':g.lhab,'age':g.age,'age2':g.age**2/100,'fem':(g.sex==2)*1.0,'rel':g.F028},index=g.index)
    X=X.join(pd.get_dummies(g.X007.astype(int),prefix='m',drop_first=True).astype(float)).join(pd.get_dummies(g.cw,prefix='cw',drop_first=True).astype(float))
    X=X[['c']+[k for k in X.columns if k!='c' and X[k].std()>0]]; j=list(X.columns).index('hab')
    bi=ppml(g.D017c.values,X.values,g.wt.values)[j]; br=ppml(g.X011c.values,X.values,g.wt.values)[j]; return [bi,br,bi/br]
for lab,g in [("WVS only",a[a.study=='W']),("pooled WVS + EVS",a),("EVS only",a[a.study=='E'])]:
    g=g.dropna(subset=['F028','X007','sex']).reset_index(drop=True)
    est,bs=boot(afn,g,'cw',B=120); se=bs.std(0); ci=np.percentile(bs[:,2],[2.5,97.5])
    print(f"  {lab:18s} N={len(g):6d} cw={g.cw.nunique():3d} | ideal {est[0]:+.4f}({se[0]:.4f}) realised {est[1]:+.4f}({se[1]:.4f}) tau {est[2]:+.3f} [{ci[0]:+.3f},{ci[1]:+.3f}]",flush=True)
