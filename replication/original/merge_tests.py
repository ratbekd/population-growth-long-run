import pandas as pd, numpy as np, pycountry, warnings; warnings.filterwarnings("ignore")
exec(open('taste_tests.py').read().split("# ================= TEST 1")[0])
# ---------- robustness on the two direct tests ----------
print("ROBUSTNESS")
a=pd.read_pickle("gss_sub.pkl"); e=pd.read_pickle("gss_extra.pkl")[['year','id','sibs','paeduc','maeduc','mobile16']]
d=a.merge(e,on=['year','id'],how='left'); d=d[(d.year>=1975)&d['size'].notna()&d.vpsu.notna()&d.res16.between(1,6)].copy()
d['hab']=-np.log(d['size'].clip(lower=0.5)); d['hab16']=-np.log(d.res16.map({2:0.5,1:2.0,3:20.0,4:120.0,5:250.0,6:1000.0}))
d['psu']=d.year.astype(int).astype(str)+"_"+d.vstrat.astype(int).astype(str)+"_"+d.vpsu.astype(int).astype(str); d['wt']=d.wtssps
d['ideal']=d.chldidel.where(d.chldidel<=7); d['lsibs']=np.log1p(d.sibs)
def X_of(g,extra):
    X=pd.DataFrame({'c':1.0,'age':g.age,'age2':g.age**2/100,'fem':(g.sex==2)*1.0,'att':g.attend},index=g.index)
    for k in extra: X[k]=g[k]
    for v in ['marital','race','region','year']: X=X.join(pd.get_dummies(g[v].astype(int),prefix=v,drop_first=True).astype(float))
    return X[['c']+[k for k in X.columns if k!='c' and X[k].std()>0]]
def spec(g,y,extra):
    X=X_of(g,extra); b=ppml(g[y].values,X.values,g.wt.values); cols=list(X.columns); return [b[cols.index(k)] for k in extra]
q=d[d.age.between(25,64)].dropna(subset=['ideal','hab','hab16','age','sex','attend','marital','race','region','lsibs','paeduc','maeduc'])
for lab,ex in [("ideal, both, same sample",['hab','hab16']),("ideal + siblings + parents' educ",['hab','hab16','lsibs','paeduc','maeduc'])]:
    est,bs=boot(lambda gg: spec(gg,'ideal',ex),q,'psu',B=120); se=bs.std(0)
    print(f"  GSS {lab:34s} N={len(q)} hab={est[0]:+.4f}({se[0]:.4f}) hab16={est[1]:+.4f}({se[1]:.4f})",flush=True)

# ---------- WVS x WPP x World Bank ----------
w=pd.read_pickle("/home/claude/work/ivs_clean.pkl")
for c in ['D017','X011','X049','F028','X007','X047_WVS','X025']: w[c]=w[c].where(w[c]>=0)
w['age']=w.X003.where(w.X003.between(15,99)); w['sex']=w.X001.where(w.X001.isin([1,2]))
w['cw']=w.S003.astype(str)+"_"+w.S002.astype(str); w['wt']=w.S017.where(w.S017>0,1.0)
def iso3(n):
    n=int(n)
    if n==158: return 'TWN'
    try: return pycountry.countries.get(numeric=f"{n:03d}").alpha_3
    except Exception: return None
w['iso3']=w.S003.map({k:iso3(k) for k in w.S003.dropna().unique()}); w['year']=w.S020.astype(int)
p=pd.read_csv("/home/claude/work2/wpp.csv"); 
for k in ['year','pop_k','density','tfr']: p[k]=pd.to_numeric(p[k],errors='coerce')
gd=pd.read_csv("/home/claude/data/gdp/gdp.csv").rename(columns={'Country Code':'iso3','Year':'year','Value':'gdp'})
w=w.merge(p[['iso3','year','density','tfr','pop_k']],on=['iso3','year'],how='left').merge(gd[['iso3','year','gdp']],on=['iso3','year'],how='left')
w['nhab']=-np.log(w.density); w['lgdppc']=np.log(w.gdp/(w.pop_k*1e3))
mid={1:1000,2:3162,3:7071,4:14142,5:31623,6:70711,7:223607,8:1000000}; w['lhab']=-np.log(w.X049.map(mid)/1e3)
print(f"\nMERGE: WVS rows {len(w)}; with WPP density {w.density.notna().mean():.3f}; with GDP {w.gdp.notna().mean():.3f}; unmatched codes: {sorted(w.loc[w.iso3.isna(),'S003'].unique())[:10]}")
def Xn(g,extra,fe):
    X=pd.DataFrame({'c':1.0,'age':g.age,'age2':g.age**2/100,'fem':(g.sex==2)*1.0,'rel':g.F028},index=g.index)
    for k in extra: X[k]=g[k]
    X=X.join(pd.get_dummies(g.X007.astype(int),prefix='m',drop_first=True).astype(float))
    for f in fe: X=X.join(pd.get_dummies(g[f].astype(str),prefix=f,drop_first=True).astype(float))
    return X[['c']+[k for k in X.columns if k!='c' and X[k].std()>0]]
def nspec(g,y,extra,fe):
    X=Xn(g,extra,fe); b=ppml(g[y].values,X.values,g.wt.values); cols=list(X.columns); return [b[cols.index(k)] for k in extra]
def show(lab,g,y,extra,fe,B=120):
    g=g.dropna(subset=[y]+extra+['age','sex','F028','X007']).reset_index(drop=True)
    est,bs=boot(lambda gg: nspec(gg,y,extra,fe),g,'iso3',B=B); se=bs.std(0)
    print(f"  {lab:50s} N={len(g):6d} countries={g.iso3.nunique():3d} cw={g.cw.nunique():3d} | "+"  ".join(f"{k}={v:+.4f}({s:.4f})" for k,v,s in zip(extra,est,se)),flush=True)
print("\nIDEAL fertility (D017, waves 1-4) on NATIONAL habitat per capita (WPP), clustered by country")
s=w[w.S002<=4]; s=s[s.D017<=8]
show("national hab + ln GDPpc, wave FE",s,'D017',['nhab','lgdppc'],['S002'])
show("national + local hab + ln GDPpc, wave FE",s[s.lhab.notna()],'D017',['nhab','lhab','lgdppc'],['S002'])
show("national hab + ln GDPpc, country + wave FE",s,'D017',['nhab','lgdppc'],['iso3','S002'])
print("\nREALISED completed fertility (X011, 45-70, waves 1-7) on national habitat")
r=w[w.age.between(45,70)]
show("national hab + ln GDPpc, wave FE",r,'X011',['nhab','lgdppc'],['S002'])
show("national + local hab + ln GDPpc, wave FE",r[r.lhab.notna()],'X011',['nhab','lhab','lgdppc'],['S002'])
show("national hab + ln GDPpc, country + wave FE",r,'X011',['nhab','lgdppc'],['iso3','S002'])
print("\nIDEAL fertility by cohort (rationalisation check), WVS waves 2-4, local hab, country-wave FE")
for lo,hi in [(18,34),(35,44),(45,70)]:
    g=w[w.S002.isin([2,3,4])&w.age.between(lo,hi)&(w.D017<=8)]
    show(f"ages {lo}-{hi}",g[g.lhab.notna()],'D017',['lhab'],['cw'],B=80)
