import pandas as pd, numpy as np, statsmodels.api as sm
gam,eta,phi,th=0.26,0.07,0.138,0.83; G=gam*(eta+th); a0=0.39
A=lambda n:(phi*n+G)/(1-phi*n)
w=pd.read_csv("wpp.csv"); 
for k in ['year','pop_k','density','tfr']: w[k]=pd.to_numeric(w[k],errors='coerce')
w=w.dropna(subset=['tfr','density','year'])
print("WPP rows",len(w),"countries",w.iso3.nunique(),"years",int(w.year.min()),int(w.year.max()))
t=w[w.year==2023].set_index('iso3')
D0=float(t.loc['FRA','density']); print(f"\nFrance 2023: TFR={t.loc['FRA','tfr']:.2f} density={D0:.1f}")
print(f"{'':12}{'TFR':>6}{'alpha':>8}", "".join(f"{'eps@'+str(x):>10}{'Dmodel':>8}" for x in (0.146,0.123)), f"{'D obs':>8}{'ratio.146':>10}{'ratio.123':>10}")
for iso,nm in [('KOR','Korea'),('TWN','Taiwan'),('JPN','Japan'),('ITA','Italy'),('FRA','France'),('SGP','Singapore')]:
    tfr,d=t.loc[iso,'tfr'],t.loc[iso,'density']; a=A(tfr/2.1)
    out=f"{nm:12}{tfr:6.2f}{a:8.3f}"; rat=[]
    for xi in (0.146,0.123):
        e=(a/a0)**(1/xi); out+=f"{e:10.3f}{D0/e:8.0f}"; rat.append(d/(D0/e))
    print(out+f"{d:8.0f}{rat[0]:10.2f}{rat[1]:10.2f}")
k=w[w.iso3=='KOR'].set_index('year')
print(f"\nKorea: 1960 pop {k.loc[1960,'pop_k']/1e3:.1f}M TFR {k.loc[1960,'tfr']:.2f}; 2023 pop {k.loc[2023,'pop_k']/1e3:.1f}M TFR {k.loc[2023,'tfr']:.2f}")
# macro panel
p=w[(w.pop_k>=1000)&(w.tfr<2*7.0)].copy(); p['n']=p.tfr/2.1; p=p[p.n<1/phi-0.5]
p['la']=np.log(A(p.n)); p['le']=-np.log(p.density)
def fe(df,fes,lab):
    X=pd.DataFrame({'le':df['le'].values},index=df.index)
    for f in fes: X=X.join(pd.get_dummies(df[f],prefix=f,drop_first=True).astype(float))
    m=sm.OLS(df['la'],sm.add_constant(X)).fit(cov_type='cluster',cov_kwds={'groups':df.iso3})
    print(f"{lab:46s} slope={m.params['le']:+.3f} ({m.bse['le']:.3f})  N={len(df)} countries={df.iso3.nunique()}")
print("\nMacro panel: ln alpha-hat on ln(1/density)")
fe(p[p.year==2023],[],"cross-section 2023")
fe(p,['year'],"pooled, year FE (between-country)")
fe(p,['iso3','year'],"country + year FE (within-country)")
fe(p[p.year>=1990],['iso3','year'],"country + year FE, 1990-2023")
