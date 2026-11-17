import pandas as pd, numpy as np, statsmodels.api as sm
gam,eta,phi,th=0.26,0.07,0.138,0.83; Gam=gam*(eta+th)
d=pd.read_pickle("ivs_clean.pkl")
for c in ['X011','X049']: d[c]=d[c].where(d[c]>=0)
d['age']=d.X003.where(d.X003.between(15,99)); d['wt']=d.S017.where(d.S017>0,1.0)
d=d.dropna(subset=['X011','X049','age']); d['cw']=d.S003.astype(str)+"_"+d.S002.astype(str)
print("X011 top-coding: share at the max value =", f"{(d.X011>=5).mean():.3%}")

def slope(top=1_000_000, minN=30, agelo=45, agehi=70, sub=None):
    mid={1:1000,2:3162,3:7071,4:14142,5:31623,6:70711,7:223607,8:top}
    x=d[d.age.between(agelo,agehi)].copy()
    if sub is not None: x=x[sub(x)]
    x['pop']=x.X049.map(mid)
    c=(x.groupby(['cw','X049']).apply(lambda g: pd.Series(
        {'n':np.average(g.X011,weights=g.wt)/2,'N':len(g),'pop':g['pop'].iloc[0]}),
        include_groups=False).reset_index())
    c=c[c.N>=minN]
    c['lnalpha']=np.log((phi*c.n+Gam)/(1-phi*c.n)); c['lneps']=-np.log(c['pop'])
    X=pd.get_dummies(c.cw,prefix='w',drop_first=True).astype(float); X['lneps']=c['lneps']
    m=sm.WLS(c['lnalpha'],sm.add_constant(X),weights=c.N).fit(cov_type='cluster',cov_kwds={'groups':c.cw})
    return m.params['lneps'], m.bse['lneps'], len(c)

print(f"\n{'specification':38s}{'b':>9}{'se':>8}{'cells':>7}{'xi @beta=.85':>14}")
for lab,kw in [("baseline (top=1.0M, N>=30, 45-70)",{}),
               ("top bin midpoint 750k",dict(top=750_000)),
               ("top bin midpoint 2.0M",dict(top=2_000_000)),
               ("min cell 50",dict(minN=50)),
               ("min cell 15",dict(minN=15)),
               ("age 40-75",dict(agelo=40,agehi=75)),
               ("age 50-70",dict(agelo=50,agehi=70)),
               ("women only",dict(sub=lambda x:x.X001==2)),
               ("waves 5-7 only",dict(sub=lambda x:x.S002>=5)),
               ("waves 2-4 only",dict(sub=lambda x:x.S002<=4))]:
    b,se,nc=slope(**kw); print(f"{lab:38s}{b:>9.4f}{se:>8.4f}{nc:>7}{b/0.15:>14.3f}")

print("\nbeta that reproduces the calibrated xi=0.146:")
b,_,_=slope(); print(f"   1 - b/0.146 = {1-b/0.146:.3f}")
print("\nmonotonicity in alpha space (bin dummies, bin1 omitted, country-wave FE):")
mid={1:1000,2:3162,3:7071,4:14142,5:31623,6:70711,7:223607,8:1_000_000}
x=d[d.age.between(45,70)].copy(); x['pop']=x.X049.map(mid)
c=(x.groupby(['cw','X049']).apply(lambda g: pd.Series(
    {'n':np.average(g.X011,weights=g.wt)/2,'N':len(g)}),include_groups=False).reset_index())
c=c[c.N>=30]; c['lnalpha']=np.log((phi*c.n+Gam)/(1-phi*c.n))
X=pd.get_dummies(c.X049.astype(int),prefix='b',drop_first=True).astype(float)
X=X.join(pd.get_dummies(c.cw,prefix='w',drop_first=True).astype(float))
m=sm.WLS(c['lnalpha'],sm.add_constant(X),weights=c.N).fit(cov_type='cluster',cov_kwds={'groups':c.cw})
for k in range(2,9): print(f"   bin{k}: {m.params[f'b_{k}']:+.4f} ({m.bse[f'b_{k}']:.4f})")
