import pandas as pd, numpy as np, statsmodels.api as sm
gam,eta,phi,th=0.26,0.07,0.138,0.83; Gam=gam*(eta+th); a0,XI=0.39,0.146
d=pd.read_pickle("ivs_clean.pkl")
for c in ['X011','X049']: d[c]=d[c].where(d[c]>=0)
d['age']=d.X003.where(d.X003.between(15,99))
d['wt']=d.S017.where(d.S017>0,1.0)
d=d.dropna(subset=['X011','X049','age'])
d=d[d.age.between(45,70)]                      # completed fertility
d['cw']=d.S003.astype(str)+"_"+d.S002.astype(str)

mid={1:1000,2:3162,3:7071,4:14142,5:31623,6:70711,7:223607,8:1000000}
d['pop']=d.X049.map(mid)

cell=(d.groupby(['cw','X049'])
        .apply(lambda g: pd.Series({'n_bar':np.average(g.X011,weights=g.wt)/2.0,
                                    'N':len(g),'pop':g['pop'].iloc[0]}), include_groups=False)
        .reset_index())
cell=cell[cell.N>=30]
cell['alpha']=(phi*cell.n_bar+Gam)/(1-phi*cell.n_bar)
cell['lnalpha']=np.log(cell.alpha); cell['lneps_naive']=-np.log(cell['pop'])
print(f"cells={len(cell)}  country-waves={cell.cw.nunique()}  "
      f"n_bar: mean={cell.n_bar.mean():.3f} min={cell.n_bar.min():.3f} max={cell.n_bar.max():.3f}")
print(f"implied alpha: mean={cell.alpha.mean():.4f}  [calibration: a0={a0}, a*={(phi+Gam)/(1-phi):.4f}]\n")

X=pd.get_dummies(cell.cw,prefix='cw',drop_first=True).astype(float)
X['lneps']=cell.lneps_naive; X=sm.add_constant(X)
m=sm.WLS(cell.lnalpha,X,weights=cell.N).fit(cov_type='cluster',cov_kwds={'groups':cell.cw})
b,se=m.params['lneps'], m.bse['lneps']
print(f"slope of ln(alpha) on -ln(settlement population):  b = {b:.4f}  (se {se:.4f})")
print(f"   = xi * (d ln density / d ln population)\n")
print(f"{'assumed area ~ pop^beta':>26}{'d ln dens/d ln pop':>20}{'implied xi':>13}{'95% CI':>20}")
for beta,lab in [(0.0,'area fixed (naive)'),(0.85,'urban scaling 0.85'),(1.0,'area ~ pop')]:
    kD=1-beta
    if kD==0: print(f"{lab:>26}{kD:>20.2f}{'undefined':>13}{'':>20}"); continue
    xi=b/kD; lo,hi=(b-1.96*se)/kD,(b+1.96*se)/kD
    print(f"{lab:>26}{kD:>20.2f}{xi:>13.3f}   [{lo:.3f}, {hi:.3f}]")
print(f"\ncalibrated value in the paper: xi = {XI}")
cell.to_pickle("cells.pkl")
