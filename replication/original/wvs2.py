import pandas as pd, numpy as np, statsmodels.api as sm
gam,eta,phi,th=0.26,0.07,0.138,0.83; G=gam*(eta+th)
K=lambda a:a*(1+G)/((a-G)*(1+a))
c=pd.read_pickle("/home/claude/work/cells.pkl")
c['lnn']=np.log(c.n_bar); c['x']=-np.log(c['pop'])
X=pd.get_dummies(c.cw,prefix='w',drop_first=True).astype(float); X['x']=c['x']; X=sm.add_constant(X)
for y,lab in [('lnn','ln n (direct fertility elasticity)'),('lnalpha','ln alpha-hat (inversion)')]:
    m=sm.WLS(c[y],X,weights=c.N).fit(cov_type='cluster',cov_kwds={'groups':c.cw})
    b,se=m.params['x'],m.bse['x']
    print(f"{lab:38s} slope x(1-beta) = {b:.4f} ({se:.4f});  at beta=.85: {b/.15:.3f} [{(b-1.96*se)/.15:.3f},{(b+1.96*se)/.15:.3f}]")
    if y=='lnn': S,Slo,Shi=b/.15,(b-1.96*se)/.15,(b+1.96*se)/.15
abar=np.average(c.alpha,weights=c.N); print(f"\nweighted mean alpha-hat={abar:.4f}  K(mean)={K(abar):.3f}  K(alpha*)={K((phi+G)/(1-phi)):.3f}")
print(f"S-hat = {S:.3f} [{Slo:.3f},{Shi:.3f}]  -> f'={1-S:.3f}, half-life={np.log(.5)/np.log(1-S):.2f} gens "
      f"[{np.log(.5)/np.log(1-Shi):.2f},{np.log(.5)/np.log(1-Slo):.2f}]")
print("implied taste elasticity for assumed chi:  xi = (S - chi)/K(mean alpha)")
for chi in (0,0.05,0.1,0.15,S): print(f"   chi={chi:.3f}: xi={(S-chi)/K(abar):.3f}")
for b_ in (0.80,0.85,0.90): 
    m=sm.WLS(c['lnn'],X,weights=c.N).fit(); print(f"   beta={b_}: S={m.params['x']/(1-b_):.3f}")
