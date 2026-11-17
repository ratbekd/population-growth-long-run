import pandas as pd, numpy as np, statsmodels.api as sm
d=pd.read_pickle("ivs_clean.pkl")
for c in ['D017','X011','X049','X025','X047_WVS','F028','X007']:
    d[c]=d[c].where(d[c]>=0)
d['age']=d.X003.where(d.X003.between(15,99)); d['sex']=d.X001.where(d.X001.isin([1,2]))
d['cw']=d.S003.astype(str)+"_"+d.S002.astype(str)
d['wt']=d.S017.where(d.S017>0,1.0)

def design(g, ctrl):
    X=pd.get_dummies(g.X049.astype(int),prefix='bin',drop_first=True).astype(float)
    X=X.join(pd.get_dummies(g.cw,prefix='cw',drop_first=True).astype(float))
    X['age']=g.age; X['age2']=g.age**2/100; X['female']=(g.sex==2).astype(float)
    if ctrl>=1:
        X=X.join(pd.get_dummies(g.X007.astype('Int64'),prefix='mar',drop_first=True).astype(float))
        X['relig']=g.F028
    if ctrl>=2:
        X['educ']=g.X025; X['inc']=g.X047_WVS
    return sm.add_constant(X)

def run(y,ctrl,lab,sample):
    g=sample.dropna(subset=[y,'X049','age','sex']+(['X007','F028'] if ctrl>=1 else [])
                    +(['X025','X047_WVS'] if ctrl>=2 else [])).copy()
    X=design(g,ctrl)
    m=sm.WLS(g[y],X,weights=g.wt).fit(cov_type='cluster',cov_kwds={'groups':g.cw})
    b=[m.params.get(f'bin_{k}',np.nan) for k in range(2,9)]
    se=[m.bse.get(f'bin_{k}',np.nan) for k in range(2,9)]
    print(f"{lab:34s} N={len(g):>6}  "+" ".join(f"{x:+.3f}" for x in b))
    print(f"{'':34s} {'':>8}  "+" ".join(f"({x:.3f})" for x in se))
    return b,se,len(g)

s=d[d.S002.isin([2,3,4])]
print("Settlement-size coefficients (bin 1 = <2,000 omitted), country x wave FE, clustered\n")
print(f"{'':34s} {'':>8}  "+" ".join(f"{'bin'+str(k):>7}" for k in range(2,9)))
res={}
for ctrl,lab in [(0,'D017: FE + age,sex'),(1,'D017: + marital, religiosity'),
                 (2,'D017: + education, income  [bad]')]:
    res[lab]=run('D017',ctrl,lab,s)
print()
for ctrl,lab in [(0,'X011: FE + age,sex'),(1,'X011: + marital, religiosity'),
                 (2,'X011: + education, income  [bad]')]:
    run('X011',ctrl,lab,s[s.age.between(45,70)])
np.save("grad_res.npy",np.array(res['D017: + marital, religiosity'][:2]))
