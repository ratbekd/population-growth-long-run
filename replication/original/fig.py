import pandas as pd, numpy as np, statsmodels.api as sm
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
d=pd.read_pickle("ivs_clean.pkl")
for c in ['D017','X011','X049','F028','X007']: d[c]=d[c].where(d[c]>=0)
d['age']=d.X003.where(d.X003.between(15,99)); d['sex']=d.X001.where(d.X001.isin([1,2]))
d['cw']=d.S003.astype(str)+"_"+d.S002.astype(str); d['wt']=d.S017.where(d.S017>0,1.0)

def coefs(y,sample):
    g=sample.dropna(subset=[y,'X049','age','sex','X007','F028']).copy()
    X=pd.get_dummies(g.X049.astype(int),prefix='bin',drop_first=True).astype(float)
    X=X.join(pd.get_dummies(g.cw,prefix='cw',drop_first=True).astype(float))
    X=X.join(pd.get_dummies(g.X007.astype('Int64'),prefix='mar',drop_first=True).astype(float))
    X['age']=g.age; X['age2']=g.age**2/100; X['female']=(g.sex==2).astype(float); X['rel']=g.F028
    m=sm.WLS(g[y],sm.add_constant(X),weights=g.wt).fit(cov_type='cluster',cov_kwds={'groups':g.cw})
    b=np.array([0]+[m.params[f'bin_{k}'] for k in range(2,9)])
    se=np.array([0]+[m.bse[f'bin_{k}'] for k in range(2,9)])
    return b,se,len(g),g.cw.nunique()

s=d[d.S002.isin([2,3,4])]
A=coefs('D017',s); B=coefs('X011',s[s.age.between(45,70)])
lab=["<2k","2–5k","5–10k","10–20k","20–50k","50–100k","0.1–0.5M",">0.5M"]
INK="#1a1a1a"; GRID="#d8d8d8"
plt.rcParams.update({"font.family":"serif","font.size":9,"axes.linewidth":0.7})
fig,axes=plt.subplots(1,2,figsize=(7.0,2.9),sharex=True)
for ax,(b,se,N,C),ttl in zip(axes,[A,B],
        ["(a)  Ideal number of children","(b)  Completed fertility, age 45–70"]):
    x=np.arange(8)
    ax.axhline(0,color=GRID,lw=0.8,zorder=1)
    ax.errorbar(x,b,yerr=1.96*se,fmt='o',ms=4.5,lw=0,elinewidth=1.3,capsize=2.6,
                color=INK,ecolor=INK,zorder=3)
    ax.plot(x,b,lw=1.0,color=INK,alpha=0.45,zorder=2)
    ax.set_title(ttl,fontsize=9,loc='left',pad=6)
    ax.set_xticks(x); ax.set_xticklabels(lab,rotation=45,ha='right',fontsize=7.5)
    ax.grid(axis='y',color=GRID,lw=0.5,alpha=0.8); ax.set_axisbelow(True)
    for sp in ('top','right'): ax.spines[sp].set_visible(False)
    ax.spines['left'].set_color('#999'); ax.spines['bottom'].set_color('#999')
    ax.annotate(f"N = {N:,}   country×wave = {C}",xy=(0.03,0.06),xycoords='axes fraction',
                fontsize=7,color='#555')
axes[0].set_ylabel("coefficient relative to\nsmallest settlement",fontsize=8)
fig.supxlabel("size of town (X049)",fontsize=8,y=0.02)
fig.tight_layout(rect=[0,0.03,1,1])
fig.savefig("/home/claude/build/fig_gradient.pdf",bbox_inches='tight')
fig.savefig("/home/claude/work/fig_gradient.png",dpi=170,bbox_inches='tight')
print("A:",np.round(A[0],3),"N",A[2]); print("B:",np.round(B[0],3),"N",B[2])
