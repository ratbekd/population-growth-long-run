import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from fullsim import simulate
T=50
runs=[("Decoupling ($\\zeta=7$)",7,0.0,dict(color="#111111",ls="-",lw=1.6)),
      ("Decoupling, remediation tax $\\tau=0.1$",7,0.1,dict(color="#111111",ls="--",lw=1.4)),
      ("Weak decoupling ($\\zeta=4$)",4,0.0,dict(color="#888888",ls="-",lw=1.6)),
      ("No decoupling ($\\zeta=0$)",0,0.0,dict(color="#888888",ls=":",lw=1.8))]
plt.rcParams.update({"font.family":"serif","font.size":9,"axes.linewidth":0.7})
fig,ax=plt.subplots(2,2,figsize=(7.0,5.0),sharex=True); ax=ax.ravel()
for lab,z,tau,st in runs:
    o=simulate(z,tau,T=T); t=o[:,0]; L,n,eps,y=o[:,1],o[:,2],o[:,3],o[:,6]
    gy=np.r_[np.nan,(y[1:]/y[:-1])**(1/30)-1]*100
    ax[0].plot(t,L,label=lab,**st); ax[1].plot(t,n,**st); ax[2].plot(t,eps,**st); ax[3].plot(t[2:],gy[2:],**st)
ax[0].set_yscale('log'); ax[1].axhline(1,color="#cccccc",lw=0.8,zorder=0)
for a,tt in zip(ax,["(a) Population $L_t$ (log scale)","(b) Fertility $n_t$ (replacement $=1$)",
                    "(c) Habitat per capita $\\epsilon_t$","(d) Income growth, % per year"]):
    a.set_title(tt,loc='left',fontsize=9)
    for sp in ('top','right'): a.spines[sp].set_visible(False)
    a.grid(axis='y',color="#e3e3e3",lw=0.5); a.set_axisbelow(True)
for a in ax[2:]: a.set_xlabel("generations (30 years)",fontsize=8)
ax[3].set_ylim(0,2.4)
fig.legend(*ax[0].get_legend_handles_labels(),loc='lower center',ncol=2,frameon=False,fontsize=8,bbox_to_anchor=(0.5,-0.01))
fig.tight_layout(rect=[0,0.07,1,1])
fig.savefig("/home/claude/fig_decoupling.pdf",bbox_inches='tight'); fig.savefig("fig_decoupling.png",dpi=160,bbox_inches='tight')
