import pandas as pd, numpy as np
src=open('cex_est2.py').read().split("w=c[c.wage.between(3,200)")[0]
src=src.replace("bs=np.array(bs); bs=bs[(np.abs(bs-est)<0.5).all(1)]; se=bs.std(0)","bs=np.array(bs); bs=bs[(np.abs(bs-est)<0.5).all(1)]; se=(np.percentile(bs,75,axis=0)-np.percentile(bs,25,axis=0))/1.349; print('   coef 95% pct CI:',np.round(np.percentile(bs,[2.5,97.5],axis=0).T,4).tolist())")
exec(src)
c['hrs_']=c.hrs
run2(c,[('fertility','PERSLT18','ppml'),('annual hours incl. zeros (PPML)','hrs_','ppml')],"LABOUR SUPPLY: all 1-2 adult households (robust se)")
print("S at beta=.85 from headline pi_n=0.0370(0.0053):",round(0.0370/0.15,3),[round((0.0370-1.96*0.0053)/0.15,3),round((0.0370+1.96*0.0053)/0.15,3)])
