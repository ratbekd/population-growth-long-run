import pyreadr, gc, time
t=time.time(); r=pyreadr.read_r("/home/claude/data/gss/gss_all.rda"); k=list(r.keys()); print(k, flush=True)
df=r[k[0]]; del r; gc.collect(); print(df.shape, round(time.time()-t), flush=True)
want=['year','id','size','xnorcsiz','srcbelt','childs','chldidel','hrs1','hrs2','wrkstat','age','sex','educ','degree',
      'marital','attend','realinc','realrinc','conrinc','race','region','wtssall','wtssps','wtssnrps','wtsscomp','vstrat','vpsu',
      'sphrs1','spwrksta','hompop','babies','preteen','teens','agekdbrn','born','res16','reg16']
have=[c for c in want if c in df.columns]; print("missing:", set(want)-set(have), flush=True)
df[have].to_pickle("/home/claude/work3/gss_sub.pkl"); print("saved", flush=True)
