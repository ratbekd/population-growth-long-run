"""Run reviewed original WVS analyses on the documented v5 preprocessing."""
import os,sys,zipfile,io,contextlib,json
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1'; os.environ['OMP_NUM_THREADS']='1'; os.environ['MKL_NUM_THREADS']='1'
LOCAL=Path(r'C:/Users/ratbekd/.codex/visualizations/2026/11/11/01a1c8fa-6b21-7f31-a0ee-8ee897cb4798')
sys.path.insert(0,str(LOCAL/'cex_deps'))
import pandas as pd
HERE=Path(__file__).resolve().parent; OUT=HERE/'results'/'wvs'
raw=HERE/'data'/'wvs_v5_selected.pkl'
d=pd.read_pickle(raw)
assert set(d.S001.unique())=={2},'Only WVS observations may use this wave alias'
assert set(d.S002.unique())<=set(range(1,8))
assert set(d.X049[d.X049>=0].unique())<=set(range(1,9))
coverage={}
for y in ['D017','X011']:
 g=d[(d[y]>=0)&(d.X049>=0)]
 coverage[y]={'N':len(g),'countries':int(g.S003.nunique()),'country_waves':len(g[['S003','S002']].drop_duplicates())}
(OUT/'analysis_coverage.json').write_text(json.dumps(coverage,indent=2))
with zipfile.ZipFile(r'C:/Users/ratbekd/OneDrive - Monash University/Reserach/Population growth in the LR/compute/analysis_code.zip') as z:
 for name in ['grad.py','xi.py','wvs2.py']:
  source=z.read(name).decode('utf-8')
  (HERE/'original'/name).write_text(source,encoding='utf-8')
  source=source.replace('"ivs_clean.pkl"',repr(str(raw))).replace('"cells.pkl"',repr(str(OUT/'cells.pkl'))).replace('"/home/claude/work/cells.pkl"',repr(str(OUT/'cells.pkl'))).replace('"grad_res.npy"',repr(str(OUT/'grad_res.npy')))
  log=io.StringIO()
  with contextlib.redirect_stdout(log): exec(compile(source,name,'exec'),{'__name__':'__main__'})
  (OUT/(name+'.log')).write_text(log.getvalue(),encoding='utf-8')
  print(name,log.getvalue(),flush=True)
