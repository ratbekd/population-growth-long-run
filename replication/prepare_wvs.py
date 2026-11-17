"""Read official WVS v5 data without guessing undocumented recodes."""
from pathlib import Path
import zipfile, json, pandas as pd
HERE=Path(__file__).resolve().parent
OUT=HERE/'results'/'wvs'; OUT.mkdir(exist_ok=True)
cols=['S001','S002VS','S003','S017','D017','X011','X049','X025','X047_WVS','F028','X007','X003','X001']
archive=HERE/'data'/'WVS_Time_Series_1981-2022_Stata_v5_0.zip'
with zipfile.ZipFile(archive) as z:
 with z.open('WVS_Time_Series_1981-2022_stata_v5_0.dta') as f:
  reader=pd.read_stata(f,iterator=True,convert_categoricals=False)
  labels=reader.variable_labels()
  missing=[c for c in cols if c not in labels]
  if missing: raise ValueError(f'Missing required columns: {missing}')
  d=reader.read(columns=cols)
  values=reader.value_labels()
  relevant={k:{str(code):label for code,label in v.items()} for k,v in values.items() if any(c.lower() in k.lower() for c in cols)}
(OUT/'metadata.json').write_text(json.dumps({'version':'5.0','rows':len(d),'variable_labels':{c:labels[c] for c in cols},'value_labels':relevant},indent=2),encoding='utf-8')
d=d.rename(columns={'S002VS':'S002'})
d.to_pickle(HERE/'data'/'wvs_v5_selected.pkl')
rows=[]
for wave,g in d.groupby('S002'):
 rows.append({'wave':int(wave),'N':len(g),'countries':g.S003.nunique(),**{c:int((g[c]>=0).sum()) for c in ['D017','X011','X049']}})
pd.DataFrame(rows).to_csv(OUT/'wave_coverage.csv',index=False)
for c in ['S001','D017','X011','X049']:
 d[c].value_counts(dropna=False).sort_index().rename('N').to_csv(OUT/f'codes_{c}.csv')
print('ROWS',len(d), 'COUNTRIES',d.S003.nunique(),flush=True)
print(pd.DataFrame(rows).to_string(index=False),flush=True)
print({c:labels[c] for c in cols},flush=True)
print('X011',d.X011.value_counts().sort_index().to_dict(),flush=True)
