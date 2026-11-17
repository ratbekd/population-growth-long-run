from pathlib import Path
import hashlib, json, zipfile, io, os, runpy
import pandas as pd
SOURCE=Path(r'C:\Users\ratbekd\OneDrive - Monash University\Reserach\Population growth in the LR')
CACHE=Path(r'C:\Users\ratbekd\.codex\visualizations\2026\11\11\01a1c8fa-6b21-7f31-a0ee-8ee897cb4798\cex_replication')
HERE=Path(__file__).resolve().parent
RESULTS=HERE/'results'
for p in [CACHE,RESULTS,HERE/'original']: p.mkdir(parents=True,exist_ok=True)
with zipfile.ZipFile(SOURCE/'compute/analysis_code.zip') as z:
 for name in z.namelist():
  if name.startswith('cex_') and name.endswith('.py'): (HERE/'original'/name).write_bytes(z.read(name))
keep=['NEWID','QINTRVMO','QINTRVYR','FINLWT21','EDUCAPQ','EDUCACQ','TOTEXPPQ','TOTEXPCQ','FSALARYX','FSALARYM','FINCBTAX','FINCBTXM','FAM_SIZE','PERSLT18','AGE_REF','AGE2','SEX_REF','MARITAL1','EDUC_REF','REF_RACE','REGION','BLS_URBN','POPSIZE','SMSASTAT','STATE','PSU','INC_HRS1','INC_HRS2','INCWEEK1','INCWEEK2','FAM_TYPE','NO_EARNR','CUTENURE','SHELTPQ','SHELTCQ','ROOMSQ']
entries=[]
for y in range(26):
 archive=SOURCE/f'intrvw{y:02d}.zip'
 with zipfile.ZipFile(archive) as z:
  for e in z.infolist():
   name=Path(e.filename).name.lower()
   if name.startswith('fmli') and name.endswith('.csv'): entries.append(dict(archive=str(archive),member=e.filename,file=name,size=e.file_size))
names={e['file'] for e in entries}; frames=[]; manifest=[]; missing=[]
for e in sorted(entries,key=lambda e:e['file']):
 name=e['file']; b=name[4:-4]
 skipped=b[2]=='1' and b[3:]=='' and f'fmli{b[:2]}1x.csv' in names
 e['selected']=not skipped
 if skipped:
  e['reason']='Prefer x file, matching original loader'; manifest.append(e); continue
 with zipfile.ZipFile(e['archive']) as z: raw=z.read(e['member'])
 e['sha256_uncompressed']=hashlib.sha256(raw).hexdigest()
 d=pd.read_csv(io.BytesIO(raw),encoding='latin-1',dtype=str,usecols=lambda c:c.strip().upper() in keep,low_memory=False)
 d.columns=[c.strip().upper() for c in d.columns]; e['rows']=len(d); e['missing_columns']=sorted(set(keep)-set(d.columns))
 missing.append({'file':name,**{c:int(d[c].isna().sum()) if c in d else None for c in keep}})
 d['file']=b; frames.append(d); manifest.append(e)
 if b[2]=='4': print(name,len(d),flush=True)
d=pd.concat(frames,ignore_index=True)
for c in d.columns:
 if c not in ('file','STATE','PSU'): d[c]=pd.to_numeric(d[c],errors='coerce')
d.to_pickle(CACHE/'cex_raw.pkl')
pd.DataFrame(missing).to_csv(RESULTS/'raw_missing_counts.csv',index=False)
d.groupby(['QINTRVYR','QINTRVMO']).size().rename('interviews').to_csv(RESULTS/'interviews_by_month.csv')
(RESULTS/'input_manifest.json').write_text(json.dumps(manifest,indent=2))
print('RAW',d.shape,'YEARS',d.QINTRVYR.min(),d.QINTRVYR.max(),'DUPLICATE NEWID',d.duplicated(['NEWID','QINTRVYR','QINTRVMO']).sum(),flush=True)
os.chdir(CACHE)
runpy.run_path(str(HERE/'original/cex_build.py'),run_name='__main__')
