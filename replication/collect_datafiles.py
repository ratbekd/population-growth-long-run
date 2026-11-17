"""Consolidate local research data and fetch official GSS/WPP releases."""
from pathlib import Path
import shutil,hashlib,json,urllib.request,zipfile,concurrent.futures
ROOT=Path(__file__).resolve().parents[1]; DEST=ROOT/'Datafiles'
SOURCE=Path(r'C:/Users/ratbekd/OneDrive - Monash University/Reserach/Population growth in the LR')
DOWNLOADS=Path(r'C:/Users/ratbekd/Downloads')
LOCAL=Path(r'C:/Users/ratbekd/.codex/visualizations/2026/11/11/01a1c8fa-6b21-7f31-a0ee-8ee897cb4798')
DEST.mkdir(exist_ok=True)
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
 return h.hexdigest()
rows=[]
def copy(src,relative):
 target=DEST/relative; target.parent.mkdir(parents=True,exist_ok=True)
 if not src.exists(): print('MISSING',src,flush=True); return
 checksum=digest(src)
 if not target.exists() or digest(target)!=checksum: shutil.copy2(src,target)
 assert digest(target)==checksum,str(target)
 rows.append({'file':relative,'source':str(src),'bytes':target.stat().st_size,'sha256':checksum})
 print('COPIED',relative,flush=True)
def download(item):
 relative,url=item; target=DEST/relative; target.parent.mkdir(parents=True,exist_ok=True)
 if not target.exists():
  partial=target.with_suffix(target.suffix+'.partial')
  request=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'})
  with urllib.request.urlopen(request,timeout=90) as response, partial.open('wb') as f: shutil.copyfileobj(response,f,4*1024*1024)
  partial.replace(target)
 with zipfile.ZipFile(target) as z:
  assert z.testzip() is None
  members=z.namelist()
 row={'file':relative,'source':url,'bytes':target.stat().st_size,'sha256':digest(target),'archive_members':members,'retrieved':'2026-11-17'}
 print('DOWNLOADED',relative,row['bytes'],flush=True)
 return row
urls=[('raw/GSS/GSS_1972_2024_release3a_stata.zip','https://gss.norc.org/content/dam/gss/get-the-data/documents/stata/GSS_stata.zip'),('raw/WPP/WPP2024_GEN_F01_DEMOGRAPHIC_INDICATORS_FULL.xlsx','https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/EXCEL_FILES/1_General/WPP2024_GEN_F01_DEMOGRAPHIC_INDICATORS_FULL.xlsx')]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 futures=[pool.submit(download,item) for item in urls]
 for f in futures: rows.append(f.result())
for y in range(26): copy(SOURCE/f'intrvw{y:02d}.zip',f'raw/CEX/intrvw{y:02d}.zip')
copy(ROOT/'replication/data/WVS_Time_Series_1981-2022_Stata_v5_0.zip','raw/WVS/WVS_Time_Series_1981-2022_Stata_v5_0.zip')
copy(DOWNLOADS/'ZA7503_v3-0-0.dta.zip','raw/EVS/ZA7503_v3-0-0.dta.zip')
for name in ['ZA7503_Common_EVS_WVS_Dictionary_IVS.xlsx','ZA7503_EVS_WVS_Merge_Syntax_stata.zip']:
 copy(DOWNLOADS/name,'documentation/'+name)
for name in ['cex_raw.pkl','cex_cu.pkl']: copy(LOCAL/'cex_replication'/name,'prepared/CEX/'+name)
copy(ROOT/'replication/data/wvs_v5_selected.pkl','prepared/WVS/wvs_v5_selected.pkl')
copy(ROOT/'replication/results/wvs/cells.pkl','prepared/WVS/cells.pkl')
for p in (SOURCE/'compute/empty-planet-replication/data').rglob('*'):
 if p.is_file(): copy(p,'other_project_data/'+str(p.relative_to(SOURCE/'compute/empty-planet-replication/data')).replace('\\','/'))
copy(SOURCE/'compute/analysis_code.zip','code/analysis_code.zip')
for p in (ROOT/'replication').glob('*.py'): copy(p,'code/'+p.name)
for p in (ROOT/'replication/results').rglob('*'):
 if p.is_file() and p.suffix!='.pkl': copy(p,'results/'+str(p.relative_to(ROOT/'replication/results')).replace('\\','/'))
(DEST/'manifest.json').write_text(json.dumps(rows,indent=2),encoding='utf-8')
(DEST/'.gitignore').write_text('*\n!.gitignore\n!README.md\n!manifest.json\n')
print('COMPLETE',len(rows),'files',sum(r['bytes'] for r in rows),'bytes',flush=True)
