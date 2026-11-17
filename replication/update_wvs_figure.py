"""Refresh the WVS figure as PNG only; never compile or generate a PDF."""
import os,sys,zipfile
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1'; os.environ['OMP_NUM_THREADS']='1'
LOCAL=Path(r'C:/Users/ratbekd/.codex/visualizations/2026/11/11/01a1c8fa-6b21-7f31-a0ee-8ee897cb4798')
sys.path.insert(0,str(LOCAL/'cex_deps'))
HERE=Path(__file__).resolve().parent
with zipfile.ZipFile(r'C:/Users/ratbekd/OneDrive - Monash University/Reserach/Population growth in the LR/compute/analysis_code.zip') as z:
 source=z.read('fig.py').decode('utf-8')
source='\n'.join(line for line in source.splitlines() if 'fig.savefig' not in line or '.pdf' not in line)
source=source.replace('"ivs_clean.pkl"',repr(str(HERE/'data'/'wvs_v5_selected.pkl')))
source=source.replace('"/home/claude/work/fig_gradient.png"',repr(str(HERE.parent/'revised'/'fig_gradient_v5.png')))
exec(compile(source,'fig.py','exec'),{'__name__':'__main__'})
