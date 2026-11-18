"""Refresh the WVS figure. By default PNG only; with --pdf also write revised/fig_gradient.pdf from the same
original fig.py (same style), and print the largest-bin coefficients with their clustered standard errors."""
import os,sys,zipfile
from pathlib import Path
os.environ['OPENBLAS_NUM_THREADS']='1'; os.environ['OMP_NUM_THREADS']='1'
LOCAL=Path(r'C:/Users/ratbekd/.codex/visualizations/2026/11/11/01a1c8fa-6b21-7f31-a0ee-8ee897cb4798')
sys.path.append(str(LOCAL/'cex_deps'))   # fallback only: its wheels are for Python 3.12
HERE=Path(__file__).resolve().parent
with zipfile.ZipFile(r'C:/Users/ratbekd/OneDrive - Monash University/Reserach/Population growth in the LR/compute/analysis_code.zip') as z:
 source=z.read('fig.py').decode('utf-8')
if '--pdf' in sys.argv:
 source=source.replace('"/home/claude/build/fig_gradient.pdf"',repr(str(HERE.parent/'revised'/'fig_gradient.pdf')))
else:
 source='\n'.join(line for line in source.splitlines() if 'fig.savefig' not in line or '.pdf' not in line)
source+='\nprint(f"largest bin: D017 {A[0][-1]:.3f} ({A[1][-1]:.3f})   X011 {B[0][-1]:.3f} ({B[1][-1]:.3f})")\n'
source=source.replace('"ivs_clean.pkl"',repr(str(HERE/'data'/'wvs_v5_selected.pkl')))
source=source.replace('"/home/claude/work/fig_gradient.png"',repr(str(HERE.parent/'revised'/'fig_gradient_v5.png')))
exec(compile(source,'fig.py','exec'),{'__name__':'__main__'})
