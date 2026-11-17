# Population growth in the long run

Private working repository containing the manuscript source, empirical analysis code and replication results prepared in this project.

## Contents

- `revised/`: Manuscript `Population_growth_22.tex` (author-supplied version 22, revised for consistency), matching `resources.bib`, and the three existing PDF figures referenced by that manuscript. The older version-3.8 manuscript remains available in Git history. The PNG figure is retained as an earlier replication output.
- `replication/`: WVS and CEX preparation/audit programs, original analysis scripts and reports.
- `replication/results/`: aggregate estimates, sample counts, bootstrap diagnostics and simulated coefficient draws. No respondent-level extracts are included.
- `data/`: acquisition instructions and input checksums. Raw and prepared datasets remain local.

## Status

WVS version-5 headline gradients and structural-cell regressions have been rerun. CEX sample and point estimates reproduce the supplied code; full bootstrap outputs are included for review. GSS and WPP inputs have been downloaded locally, but their estimates have not yet been replicated here. Not every manuscript claim or robustness exercise has been independently verified.

The scripts currently contain machine-specific paths. Configure those paths before running them; see replication/README.md. This is a research work-in-progress, not yet a portable one-command replication package.

PDF compilation is performed only on explicit request. Existing PDF figure assets are included; no compiled manuscript PDF or automatic compilation workflow is included. No open-source or data redistribution license is granted by this upload.

