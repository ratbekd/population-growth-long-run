# Population growth in the long run

Private working repository containing the manuscript source, empirical analysis code and replication results prepared in this project.

## Contents

- `revised/`: LaTeX source, bibliography and WVS figure. This is the working revision developed in this chat, not a claim to be the latest manuscript elsewhere on the author's computer.
- `replication/`: WVS and CEX preparation/audit programs, original analysis scripts and reports.
- `replication/results/`: aggregate estimates, sample counts, bootstrap diagnostics and simulated coefficient draws. No respondent-level extracts are included.
- `data/`: acquisition instructions and input checksums. Raw and prepared datasets remain local.

## Status

WVS version-5 headline gradients and structural-cell regressions have been rerun. CEX sample and point estimates reproduce the supplied code; full bootstrap outputs are included for review. GSS and WPP inputs have been downloaded locally, but their estimates have not yet been replicated here. Not every manuscript claim or robustness exercise has been independently verified.

The scripts currently contain machine-specific paths. Configure those paths before running them; see replication/README.md. This is a research work-in-progress, not yet a portable one-command replication package.

PDF compilation is performed only on explicit request. No PDFs or automatic compilation workflow are included. No open-source or data redistribution license is granted by this upload.
