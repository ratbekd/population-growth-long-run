# Empirical replication audit

This directory reproduces the CEX sample and point estimates from the supplied code, then audits inference. It does not yet constitute replication of the WVS, GSS or WPP analyses.

## Run

1. Set `SOURCE` and `CACHE` in `prepare_cex.py` to the archive directory and a local working directory. The input directory must contain `intrvw00.zip` through `intrvw25.zip` and `compute/analysis_code.zip`.
2. Set `LOCAL` in `audit_cex.py` to the local dependency/cache directory. The current paths are configured for this machine.
3. Run `python prepare_cex.py`, then `python audit_cex.py --bootstrap 999 --legacy-bootstrap 200`.

Package versions and platform are recorded in `results/environment.json`. The installed extra packages are in the configured `cex_deps` directory. The preparation step uses pandas and preserves the supplied construction code in `original/`. Only run the preparation step with the reviewed, supplied archive: it executes its `cex_build.py`.

## What is repaired

The audit retains the original sample, outcome definitions, controls and survey weights. It uses converged PPML with a line search and checks the normalized score; point estimates are independently compared with statsmodels GLM/WLS. Bootstrap sampling resamples the original location clusters, with seed 3. Repeated clusters enter as multiplicities in observation weights.

The new bootstrap does not reject draws because their coefficients are distant from the point estimate. Every attempted draw is logged; failed fits are reported by outcome. A failed fit receives a second attempt with a fresh start and a higher iteration limit. Rank deficiencies in nuisance controls use a least-squares solution; the habitat coefficient must remain estimable in the weighted design. This matters when cluster sampling removes all observations in a categorical control level. Residual numerical failures must be inspected before publishing inference; successful draws alone are not automatically an adequate basis for inference.

The first 200 bootstrap resamples also run the original ten-iteration solver and coefficient filter for comparison. The new standard errors use the sample standard deviation, not the IQR divided by 1.349 used in the supplied panel-C variant. Ratio intervals are percentile intervals of paired numerator/children-coefficient draws, with nonpositive denominators counted.

The legacy comparison uses algebraically equivalent sparse weighted cross-products for speed. Floating-point summation may differ from the dense original, so it is a reconstruction of the algorithm rather than a claim of bit-for-bit reproduction.

## Input and sample documentation

`results/input_manifest.json` records archive members, selected files, their decompressed SHA-256 hashes and row counts. First-quarter overlap follows the original preference for `fmliYY1x.csv`. The raw cache has 662,987 interviews; the constructed cache has 236,529 consumer units. Records include interviews in early 2026; this should not be described as a complete 2026 survey year.

`results/sample_flow.csv` documents sequential restrictions. Household characteristics are averages across interviews, with categorical values taken as the first nonmissing observation. Reference age is restricted to 25–44; the rounded mean adult count must be one or two. The main wage sample additionally requires hourly wages from 3 to 200, mean reported earners at most two, and positive total and noneducation expenditure. These filters produce 64,853 units and 139 clusters. The broader panel-C base has 74,902 units.

Education divided by hourly wage is capped at its sample 99th percentile. It is not education spending alone that is winsorized. Missing expenditure components and hours inputs are set to zero in the original construction. These choices are preserved for reproduction and need separate sensitivity analysis. A match to the original estimates does not validate the wage measure, the clustering design or a causal interpretation.

## Remaining data work

The supplied WVS/EVS dictionary and merge syntax are metadata, not respondent data. The merge syntax names WVS version 4.0 and EVS version 3.0. The WVS download page currently offers time-series version 5.0. A replication using version 5.0 must disclose the change and compare coverage; it cannot be labelled exact version-4 replication. WVS version 5.0 was downloaded on 12 November 2026 with user authorization; its ZIP integrity check passed. The archive and checksum manifest are in data/. Do not redistribute restricted raw data in a replication package.

The dictionary places ideal number of children (`D017`) in WVS waves 1–4. The original plotting code restricts both ideal and realized fertility to waves 2–4; the text should state both restrictions. EVS should not be appended unless the analysis explicitly calls for it. GSS and WPP input datasets and the original WVS preprocessing cache remain unverified.


## Our World in Data checks

`owid_checks.py` reproduces every result in the manuscript that uses the Our World in Data CO2 file (`owid-co2-data.csv`, which carries Global Carbon Budget emissions and Maddison Project GDP): the compound emission-intensity rates and required clean effort of Table 5, the simulated country paths of Table D.1, the Section 5.2 simulation checks, the asymptotic decline rates for Korea and Japan, Table E.1, and the sample of Table 3, panel A. Run `python owid_checks.py path/to/owid-co2-data.csv`; it needs numpy and pandas only. The script uses the unrounded taste share (0.539) and pollution elasticity (0.068/0.289); with them it also reproduces the earlier version of Table 5 exactly and Table D.1 to the reported precision.

`fig_transitions.py` draws Figure 3 (`revised/fig_transitions.pdf`) from the same file and calibration: `python fig_transitions.py path/to/owid-co2-data.csv ../revised/fig_transitions.pdf`. It needs matplotlib in addition.
