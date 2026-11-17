# Obtaining the data

Download data directly from the official providers and comply with their use conditions. Raw files and respondent-level prepared extracts are excluded from this repository, including private uploads. WVS expressly prohibits redistribution.

| Dataset | Version used or obtained | Official source |
|---|---|---|
| CEX Interview Survey | intrvw00.zip through intrvw25.zip | https://www.bls.gov/cex/pumd_data.htm |
| WVS | 1981–2022 time series, Stata v5.0; downloaded 12 November 2026 | https://www.worldvaluessurvey.org/WVSDocumentationWVL.jsp |
| EVS | ZA7503 v3.0.0, supplied locally; not yet merged | https://doi.org/10.4232/1.14021 |
| GSS | 1972–2024 cumulative, Stata release 3a, July 2026; downloaded 17 November 2026 | https://gss.norc.org/get-the-data/stata.html |
| UN WPP | WPP2024_GEN_F01_DEMOGRAPHIC_INDICATORS_FULL.xlsx; downloaded 17 November 2026 | https://population.un.org/wpp/ |

`source_checksums.json` records exact archive/file sizes and SHA-256 hashes. WVS download requires registration and acceptance of the provider's terms; this repository does not bypass that process. The supplied EVS/WVS merge syntax references WVS v4.0, whereas the replicated WVS analysis uses v5.0. The original GSS code reads an R-format cumulative file; the newly obtained official Stata release requires a documented adapter before analysis.

Place WVS in `replication/data/` using the basename expected by `prepare_wvs.py`. For CEX, set SOURCE and CACHE in `prepare_cex.py`; set LOCAL in `audit_cex.py`. Run WVS preparation before `run_wvs.py`. These directories and respondent-level caches must remain untracked. See the replication reports for exact sample restrictions and limitations.
