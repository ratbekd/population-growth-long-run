# WVS version-5 replication

The downloaded respondent data are usable. The main settlement-gradient and structural-cell regressions have now been rerun using the supplied analysis scripts and an explicit preprocessing step.

## Data and construction

The official WVS 1981–2022 Stata v5.0 archive contains 443,488 respondents from 108 country codes. All selected observations belong to WVS (S001 = 2). The file calls its chronology variable S002VS; `prepare_wvs.py` aliases it to S002 for compatibility with the original scripts, after retaining the original metadata. Its labels identify waves 1–7. No EVS observations are appended.

The selected variables are saved locally in `data/wvs_v5_selected.pkl`. Raw negative missing codes are retained there; the supplied estimation scripts convert them to missing, apply their age restrictions and use S017 weights (replacing nonpositive or missing weights with one, as in the original code). X011 code 5 is explicitly labelled “5 children or more”; it is not an uncensored count. X049 has the eight expected settlement categories.

Version 5.0 differs from the file described in the paper and from the version-4 input named in the supplied merge syntax. The following is a replication on v5.0, not proof of exact reproduction of the earlier dataset.

| Coverage | Paper | WVS v5.0 |
|---|---:|---:|
| All respondents | 442,473 | 443,488 |
| Ideal fertility and settlement size observed | 98,512 | 99,126 |
| Countries in that overlap | 51 | 52 |
| Country–waves in that overlap | 77 | 78 |
| Realized fertility and settlement size observed | 316,751 | 317,428 |
| Countries in that overlap | 100 | 100 |
| Country–waves in that overlap | 231 | 232 |

## Main estimates

Both gradient regressions below use waves 2–4, country–wave fixed effects, age, age squared, sex, marital status and religiosity, with standard errors clustered by country–wave. The realized-fertility regression additionally restricts age to 45–70. Coefficients compare the largest settlement category with settlements under 2,000 people.

| Estimate | Paper | WVS v5.0 |
|---|---:|---:|
| Ideal fertility, largest-settlement coefficient | −0.330 (SE 0.051) | −0.330 (SE 0.050) |
| Realized fertility, largest-settlement coefficient | −0.504 (SE 0.051) | −0.499 (SE 0.051) |
| Ideal fertility, adding education and income | −0.275 | −0.275 |
| Realized fertility, adding education and income | −0.393 | −0.390 |

The headline ideal-fertility regression uses 95,985 respondents, and the headline realized-fertility regression uses 31,842. The declining gradients remain monotone in the headline specifications. The caption should explicitly state that **both panels use waves 2–4**.

The structural-cell analysis uses all available waves with settlement data, ages 45–70, and a minimum cell size of 30. It yields 1,034 cells in 226 country–waves, compared with 1,031 and 225 in the paper. Its principal estimates remain unchanged at the paper's rounding:

- Composite slope of log implied taste on minus log settlement population: **0.0185**, SE **0.0012**.
- Implied taste elasticity under the maintained area-scaling assumption beta = 0.85: **0.123**, 95% interval **[0.108, 0.139]**.
- Direct fertility feedback estimate under the same scaling assumption: **0.209**, 95% interval **[0.182, 0.236]**.

The unweighted mean implied taste across cells is 0.5027; the respondent-count-weighted mean is 0.4928. These should not be conflated.

## Scope and reproducibility

Run `prepare_wvs.py`, then `run_wvs.py`, using the configured Python dependencies. The wrapper preserves original scripts in `original/` and changes only input/output paths. Results are in `results/wvs/grad.py.log`, `xi.py.log`, and `wvs2.py.log`; metadata, code frequencies and coverage counts are in the same directory. The archive checksum is recorded in `data/WVS_download_manifest.json`.

This confirms the main gradients and cell-regression estimates on the current release. It does not yet verify every robustness exercise, resolve country-specific settlement-boundary differences, validate the area-scaling assumption or establish causal identification. The manuscript has not yet been updated with these new counts and estimates. Raw WVS data must not be redistributed; raw archives and caches in data/ are excluded from Git.
