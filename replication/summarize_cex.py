"""Build a readable results report only after all three panels finish."""
from pathlib import Path
import csv, json
HERE=Path(__file__).resolve().parent
lines=['# CEX replication and inference audit','',
'The supplied archives reproduce the main sample (64,853 consumer units, 139 clusters) and the original point estimates. This audit holds data construction and sample definitions fixed while repairing numerical convergence and removing coefficient-based bootstrap rejection.','',
'Standard errors below use the standard deviation across successful cluster-bootstrap draws (999 attempted, seed 3). Ratio intervals are paired percentile intervals. All failures and retries are retained in the diagnostics; the counts below must be read with the estimates.','']
for panel in 'ABC':
 with (HERE/'results'/f'results_{panel}.csv').open() as f: rows=list(csv.DictReader(f))
 diag=json.loads((HERE/'results'/f'diagnostics_{panel}.json').read_text())
 failed=sum(bool(x['failures']) for x in diag['bootstrap'])
 retries=sum(x['retry_count'] for x in diag['bootstrap'])
 lines += [f'## Panel {panel}','',f"N = {rows[0]['N']}; clusters = {rows[0]['clusters']}. Draws with any failed fit: {failed}; fit retries: {retries}.",'',
'| Outcome | Coefficient | New SE | Successful draws | Ratio to children | 95% ratio interval |',
'|---|---:|---:|---:|---:|---:|']
 for r in rows:
  ratio=f"{float(r['ratio']):.3f}" if r.get('ratio') else '—'
  interval=f"[{float(r['ratio_lo']):.3f}, {float(r['ratio_hi']):.3f}]" if r.get('ratio_lo') else '—'
  lines.append(f"| {r['outcome']} | {float(r['repaired_coefficient']):.5f} | {float(r['bootstrap_se']):.5f} | {r['bootstrap_valid']} | {ratio} | {interval} |")
 lines += ['',f"The legacy routine retained {rows[0]['legacy_retained']} of the first 200 draws after its coefficient filter; {rows[0]['legacy_valid_before_trim']} were numerically returned before that filter. Nonpositive children-coefficient draws: {rows[0]['denominator_nonpositive_draws']}.",'']
 if failed: lines += ['**Inference remains provisional:** inspect the recorded failures and assess their effect before substituting these intervals into the paper.','']
lines += ['## What the paper must disclose','',
'The sample uses reference age 25–44, rounded mean adult count one or two, hourly wage 3–200, mean earners at most two, and positive total and noneducation expenditure. The capped variable is education expenditure divided by hourly wage (99th percentile = 153.3582), not education expenditure alone. The broad panel-C sample omits the wage and expenditure restrictions.','',
'Consumer-unit means across interviews and first nonmissing categorical characteristics are used. Missing spending and hours components are zero-filled in the original construction. These remain substantive measurement choices requiring sensitivity checks. The use of location-cell clusters is reproduced rather than validated as the survey sampling design.','',
'Panel C previously used IQR/1.349 as a scale estimate; this report consistently uses the bootstrap standard deviation. Removing coefficient censoring and changing the number of draws can also alter intervals. The CSV files retain a same-resample legacy comparison for the first 200 draws.','',
'## Limits','',
'This verifies the CEX point estimates and audits their inference; it does not establish causal identification or validate the structural interpretation. WVS, GSS and WPP replication is not complete. WVS respondent-data download awaits authorization to submit registration details and accept use conditions. See README.md for paths and reproducible commands.','']
(HERE/'CEX_replication_report.md').write_text('\n'.join(lines),encoding='utf-8')
print('\n'.join(lines))
