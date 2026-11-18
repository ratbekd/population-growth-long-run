"""Section 4.5: sensitivity of the taste share varpi-hat (Table 2, panel A) to the top-code on fertility.

The WVS records children ever had (X011) as 0-4 and "5 or more"; the EVS records exact counts up to 8.
Panel A caps both measures at five. This script rebuilds the panel A sample exactly as evs_pool.py does and
reports, for the pooled WVS + EVS sample of respondents aged 45-70:
  1. the share of respondents at or above each cap, by measure and study;
  2. varpi-hat with both measures capped at 5 (baseline), 4 and 3, and with ideal fertility left uncapped;
  3. varpi-hat on the country-waves where fewer than 5 per cent of respondents sit at the top code (X011 >= 5);
  4. interval-regression checks that treat the top category as right-censored, [c, inf), instead of capping:
     a censored Poisson model (same mean function as PPML, so its coefficient is the same elasticity) and a
     Gaussian interval regression in levels (elasticity = coefficient / weighted mean of the outcome).
All estimates: survey weights, controls for age, age^2, sex, religiosity, marital status, country x wave
effects; intervals are percentiles of a bootstrap over country-waves (B=120, seed 1), as in evs_pool.py.
Inputs (not redistributed; see data/README.md): replication/data/WVS_Time_Series_1981-2022_Stata_v5_0.zip and
replication/data/ZA7503_v3-0-0.dta (EVS trend file). The panel A sample follows the supplied evs_pool.py
(replication/original/evs_pool.py); the WVS wave is S002VS, as in prepare_wvs.py.
Output: results/topcode/topcode.log (also printed).
Run: python topcode.py   (about 40 minutes on 14 cores; the 12 specifications run in parallel, one thread each)
"""
import os, zipfile
from pathlib import Path
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')     # one thread per worker: the specifications run in parallel instead
import numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy.special import gammaln, gammainc, log_ndtr
from multiprocessing import Pool

HERE = Path(__file__).resolve().parent
DATA = HERE / 'data'; OUT = HERE / 'results' / 'topcode'


def ppml(y, X, w, it=80, b0=None):
    """Poisson pseudo-maximum likelihood by IRLS; identical to common.ppml of the supplied code."""
    b = np.zeros(X.shape[1]) if b0 is None else b0.copy()
    if b0 is None: b[0] = np.log(max(np.average(y, weights=w), 1e-9))
    for _ in range(it):
        mu = np.exp(np.clip(X@b, -30, 30)); W = w*mu; z = X@b + (y-mu)/mu
        bn = np.linalg.solve((X*W[:, None]).T@X + 1e-8*np.eye(X.shape[1]), (X*W[:, None]).T@z)
        if np.max(np.abs(bn-b)) < 1e-8: return bn
        b = bn
    return b

B = 120
def panel_a():
    """Panel A sample, built as in the supplied evs_pool.py and load_evs.py."""
    wcols = ['S002VS', 'S003', 'S017', 'S020', 'D017', 'F028', 'X001', 'X003', 'X007', 'X011', 'X049']
    with zipfile.ZipFile(DATA / 'WVS_Time_Series_1981-2022_Stata_v5_0.zip') as z:
        with z.open('WVS_Time_Series_1981-2022_stata_v5_0.dta') as f:
            w = pd.read_stata(f, columns=wcols, convert_categoricals=False)
    w = w.rename(columns={'S002VS': 'S002'}).astype(float); w['study'] = 'W'; w['wv'] = w.S002
    ecols = ['S002EVS', 'S003', 'S017', 'S020', 'D017', 'F028', 'X001', 'X003', 'X007', 'X011', 'X049']
    e = pd.read_stata(DATA / 'ZA7503_v3-0-0.dta', columns=ecols, convert_categoricals=False)
    e = e.rename(columns={'S002EVS': 'evsw'}).apply(pd.to_numeric, errors='coerce'); e['study'] = 'E'; e['wv'] = e.evsw
    cols = ['study', 'wv', 'S003', 'S017', 'S020', 'D017', 'F028', 'X001', 'X003', 'X007', 'X011', 'X049']
    d = pd.concat([w[cols], e[cols]], ignore_index=True)
    for c in ['D017', 'X011', 'X049', 'F028', 'X007']: d[c] = d[c].where(d[c] >= 0)
    dup = set(zip(e.S003, e.S020)) & set(zip(w.S003, w.S020))
    d = d[~((d.study == 'E') & pd.Series(list(zip(d.S003, d.S020))).isin(dup).values)]
    d['age'] = d.X003.where(d.X003.between(15, 99)); d['sex'] = d.X001.where(d.X001.isin([1, 2]))
    d['wt'] = d.S017.where(d.S017 > 0, 1.0)
    d['cw'] = d.study + "_" + d.S003.astype(int).astype(str) + "_" + d.wv.astype(int).astype(str)
    mid = {1: 1000, 2: 3162, 3: 7071, 4: 14142, 5: 31623, 6: 70711, 7: 223607, 8: 1000000}
    d['lhab'] = -np.log(d.X049.map(mid) / 1e3)
    s = d[d.D017.notna() & (d.D017 <= 8)]
    a = s[s.age.between(45, 70) & s.X011.notna() & s.lhab.notna()]
    a = a.dropna(subset=['F028', 'X007', 'sex']).reset_index(drop=True)
    return a


def design(g):
    """Covariates with a full set of country x wave dummies and no separate constant (the same model as
    evs_pool.py, which uses a constant and drops one dummy). Built once per sample; a bootstrap draw takes rows,
    and the dummy of a country-wave that is not drawn is then a zero column, which the small ridge in the
    estimators sets to zero without affecting the other coefficients."""
    X = pd.DataFrame({'hab': g.lhab, 'age': g.age, 'age2': g.age**2/100, 'fem': (g.sex == 2)*1.0,
                      'rel': g.F028}, index=g.index)
    X = X.join(pd.get_dummies(g.X007.astype(int), prefix='m', drop_first=True).astype(float))
    X = X.join(pd.get_dummies(g.cw, prefix='cw').astype(float))
    return np.ascontiguousarray(X.values), 0, X.columns.str.startswith('cw_')


# ---- estimators ----
def newton(f, p0, it=100, tol=1e-9):
    """Maximise a concave log-likelihood; f(p) -> (ll, gradient, Hessian). Step-halving line search."""
    p = p0.copy(); ll, g, H = f(p)
    for _ in range(it):
        step = np.linalg.solve(H - 1e-10*np.eye(len(p)), -g); t = 1.0
        while True:
            ln, gn, Hn = f(p + t*step)
            if np.isfinite(ln) and ln >= ll - 1e-12 or t < 1e-6: break
            t /= 2
        p = p + t*step; done = abs(ln - ll) < tol*max(1, abs(ll)) and np.max(np.abs(t*step)) < 1e-7
        ll, g, H = ln, gn, Hn
        if done: break
    return p


def cens_poisson(y, X, w, c, b0):
    """Poisson MLE with y >= c treated as right-censored: P(Y >= c) = P(c, mu), regularised lower gamma."""
    cen = y >= c; yy = np.where(cen, c, y); w = w / w.mean()
    def f(b):
        eta = np.clip(X @ b, -30, 30); mu = np.exp(eta)
        S = np.clip(gammainc(c, mu), 1e-300, 1)
        ll = np.where(cen, np.log(S), yy*eta - mu - gammaln(yy+1))
        p = np.exp((c-1)*eta - mu - gammaln(c))            # Poisson pmf at c-1; dS/dmu = p
        r = p/S
        # d ll/d eta: uncensored y - mu; censored mu p/S.  d2: -mu; mu^2 r [(c-mu)/mu - r]
        g = np.where(cen, mu*r, yy - mu)
        h = np.where(cen, mu*r*((c - mu) - mu*r), -mu)
        return (w*ll).sum(), X.T @ (w*g), (X*(w*h)[:, None]).T @ X
    return newton(f, b0)


def intreg(y, X, w, c, b0):
    """Gaussian interval regression in levels, y >= c right-censored; parameters (beta, log sigma)."""
    cen = y >= c; w = w / w.mean()
    def f(q):
        b, ls = q[:-1], q[-1]; s = np.exp(ls); z = (np.where(cen, c, y) - X @ b)/s
        lnS = log_ndtr(-z); lam = np.exp(-0.5*z**2 - 0.5*np.log(2*np.pi) - lnS)   # phi(z)/(1-Phi(z))
        ll = np.where(cen, lnS, -0.5*z**2 - ls)
        gx = np.where(cen, lam/s, z/s); gs = np.where(cen, lam*z, z**2 - 1)
        d = lam*(lam - z)
        hxx = np.where(cen, -d/s**2, -1/s**2)
        hxs = np.where(cen, (-d*z - lam)/s, -2*z/s)
        hss = np.where(cen, -d*z**2 - lam*z, -2*z**2)
        H = np.empty((len(q), len(q)))
        H[:-1, :-1] = (X*(w*hxx)[:, None]).T @ X; H[:-1, -1] = H[-1, :-1] = X.T @ (w*hxs); H[-1, -1] = (w*hss).sum()
        return (w*ll).sum(), np.r_[X.T @ (w*gx), (w*gs).sum()], H
    return newton(f, b0)


def fit(kind, k, X, j, cwm, yi, yr, wt):
    """(ideal elasticity, realised elasticity, ratio) for one sample. kind: cap | cens | intreg."""
    out = []
    def pp(y):   # PPML started from log mean on the country-wave dummies (there is no constant column)
        return ppml(y, X, wt, b0=np.where(cwm, np.log(max(np.average(y, weights=wt), 1e-9)), 0.0))
    for y in (yi, yr):
        if kind == 'cap':
            out.append(pp(np.minimum(y, k[len(out)]))[j]); continue
        y = np.minimum(y, k)
        if kind == 'cens':
            out.append(cens_poisson(y, X, wt, k, pp(y))[j])
        else:
            sw = np.sqrt(wt)
            b0 = np.r_[np.linalg.lstsq(X*sw[:, None], y*sw, rcond=None)[0], np.log(y.std())]
            out.append(intreg(y, X, wt, k, b0)[j]/np.average(y, weights=wt))
    return [out[0], out[1], out[0]/out[1]]


def run_spec(spec):
    """Point estimate and cluster bootstrap over country-waves; draws identical to common.boot (seed 1)."""
    lab, g, kind, k = spec
    X, j, cwm = design(g); yi, yr, wt = g.D017.values, g.X011.values, g.wt.values
    est = fit(kind, k, X, j, cwm, yi, yr, wt)
    codes, u = pd.factorize(g.cw); idx = [np.where(codes == c)[0] for c in range(len(u))]
    rng = np.random.default_rng(1); bs = []
    for _ in range(B):
        ii = np.concatenate([idx[c] for c in rng.integers(0, len(u), len(u))])
        bs.append(fit(kind, k, X[ii], j, cwm, yi[ii], yr[ii], wt[ii]))
    return lab, len(g), g.cw.nunique(), np.array(est), np.array(bs)


lines = []
def say(t=""):
    print(t, flush=True); lines.append(t)


def show(r):
    lab, n, ncw, est, bs = r; se = bs.std(0); ci = np.percentile(bs[:, 2], [2.5, 97.5])
    say(f"  {lab:44s} N={n:6d} cw={ncw:3d} | ideal {est[0]:+.4f}({se[0]:.4f}) "
        f"realised {est[1]:+.4f}({se[1]:.4f}) varpi {est[2]:.3f} [{ci[0]:.3f},{ci[1]:.3f}]")


def main():
    a = panel_a()

    # ---- 1. top-code shares ----
    say("1. SHARE AT OR ABOVE THE CAP, panel A sample, ages 45-70 (weighted)")
    for lab, g in [("pooled", a), ("WVS", a[a.study == 'W']), ("EVS", a[a.study == 'E'])]:
        say(f"  {lab:7s} N={len(g):6d} | " + "  ".join(
            f"{v}>={k}: {np.average(g[v] >= k, weights=g.wt):.3f}" for v in ['X011', 'D017'] for k in (5, 4, 3)))
    say("  WVS codes X011 as 0-4 and '5 or more'; EVS records exact counts to 8.")
    cwsh = a.groupby('cw').apply(lambda g: pd.Series({'x': np.average(g.X011 >= 5, weights=g.wt),
                                                      'i': np.average(g.D017 >= 5, weights=g.wt)}))
    say(f"  country-waves: {len(cwsh)}; median share X011>=5 {cwsh.x.median():.3f}, max {cwsh.x.max():.3f}; "
        f"median share D017>=5 {cwsh.i.median():.3f}, max {cwsh.i.max():.3f}")
    keep = cwsh.index[cwsh.x < 0.05]; sub = a[a.cw.isin(keep)].reset_index(drop=True)
    keep2 = cwsh.index[(cwsh.x < 0.05) & (cwsh.i < 0.05)]; sub2 = a[a.cw.isin(keep2)].reset_index(drop=True)
    wvs = a[a.study == 'W'].reset_index(drop=True)

    # every specification: (label, sample, estimator, cap or censoring point)
    sections = [
        ("\n2. VARPI-HAT UNDER ALTERNATIVE CAPS (PPML, pooled WVS + EVS)", [], [
            ("both measures capped at 5  (baseline)", a, 'cap', (5, 5)),
            ("both measures capped at 4", a, 'cap', (4, 4)),
            ("both measures capped at 3", a, 'cap', (3, 3)),
            ("ideal uncapped (<=8), realised capped at 5", a, 'cap', (99, 5)),
            ("WVS only, both capped at 5", wvs, 'cap', (5, 5)),
            ("WVS only, both capped at 4", wvs, 'cap', (4, 4))]),
        ("\n3. COUNTRY-WAVES WITH < 5% AT THE TOP CODE, both capped at 5",
         [f"  X011 >= 5 below 5%: {len(keep)} of {len(cwsh)} country-waves "
          f"({sub.loc[sub.study == 'W', 'cw'].nunique()} WVS, {sub.loc[sub.study == 'E', 'cw'].nunique()} EVS)",
          f"  X011 >= 5 and D017 >= 5 both below 5%: {len(keep2)} country-waves"], [
            ("X011 top-code share < 5%", sub, 'cap', (5, 5)),
            ("X011 and D017 top-code shares < 5%", sub2, 'cap', (5, 5))]),
        ("\n4. TOP CATEGORY TREATED AS RIGHT-CENSORED, [c, inf), instead of capped", [], [
            ("censored Poisson, censored at 5", a, 'cens', 5),
            ("censored Poisson, censored at 4", a, 'cens', 4),
            ("censored Poisson, censored at 3", a, 'cens', 3),
            ("Gaussian interval regression, censored at 5", a, 'intreg', 5)]),
    ]
    specs = [sp for _, _, ss in sections for sp in ss]
    with Pool(min(len(specs), max(1, (os.cpu_count() or 2) - 2))) as pool:
        res = dict((r[0], r) for r in pool.imap_unordered(run_spec, specs))
    for head, notes, ss in sections:
        say(head)
        for t in notes: say(t)
        for sp in ss: show(res[sp[0]])

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "topcode.log").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
