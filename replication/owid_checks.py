"""Reproduce the results that use the Our World in Data CO2 file.

Usage: python owid_checks.py owid-co2-data.csv

Prints, for the manuscript in revised/:
  - compound emission-intensity rates g_int, 2005-2022 (Tables 5 and D.1);
  - required clean effort m* (Table 5);
  - simulated country paths without a scale effect (Table D.1);
  - the Section 5.2 simulation checks (transition premium, 2143 populations);
  - the asymptotic decline rates for Korea and Japan (Table D.1 note);
  - the sensitivity of long-run population, Table 6 (all rows but quality in costs only), including the
    taste share of the low-fertility European samples, varpi = 0.34, and its required clean effort;
  - the cross-country scale-effect regressions (Table E.1);
  - the sample of Table 3, panel A (country-years and countries).

Requires numpy and pandas only.
"""
import sys

import numpy as np
import pandas as pd

# Calibration (Sections 3.4 and 5.1)
GAMMA, PHI0 = 0.234, 0.138
S_HAT, VARPI = 0.209, 0.539          # unrounded taste share (Table 2, panel A)
E_Q, C0 = 0.068 / 0.289, 4.0           # within-country coefficient ratio, unrounded
G_BAR, THETA, SIGMA = 0.0975, 0.83, 1.0

K = lambda a: a * (1 + GAMMA) / ((a - GAMMA) * (1 + a))
ALPHA_STAR = (PHI0 + GAMMA) / (1 - PHI0)
XI = VARPI * S_HAT / K(ALPHA_STAR)
CHI = (1 - VARPI) * S_HAT
MU = E_Q * (1 + C0) / C0
ZETA_STAR = SIGMA / (1 - THETA)
GY = (1 + G_BAR) ** ZETA_STAR          # output growth per generation at a stationary population
Q = lambda x: ((1 + C0) / (1 + C0 * x)) ** MU
Q_INF = (1 + C0) ** MU
INVERT = lambda n: (PHI0 * n + GAMMA) / (1 - PHI0 * n)

COUNTRIES = {  # name: (ISO code, TFR 2023 from UN WPP 2024)
    "Korea": ("KOR", 0.72), "China": ("CHN", 1.00), "Japan": ("JPN", 1.21),
    "Italy": ("ITA", 1.20), "Sweden": ("SWE", 1.43), "Germany": ("DEU", 1.44),
    "United Kingdom": ("GBR", 1.56), "United States": ("USA", 1.62), "France": ("FRA", 1.64),
}


def bisect(f, lo, hi, tol=1e-12):
    flo = f(lo)
    for _ in range(300):
        mid = (lo * hi) ** 0.5 if lo > 0 and hi / lo > 10 else (lo + hi) / 2
        fm = f(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
        if abs(hi - lo) < tol * max(1.0, abs(mid)):
            break
    return (lo + hi) / 2


def fertility(eps, a_i):
    a = a_i * eps ** XI
    phi = PHI0 * eps ** (-CHI)
    return (a - GAMMA) / (phi * (1 + a))


def eps_star(a_i):
    return bisect(lambda e: fertility(e, a_i) - 1, 1e-3, 1e6)


def simulate(tfr, g_b, m=1.0, nu=0.0, generations=40, tau=0.0, r=2.0):
    """Population path from 2023; habitat per capita and pollution load normalised to one."""
    a_i = INVERT(tfr / 2.1)
    L, x, out = 1.0, 1.0, []
    for t in range(generations + 1):
        eps = Q(x) * (1 + r * tau) / L
        n = fertility(eps, a_i)
        out.append((2023 + 30 * t, L, n))
        x = x * GY * n / (1 + m * g_b * L ** nu)
        L = L * n
    return np.array(out)


def table6(g_b):
    """Table 6: L* = Q_inf / eps* (per cent of 2023) under alternative parameters; g_b by country for m*."""
    global XI, CHI, MU, Q_INF
    saved = XI, CHI, MU, Q_INF
    order = ["Korea", "China", "Japan", "Italy", "Germany", "United States", "France"]
    k = K(ALPHA_STAR)

    def setp(xi=saved[0], chi=saved[1], e_q=E_Q, c0=C0):
        global XI, CHI, MU, Q_INF
        XI, CHI, MU = xi, chi, e_q * (1 + c0) / c0
        Q_INF = (1 + c0) ** MU

    rows = [("Baseline", S_HAT, {}), ("Taste only, xi=0.146", k * 0.146, dict(xi=0.146, chi=0.0)),
            ("Taste share varpi=0.34", S_HAT, dict(xi=0.34 * S_HAT / k, chi=0.66 * S_HAT)),
            ("Feedback doubled", 2 * S_HAT, dict(xi=2 * saved[0], chi=2 * saved[1])),
            ("Pollution elasticity e_Q=1", S_HAT, dict(e_q=1.0)),
            ("Maximum quality gain, c0=1", S_HAT, dict(c0=1.0)), ("Maximum quality gain, c0=9", S_HAT, dict(c0=9.0))]
    print("Table 6: long-run population, per cent of 2023 (" + ", ".join(order) + ")")
    for label, s, kw in rows:
        setp(**kw)
        L = [100 * Q_INF / eps_star(INVERT(COUNTRIES[c][1] / 2.1)) for c in order]
        print(f"  {label:30} S {s:.2f}  " + " ".join(f"{v:5.1f}" for v in L))
    setp(xi=0.34 * S_HAT / k, chi=0.66 * S_HAT)
    print("  varpi=0.34: xi %.3f chi %.3f; required effort m* at nu = 0, 0.5, 1:" % (XI, CHI))
    for name, (iso, tfr) in COUNTRIES.items():
        L_star = Q_INF / eps_star(INVERT(tfr / 2.1))
        m0 = (GY - 1) / g_b[name]
        print(f"    {name:15} L* {L_star:.3f}  m* {m0:.2f} {m0 * L_star ** -0.5:.2f} {m0 / L_star:.2f}")
    XI, CHI, MU, Q_INF = saved


def ols_hc1(y, X):
    X = np.column_stack([np.ones(len(X)), X])
    y = np.asarray(y, dtype=float)
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    e = y - X @ b
    n, k = X.shape
    XtXi = np.linalg.inv(X.T @ X)
    V = XtXi @ (X.T * e ** 2) @ X @ XtXi * n / (n - k)
    return b, np.sqrt(np.diag(V))


def main(path):
    raw = pd.read_csv(path, usecols=["iso_code", "year", "population", "gdp", "co2"])
    owid = raw[raw.iso_code.notna() & ~raw.iso_code.astype(str).str.startswith("OWID")]
    g_b_star = GY - 1

    print("Tables 5 and D.1: compound rates, required effort and simulated paths")
    g_b_all = {}
    for name, (iso, tfr) in COUNTRIES.items():
        d = owid[owid.iso_code == iso].set_index("year")
        intensity = d.co2 / d.gdp
        g_int = (intensity[2022] / intensity[2005]) ** (1 / 17) - 1
        g_b = (1 + g_int) ** -30 - 1
        g_b_all[name] = g_b
        e_star = eps_star(INVERT(tfr / 2.1))
        L_star = Q_INF / e_star
        m0 = g_b_star / g_b
        m_star = {nu: m0 * L_star ** -nu for nu in (0.0, 0.5, 1.0)}
        m_tau = m0 * (1.2 * L_star) ** -0.5
        zeta_hat = -ZETA_STAR * 100 * g_int / 1.84
        path = simulate(tfr, g_b, generations=60)
        years = list(path[:, 0])
        first = next(int(path[i, 0]) for i in range(len(path)) if round(2.1 * path[i, 2], 2) >= 2.00)
        pop = {y: path[years.index(y), 1] for y in (2083, 2143, 2383)}
        print(f"  {name:15} g_int {100 * g_int:6.2f}  eps* {e_star:5.1f}  L* {L_star:.3f}  "
              f"m* {m_star[0.0]:.2f} {m_star[0.5]:.2f} {m_star[1.0]:.2f} {m_tau:.2f}  "
              f"zeta-hat {zeta_hat:4.1f}  TFR>=2.0 by {first}  "
              f"L {pop[2083]:.2f} {pop[2143]:.2f} {pop[2383]:.3f}")

        # Section 5.2 checks with nu = 0.5
        settles = lambda m: abs(simulate(tfr, g_b, m=m, nu=0.5, generations=400)[-1, 1] / L_star - 1) < 0.01
        lo, hi = 0.2 * m_star[0.5], 3 * m_star[0.5]
        for _ in range(40):
            mid = (lo + hi) / 2
            lo, hi = (lo, mid) if settles(mid) else (mid, hi)
        diff_2143 = 100 * (simulate(tfr, g_b, nu=0.5, generations=4)[4, 1] / pop[2143] - 1)
        line = (f"    transition premium {100 * (hi / m_star[0.5] - 1):+.1f}%, "
                f"2143 population nu=0.5 vs nu=0 {diff_2143:+.2f}%")
        if m0 > 1:
            zeta_r = np.log(1 + g_b) / np.log(1 + G_BAR)
            rate = (1 + G_BAR) ** (MU * (zeta_r - ZETA_STAR) / (1 + MU))
            line += f", asymptotic decline {100 * (rate ** (1 / 30) - 1):.3f}% a year"
        print(line)

    table6(g_b_all)

    print("Table E.1: elasticity of clean progress with respect to population")
    for a, b in [(2005, 2022), (1995, 2022)]:
        A, B = owid[owid.year == a].set_index("iso_code"), owid[owid.year == b].set_index("iso_code")
        d = pd.DataFrame({"pa": A.population, "pb": B.population, "ia": A.co2 / A.gdp,
                          "ib": B.co2 / B.gdp, "y": A.gdp / A.population}).dropna()
        d = d[(d.pa > 1e6) & (d.ia > 0) & (d.ib > 0)]
        d["gc"] = (d.ib / d.ia) ** (1 / (b - a)) - 1
        d["gB"] = (1 + d.gc) ** -30 - 1
        d["lp"], d["ly"], d["pg"] = np.log(d.pa), np.log(d.y), np.log(d.pb / d.pa)
        for label, s in [("All", d), ("Richest 40%", d[d.ly >= d.ly.quantile(0.6)])]:
            s1 = s[s.gB > 0]
            b1, se1 = ols_hc1(np.log(s1.gB), s1[["lp", "ly", "pg"]].values)
            b2, se2 = ols_hc1(-100 * s.gc, s[["lp", "ly", "pg"]].values)
            mean = (-100 * s.gc).mean()
            print(f"  {a}-{b} {label:12} (1) {b1[1]:.2f} ({se1[1]:.2f}) N={len(s1)}  "
                  f"(2) {b2[1] / mean:.2f} ({se2[1] / mean:.2f}) N={len(s)}  mean decline {mean:.2f}")

    print("Table 3, panel A sample: country-years with population above one million and positive CO2")
    x = owid.dropna(subset=["gdp", "co2"])
    x = x[(x.population > 1e6) & (x.co2 > 0)]
    for a in (1950, 1990):
        y = x[x.year.between(a, 2022)]
        print(f"  {a}-2022: {len(y)} country-years, {y.iso_code.nunique()} countries")
    print("  countries per year:", x.groupby("year").iso_code.nunique().loc[[1950, 1990, 2022]].to_dict())


if __name__ == "__main__":
    main(sys.argv[1])
