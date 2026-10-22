# Correlation dimension of the target epicentres and its null distribution for a uniform catalogue in the same window.
# Usage: python3 R/correlation_dimension.py ascii_correct.txt      (columns: long lat year month day mag sec hour minute)
# Copyright (C) 2026  Davood Kazemi Lafmejani and Elham Shabani. GPL-3.0-or-later.
import sys, datetime as dt, numpy as np
from scipy.spatial.distance import pdist
a = np.loadtxt(sys.argv[1]); lon, lat, y, mo, d, mag, sec, h, mi = a.T
t0 = dt.datetime(2012, 8, 11)
T = np.array([(dt.datetime(int(y[i]), int(mo[i]), int(d[i]), int(h[i]), int(mi[i])) + dt.timedelta(seconds=float(sec[i])) - t0).total_seconds() / 86400 for i in range(len(a))])
m = (mag >= 3.2) & (lat >= 38.35) & (lat <= 38.5) & (lon >= 46.5) & (lon <= 47) & (T >= 0) & (T <= 4490)
lon, lat = lon[m], lat[m]; n = len(lon); print("target events:", n)
x = (lon - 46.75) * 111.19 * np.cos(np.radians(38.425)); yk = (lat - 38.425) * 111.19
r = np.logspace(-0.5, 1.6, 40)
def slope(x, y, lo, hi):
    dd = pdist(np.c_[x, y]); C = np.array([(dd < ri).mean() for ri in r]); k = (r >= lo) & (r <= hi)
    return np.polyfit(np.log10(r[k]), np.log10(C[k]), 1)[0]
W = 0.5 * 111.19 * np.cos(np.radians(38.425)); H = 0.15 * 111.19; rng = np.random.default_rng(1)
for lo, hi in [(1, 10), (1, 8), (2, 10)]:
    s = [slope(rng.uniform(0, W, n), rng.uniform(0, H, n), lo, hi) for _ in range(500)]
    print("range %g-%g km: D2 = %.2f; uniform null mean %.2f, 95%% range %.2f-%.2f" % (lo, hi, slope(x, yk, lo, hi), np.mean(s), *np.percentile(s, [2.5, 97.5])))
