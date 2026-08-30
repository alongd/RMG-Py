"""Emulate the proposed steady-state criterion on a saved RMG profile, and compare
against the variants it is being proposed over (electron-only, L2, no-arming)."""
import sys
import numpy as np

path = sys.argv[1]
TOL = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-6
raw = np.genfromtxt(path, delimiter=',', names=True, deletechars='', replace_space='_')
names = list(raw.dtype.names)
t = raw['Time_(s)']
spc = [n for n in names if n not in ('Time_(s)', 'Volume_(m^3)')]
X = np.array([raw[n] for n in spc])
FLOOR = 1e-16       # the deck's atol; moles ~= mole fraction here (N_total ~ 1 mol)


def slopes(k):
    """|d ln x_i / d ln t| at step k, per species; nan where not live."""
    out = np.full(len(spc), np.nan)
    if k == 0 or t[k - 1] <= 0.0:
        return out
    dlnt = np.log(t[k]) - np.log(t[k - 1])
    if dlnt <= 0:
        return out
    for i in range(len(spc)):
        a, b = X[i, k - 1], X[i, k]
        if a > FLOOR and b > FLOOR:
            out[i] = abs(np.log(b) - np.log(a)) / dlnt
        elif a <= FLOOR and b > FLOOR:
            out[i] = np.inf          # species just appeared: definitely not steady
    return out


def report(label, series, arm):
    """First step at which `series` < TOL, optionally after it first reached >= 1."""
    armed = not arm
    for k in range(len(t)):
        r = series[k]
        if np.isnan(r):
            continue
        if r >= 1.0:
            armed = True
        if armed and r < TOL:
            print('  %-26s fires at t = %.4e s   (R = %.3e, step %d/%d)'
                  % (label, t[k], r, k, len(t) - 1))
            return
    print('  %-26s never fires' % label)


n = len(t)
r_max = np.full(n, np.nan)
r_l2 = np.full(n, np.nan)
r_e = np.full(n, np.nan)
e_idx = [i for i, s in enumerate(spc) if s.startswith('e-')][0]
for k in range(n):
    s = slopes(k)
    live = ~np.isnan(s)
    if live.any():
        r_max[k] = np.nanmax(s)
        r_l2[k] = np.sqrt(np.nansum(s[live] ** 2))
    r_e[k] = s[e_idx]

print('file  :', path)
print('tol   :', TOL)
print()
print('Which variant fires when:')
report('max over species + ARM', r_max, True)
report('max over species, NO arm', r_max, False)
report('L2 over species + ARM', r_l2, True)
report('electron only + ARM', r_e, True)
print()
print('R_max peak = %.3e at t = %.4e s' % (np.nanmax(r_max), t[int(np.nanargmax(r_max))]))
print()
print('R_max around the cliff:')
for k in range(n):
    if not np.isnan(r_max[k]) and 1e-5 < t[k] < 1e-2:
        i = int(np.nanargmax(np.where(np.isnan(slopes(k)), -1, slopes(k))))
        if r_max[k] < 1e2:
            print('   t=%.5e  R=%.4e  worst=%s' % (t[k], r_max[k], spc[i]))
