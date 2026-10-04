#!/home/alon/anaconda3/envs/rmg_env/bin/python
"""Literal-input Rouse first-contact oracle. This module never loads met.

Exact conditional OU dynamics; adaptive observation intervals are empirically
qualified, including minimum contact step and an independent far-step audit.
All physical inputs are the bytes read and fingerprinted before computation.
"""
import os
for _key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_key] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import itertools
import json
import math
import multiprocessing
from pathlib import Path
import platform
import resource
import sys
import time

import numpy as np
import scipy
from scipy.special import erfc
from scipy.stats import norm

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = Path('/home/alon/runs/i030-met-kernel-reference/rework2/reference')


def fingerprint(data):
    return hashlib.sha256(data).hexdigest()


def read_parameters(path=HERE / 'parameters.json'):
    raw = path.read_bytes()
    digest = fingerprint(raw)
    cfg = json.loads(raw)
    assert cfg['schema'] == 3
    for field in ('D0_m2_s', 'sigma0_m', 'temperature_K',
                  'Rg_squared_per_unit_m2', 'boltzmann_J_K', 'avogadro_mol_inverse'):
        assert math.isfinite(cfg[field]) and cfg[field] > 0
        assert cfg['sources'][field]
    for n, chain in cfg['chains'].items():
        assert math.isclose(chain['D_CM_m2_s'], cfg['D0_m2_s'] / int(n), rel_tol=1e-14)
        if int(n) > 1:
            expected = chain['bond_rms_m']**2 * (int(n)**2 - 1) / (6 * int(n))
            assert math.isclose(expected, chain['Rg_m']**2, rel_tol=1e-14)
    return cfg, digest


def mapped_parameters(cfg, case):
    left, right = [cfg['chains'][str(case[key])] for key in ('units_i', 'units_j')]
    small = min(left['units'], right['units'])
    scale = cfg['sigma0_m'] if small == 1 else cfg['chains'][str(small)]['Rg_m']
    diffusion = cfg['D0_m2_s'] / small
    return {'case': case['name'], 'units_i': left['units'], 'units_j': right['units'],
            'scale_m': scale, 'diffusion_unit_m2_s': diffusion,
            'time_unit_s': scale**2 / diffusion,
            'rate_unit_m3_mol_s': cfg['avogadro_mol_inverse'] * scale * diffusion,
            'D_relative_reduced': (left['D_CM_m2_s'] + right['D_CM_m2_s']) / diffusion,
            'D_local_reduced': 2 * cfg['D0_m2_s'] / diffusion,
            'capture_reduced': cfg['sigma0_m'] / scale,
            'bead_friction_kg_s': cfg['boltzmann_J_K'] * cfg['temperature_K'] / cfg['D0_m2_s'],
            'kernel_Rg_i_m': left['kernel_Rg_m'], 'kernel_Rg_j_m': right['kernel_Rg_m'],
            'floor_active': cfg['sigma0_m'] >= 2 * min(left['kernel_Rg_m'], right['kernel_Rg_m'])}


def mode_coefficients(cfg, mapped, pair_class):
    rates, weights, variances = [], [], []
    sites = {'end/end': ('end', 'end'), 'end/mid': ('end', 'mid'),
             'mid/mid': ('mid', 'mid')}[pair_class]
    bead_d = cfg['D0_m2_s'] / mapped['diffusion_unit_m2_s']
    for index, (key, kind) in enumerate(zip(('units_i', 'units_j'), sites)):
        n = mapped[key]
        if n == 1:
            continue
        p = np.arange(1, n, dtype=float)
        lam = 4 * np.sin(np.pi * p / (2*n))**2
        site = 0 if kind == 'end' else n//2 - 1
        b2 = (cfg['chains'][str(n)]['bond_rms_m'] / mapped['scale_m'])**2
        rates.extend(3 * bead_d * lam / b2)
        weights.extend((-1 if index else 1) * np.sqrt(2/n) * np.cos(np.pi*p*(site+.5)/n))
        variances.extend(b2 / (3 * lam))
    rates, weights, variances = map(np.asarray, (rates, weights, variances))
    assert math.isclose(float(np.sum(weights**2 * variances * rates)) + mapped['D_relative_reduced'],
                        mapped['D_local_reduced'], rel_tol=1e-12)
    return rates, weights, variances


def advance_modes(q, dt, rates, variances, rng):
    """The actual transition used by both trajectories and covariance checks."""
    decay = np.exp(-dt[:, None] * rates)
    noise = np.sqrt(variances * -np.expm1(-2 * dt[:, None] * rates))
    return q * decay[:, :, None] + rng.standard_normal(q.shape) * noise[:, :, None]


def projected_site(q, weights):
    return np.einsum('bpa,p->ba', q, weights, optimize=False)


def minimum_image(vector, side):
    return vector - side * np.floor(vector / side + .5)


def guard_observation_step(q, rates, weights, gap, proposed, dt_min, safety):
    """Return an evaluated passing step, or the declared floor.

    The norm of a sum of differently relaxing modes need not be monotone.
    Halving and re-evaluating requires no monotonicity or linearization premise.
    A final time truncation may make the supplied floor smaller than dt_min.
    """
    floor = np.broadcast_to(dt_min, proposed.shape)
    dt = np.maximum(proposed, floor).copy()
    while True:
        shift = projected_site(q * np.expm1(-dt[:, None]*rates)[:, :, None], weights)
        bad = (dt > floor) & (np.linalg.norm(shift, axis=1) > gap/safety)
        if not np.any(bad):
            return dt
        dt[bad] = np.maximum(dt[bad]/2, floor[bad])


def transient_from_counts(counts, replicas, side, edges):
    """Finite-box conditional hazards; four bins are report-only."""
    survival = np.asarray(counts, dtype=float)/replicas
    bins = []
    for index, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        start, end = counts[index:index+2]
        assert start >= end > 0, 'transient needs survivors'
        rate = side**3/(hi-lo)*math.log(start/end)
        se = side**3/(hi-lo)*math.sqrt((start-end)/(start*end))
        bins.append({'t1':lo, 't2':hi, 'survivors_start':start, 'survivors_end':end,
                     'events':start-end, 'k_reduced':rate, 'se_reduced':se})
    return {'edges':list(edges), 'survivors':list(counts),
            'survival':survival.tolist(),
            'survival_se':np.sqrt(survival*(1-survival)/replicas).tolist(),
            'bins':bins, 'policy':'report-only; no acceptance'}


def hazard(hit, side, window, mapped, cfg):
    lo, hi = window
    start, end = int(np.count_nonzero(hit > lo)), int(np.count_nonzero(hit > hi))
    assert start > end > 0, 'nonzero events and survivors required'
    raw = side**3 / (hi-lo) * math.log(start/end)
    raw_se = side**3 / (hi-lo) * math.sqrt((start-end)/(start*end))
    green = cfg['finite_box_green_constant'] / (4*np.pi*mapped['D_relative_reduced']*side)
    corrected = raw / (1 + green*raw)
    se = raw_se / (1 + green*raw)**2
    return {'window': list(window), 'survivors_start': start, 'survivors_end': end,
            'events': start-end, 'k_box_reduced': raw, 'se_box_reduced': raw_se,
            'k_infinite_reduced': corrected, 'se_infinite_reduced': se,
            'box_removed_fraction': 1-corrected/raw,
            'mixing_times_at_start': lo / (side**2/(4*np.pi**2*mapped['D_relative_reduced']))}


def simulate(run, cfg, mapped, pair_class, seed):
    rates, weights, variances = mode_coefficients(cfg, mapped, pair_class)
    rng = np.random.Generator(np.random.PCG64(seed))
    side = run['box_over_scale']
    # Windows scale with relative COM diffusion so unlike pairs also mix.
    transient = run['name'] == 'transient'
    window = np.asarray(run['window']) * (1 if transient else 2 / mapped['D_relative_reduced'])
    tmax = window[1]
    dtmin, dtmax = run['dt_reduced'], run['max_dt_reduced']
    a = mapped['capture_reduced']
    numerical_a = a + cfg['boundary_shift_constant'] * math.sqrt(2*mapped['D_local_reduced']*dtmin)
    hit = np.full(run['replicas'], np.inf)
    transitions, peak = 0, 0
    for offset in range(0, len(hit), 2048):
        batch = min(2048, len(hit)-offset)
        q = rng.standard_normal((batch, len(rates), 3)) * np.sqrt(variances)[None, :, None]
        internal = projected_site(q, weights)
        initial = rng.uniform(-side/2, side/2, (batch, 3))
        inside = np.sum(initial**2, axis=1) <= a*a
        while np.any(inside):
            initial[inside] = rng.uniform(-side/2, side/2, (int(inside.sum()), 3))
            inside = np.sum(initial**2, axis=1) <= a*a
        com = initial - internal
        clock = np.zeros(batch)
        ids = np.arange(offset, offset+batch)
        while len(ids):
            position = minimum_image(com + projected_site(q, weights), side)
            gap = np.maximum(np.linalg.norm(position, axis=1)-a, 0.)
            dt = np.minimum(dtmax, (gap/run['safety'])**2 / (2*mapped['D_local_reduced']))
            dt = np.minimum(np.maximum(dt, dtmin), tmax-clock)
            dt = guard_observation_step(q, rates, weights, gap, dt,
                                        np.minimum(dtmin, tmax-clock), run['safety'])
            q = advance_modes(q, dt, rates, variances, rng)
            com += rng.standard_normal(com.shape) * np.sqrt(2*mapped['D_relative_reduced']*dt)[:, None]
            clock += dt
            position = minimum_image(com + projected_site(q, weights), side)
            contact = np.sum(position**2, axis=1) <= numerical_a**2
            hit[ids[contact]] = clock[contact]
            transitions += len(ids)
            alive = ~contact & (clock < tmax-1e-12)
            ids, q, com, clock = ids[alive], q[alive], com[alive], clock[alive]
    middle = float(np.mean(window))
    summary = {'case': mapped['case'], 'pair_class': pair_class, 'name': run['name'],
               'seed': seed, 'replicas': len(hit), 'box_over_scale': side,
               'dt_reduced': dtmin, 'max_dt_reduced': dtmax, 'safety': run['safety'],
               'numerical_capture_reduced': numerical_a,
               'hit_times_sha256': fingerprint(np.round(np.where(np.isfinite(hit), hit, -1.), 12).astype('<f8').tobytes()),
               'survival_edges': [0., window[0], middle, window[1]],
               'survivors': [int(np.count_nonzero(hit > t)) for t in [0., window[0], middle, window[1]]],
               'mode_transitions': transitions}
    if transient:
        edges = cfg['transient_edges_reduced']
        counts = [int(np.count_nonzero(hit > t)) for t in edges]
        summary['transient'] = transient_from_counts(counts, len(hit), side, edges)
        summary['survival_edges'], summary['survivors'] = list(edges), counts
    else:
        summary['long'] = hazard(hit, side, window, mapped, cfg)
        summary['late_halves'] = [hazard(hit, side, [window[0], middle], mapped, cfg),
                                  hazard(hit, side, [middle, window[1]], mapped, cfg)]
    return summary


def mean_se(values):
    return float(np.mean(values)), float(np.std(values, ddof=1)/math.sqrt(len(values)))


def sampler_checks(cfg):
    checks = []
    setup = cfg['sampler']
    for case_index, case in enumerate(cfg['cases']):
        mapped = mapped_parameters(cfg, case)
        for class_index, pair in enumerate(cfg['pair_classes']):
            rates, weights, variances = mode_coefficients(cfg, mapped, pair)
            if not len(rates):
                continue
            rng = np.random.Generator(np.random.PCG64(setup['seed'] + 100*case_index + class_index))
            q = rng.standard_normal((setup['replicas'], len(rates), 3)) * np.sqrt(variances)[None, :, None]
            initial = projected_site(q, weights)[:, 0].copy()
            dt = np.full(setup['replicas'], setup['dt_reduced'])
            for step in range(max(setup['lags'])+1):
                if step in setup['lags']:
                    observed, se = mean_se(initial*projected_site(q, weights)[:, 0])
                    theory = float(np.sum(weights**2*variances*np.exp(-rates*step*setup['dt_reduced'])))
                    passed = abs(observed-theory) <= setup['z']*se
                    checks.append({'case':case['name'], 'pair_class':pair, 'lag':step,
                                   'observed_covariance':observed, 'se':se, 'theory':theory, 'pass':passed})
                if step < max(setup['lags']):
                    q = advance_modes(q, dt, rates, variances, rng)
    assert all(item['pass'] for item in checks), 'sampled covariance failed'
    return checks


def independent_bead_contact(cfg, pair, seed):
    """Euler bead force/noise; no mode basis, rates or transition function."""
    setup = cfg['bead_crosscheck']
    n, count, dt, side = setup['N'], setup['replicas'], setup['dt_reduced'], setup['box_over_scale']
    scale = cfg['chains'][str(n)]['Rg_m']
    b2 = (cfg['chains'][str(n)]['bond_rms_m']/scale)**2
    dbead = float(n)
    a = cfg['sigma0_m']/scale
    numerical_a = a + cfg['boundary_shift_constant']*math.sqrt(4*dbead*dt)
    rng = np.random.Generator(np.random.PCG64(seed))
    # Independent equilibrium Gaussian bonds, then subtract COM.
    beads = np.zeros((count, 2, n, 3))
    beads[:, :, 1:] = np.cumsum(rng.normal(0., math.sqrt(b2/3), (count, 2, n-1, 3)), axis=2)
    beads -= beads.mean(axis=2, keepdims=True)
    positions = [0, 0] if pair == 'end/end' else ([0, n//2-1] if pair == 'end/mid' else [n//2-1]*2)
    initial = rng.uniform(-side/2, side/2, (count, 3))
    inside = np.sum(initial**2, axis=1) <= a*a
    while np.any(inside):
        initial[inside] = rng.uniform(-side/2, side/2, (int(inside.sum()), 3))
        inside = np.sum(initial**2, axis=1) <= a*a
    beads[:, 0] += (initial-beads[:, 0, positions[0]]+beads[:, 1, positions[1]])[:, None, :]
    hit = np.zeros(count, dtype=bool)
    for step in range(setup['steps']):
        bonds = beads[:, :, 1:]-beads[:, :, :-1]
        force = np.zeros_like(beads)
        force[:, :, :-1] += bonds
        force[:, :, 1:] -= bonds
        beads += 3*dbead/b2*dt*force + rng.normal(0., math.sqrt(2*dbead*dt), beads.shape)
        sep = minimum_image(beads[:, 0, positions[0]]-beads[:, 1, positions[1]], side)
        hit |= np.sum(sep**2, axis=1) <= numerical_a**2
    return int(hit.sum())


def mode_grid_contact(cfg, pair, seed):
    setup = cfg['bead_crosscheck']
    case = {'name':'small_N_contact', 'units_i':setup['N'], 'units_j':setup['N']}
    mapped = mapped_parameters(cfg, case)
    rates, weights, variances = mode_coefficients(cfg, mapped, pair)
    count, dt, side = setup['replicas'], setup['dt_reduced'], setup['box_over_scale']
    rng = np.random.Generator(np.random.PCG64(seed))
    q = rng.normal(size=(count, len(rates), 3))*np.sqrt(variances)[None, :, None]
    initial = rng.uniform(-side/2, side/2, (count, 3))
    a = mapped['capture_reduced']
    inside = np.sum(initial**2, axis=1) <= a*a
    while np.any(inside):
        initial[inside] = rng.uniform(-side/2, side/2, (int(inside.sum()), 3))
        inside = np.sum(initial**2, axis=1) <= a*a
    com = initial-projected_site(q, weights)
    numerical_a = a+cfg['boundary_shift_constant']*math.sqrt(2*mapped['D_local_reduced']*dt)
    hit = np.zeros(count, dtype=bool)
    steps = np.full(count, dt)
    for step in range(setup['steps']):
        q = advance_modes(q, steps, rates, variances, rng)
        com += rng.normal(0., math.sqrt(2*mapped['D_relative_reduced']*dt), com.shape)
        pos = minimum_image(com+projected_site(q, weights), side)
        hit |= np.sum(pos**2, axis=1) <= numerical_a**2
    return int(hit.sum())


def contact_checks(cfg):
    output = []
    setup = cfg['bead_crosscheck']
    for index, pair in enumerate(cfg['pair_classes']):
        beads = independent_bead_contact(cfg, pair, setup['seed']+index)
        modes = mode_grid_contact(cfg, pair, setup['seed']+100+index)
        size = setup['replicas']
        pb, pm = beads/size, modes/size
        se = math.sqrt((pb*(1-pb)+pm*(1-pm))/size)
        # Euler weak-error scale is quoted separately from statistical overlap.
        euler_scale = 3*setup['N']/(cfg['chains'][str(setup['N'])]['bond_rms_m']/cfg['chains'][str(setup['N'])]['Rg_m'])**2*4*setup['dt_reduced']
        bound = setup['z']*se + euler_scale * max(pb, pm)
        output.append({'pair_class':pair, 'bead_contacts':beads, 'mode_contacts':modes,
                       'replicas':size, 'bead_probability':pb, 'mode_probability':pm,
                       'difference':abs(pb-pm), 'combined_se':se, 'euler_relative_scale':euler_scale,
                       'bound':bound, 'pass':abs(pb-pm) <= bound})
    assert all(row['pass'] for row in output), 'independent bead contact cross-check failed'
    return output


def ewald_constant():
    alpha = 2.
    lattice = np.array([(i,j,k) for i in range(-6,7) for j in range(-6,7)
                        for k in range(-6,7) if (i,j,k) != (0,0,0)], dtype=float)
    radius = np.linalg.norm(lattice, axis=1)
    return -float(np.sum(erfc(alpha*radius)/radius)+np.sum(np.exp(-np.pi**2*radius**2/alpha**2)/(np.pi*radius**2))
                  -np.pi/alpha**2-2*alpha/math.sqrt(np.pi))


def log_upper_difference(left, right, z):
    a, b = left['k_infinite_reduced'], right['k_infinite_reduced']
    se = math.hypot(left['se_infinite_reduced']/a, right['se_infinite_reduced']/b)
    return abs(math.log(a/b))+z*se


def audit_linear_forms():
    """The exact 64-form representation of the observed audit-contrast sum."""
    # Coordinates: base, step, contact, adaptive, box, box2, half1, half2.
    contrasts = np.array([[1,-1,0,0,0,0,0,0], [0,1,-1,0,0,0,0,0],
                          [0,0,1,-1,0,0,0,0], [0,0,0,0,0,0,1,-1]])
    boxes = np.array([[0,0,1,0,-1,0,0,0], [0,0,0,0,1,-1,0,0]])
    return np.asarray([np.asarray(signs)@contrasts + box_sign*box
                       for signs in itertools.product((-1,1), repeat=4)
                       for box in boxes for box_sign in (-1,1)])


def audit_log_covariance(own):
    """Delta-method log-rate covariance, retaining whole/half-window dependence."""
    rates = [own[name]['long'] for name in ('base','step','contact','adaptive','box','box2')]
    halves = own['box2']['late_halves']
    samples = rates + halves
    logs = np.log([r['k_infinite_reduced'] for r in samples])
    covariance = np.diag([(r['se_infinite_reduced']/r['k_infinite_reduced'])**2 for r in samples])
    whole = rates[-1]
    # Disjoint conditional-hazard increments have zero first-order covariance.
    # Whole raw hazard is the duration-weighted mean of the two raw halves.
    full_duration = whole['window'][1]-whole['window'][0]
    raw_whole = whole['k_box_reduced']
    green = (raw_whole/whole['k_infinite_reduced']-1)/raw_whole
    derivative_whole = 1/(raw_whole*(1+green*raw_whole))
    predicted_variance = 0.
    for index, half in enumerate(halves):
        weight = (half['window'][1]-half['window'][0])/full_duration
        raw_half = half['k_box_reduced']
        derivative_half = 1/(raw_half*(1+green*raw_half))
        coefficient = weight*derivative_whole/derivative_half
        cross = coefficient*covariance[6+index,6+index]
        covariance[5,6+index] = covariance[6+index,5] = cross
        predicted_variance += coefficient**2*covariance[6+index,6+index]
    assert math.isclose(predicted_variance, covariance[5,5], rel_tol=2e-12)
    assert np.linalg.eigvalsh(covariance).min() > -1e-14
    return logs, covariance


def joint_numerical_uncertainty(own, z):
    logs, covariance = audit_log_covariance(own)
    forms = audit_linear_forms()
    contrasts = forms@logs
    observed = float(contrasts.max())
    standard_errors = np.sqrt(np.maximum(np.einsum('bi,ij,bj->b',forms,covariance,forms),0))
    tail = 2*norm.sf(z)
    joint_z = float(norm.isf(tail/len(forms)))
    upper = float(np.max(contrasts + joint_z*standard_errors))
    stat = z*math.sqrt(covariance[5,5])
    return {'statistics_log':stat, 'observed_audit_contrasts_log':observed,
            'joint_audit_sampling_log':upper-observed, 'audit_upper_log':upper,
            'numerical_log_band':stat+upper, 'linear_forms':len(forms),
            'joint_normal_quantile':joint_z, 'audit_family_tail_allowance':float(tail),
            'normal_approximation':True, 'log_rate_covariance':covariance.tolist()}


def qualify(rows, cfg, mapped):
    output = []
    z, limits = cfg['statistical_sigmas'], cfg['qualification']
    for case in mapped:
        for pair in cfg['pair_classes']:
            own = {r['name']:r for r in rows if r['case'] == case['case'] and r['pair_class'] == pair}
            ref = own['box2']['long']
            k, se = ref['k_infinite_reduced'], ref['se_infinite_reduced']
            # Individual caps stay separate from the joint acceptance envelope.
            budget = {'stat':z*se/k,
                      'step':log_upper_difference(own['base']['long'], own['step']['long'], z),
                      'contact':log_upper_difference(own['step']['long'],own['contact']['long'],z),
                      'adaptive':log_upper_difference(own['contact']['long'],own['adaptive']['long'],z),
                      'box':max(log_upper_difference(own['contact']['long'],own['box']['long'],z),
                                log_upper_difference(own['box']['long'], own['box2']['long'], z)),
                      'plateau':log_upper_difference(*own['box2']['late_halves'],z)}
            numerical = joint_numerical_uncertainty(own,z)
            total = numerical['numerical_log_band']
            checks = {key: value <= limits['max_'+key+'_log'] for key, value in budget.items()}
            checks['total'] = total <= limits['max_total_log']
            checks['correction'] = ref['box_removed_fraction'] <= limits['max_box_removed_fraction']
            checks['spatial_relaxation'] = all(own[name]['long']['mixing_times_at_start'] >= limits['minimum_mixing_times']
                                               for name in ('contact','box','box2'))
            if case['floor_active'] and case['units_i'] == case['units_j'] == 1:
                sphere = 4*np.pi*case['D_relative_reduced']*case['capture_reduced']
                checks['analytic_sphere'] = abs(math.log(k/sphere)) <= total
            checks = {key:bool(value) for key,value in checks.items()}
            output.append({'case':case['case'], 'pair_class':pair, 'reference_reduced':k,
                           'reference_se_reduced':se, 'rate_unit_m3_mol_s':case['rate_unit_m3_mol_s'],
                           'individual_qualification_log':budget,
                           'numerical_uncertainty':numerical,
                           'legacy_sum_of_marginal_bounds_log':sum(budget.values()),
                           'total_log_tolerance':total,
                           'factor_tolerance':math.exp(total), 'checks':checks,
                           'qualified':all(checks.values())})
    return output


def run_one(run, cfg, mapped, pair, seed):
    begin = time.monotonic()
    print(f"RUN {mapped['case']} {pair} {run['name']} seed={seed}", flush=True)
    row = simulate(run, cfg, mapped, pair, seed)
    if 'long' in row:
        long = row['long']
        value = f"{long['k_infinite_reduced']:.8g} +/- {long['se_infinite_reduced']:.5g}; events={long['events']}"
    else:
        value = f"S(0.2)={row['transient']['survival'][-1]:.8g}; four REPORT-ONLY bins"
    print(f"MEASURED {mapped['case']} {pair} {run['name']}: {value}; {time.monotonic()-begin:.1f}s; RSS={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.1f}MiB", flush=True)
    return row


def pooled_hazard(parts, mapped, cfg):
    first = parts[0]
    lo, hi = map(float,first['window'])
    start = sum(part['survivors_start'] for part in parts)
    end = sum(part['survivors_end'] for part in parts)
    assert start > end > 0
    # Exact same conditional survival estimator, using pooled integer counts.
    side = first['_side']
    raw = side**3/(hi-lo)*math.log(start/end)
    raw_se = side**3/(hi-lo)*math.sqrt((start-end)/(start*end))
    green = cfg['finite_box_green_constant']/(4*np.pi*mapped['D_relative_reduced']*side)
    corrected = raw/(1+green*raw)
    return {'window':[lo,hi], 'survivors_start':start, 'survivors_end':end,
            'events':start-end, 'k_box_reduced':raw, 'se_box_reduced':raw_se,
            'k_infinite_reduced':corrected, 'se_infinite_reduced':raw_se/(1+green*raw)**2,
            'box_removed_fraction':1-corrected/raw,
            'mixing_times_at_start':lo/(side**2/(4*np.pi**2*mapped['D_relative_reduced']))}


def aggregate_stages(rows, cfg, mapped):
    output = []
    for case in mapped:
        for pair in cfg['pair_classes']:
            if case['case']=='floor' and pair!='end/end':
                continue
            for run in cfg['runs']:
                parts = [r for r in rows if r['case']==case['case'] and r['pair_class']==pair and r['name']==run['name']]
                expected = next(item.get('precision_stages',cfg['precision_stages']) for item in cfg['cases'] if item['name']==case['case'])
                assert len(parts)==expected
                parts.sort(key=lambda r:r['precision_stage'])
                row = dict(parts[0])
                row.pop('precision_stage')
                row['replicas'] = sum(r['replicas'] for r in parts)
                row['seeds'] = [r['seed'] for r in parts]
                row['stage_hit_times_sha256'] = [r['hit_times_sha256'] for r in parts]
                # This is a fingerprint of actual stage hashes, never a made-up
                # trajectory or a surrogate first-contact-time array.
                row['hit_times_sha256'] = fingerprint(json.dumps(row['stage_hit_times_sha256']).encode())
                row['survivors'] = [sum(r['survivors'][i] for r in parts) for i in range(len(row['survival_edges']))]
                row['mode_transitions'] = sum(r['mode_transitions'] for r in parts)
                if run['name']=='transient':
                    row['transient'] = transient_from_counts(row['survivors'],row['replicas'],
                                                             row['box_over_scale'],row['survival_edges'])
                else:
                    row['long'] = pooled_hazard([dict(r['long'],_side=r['box_over_scale']) for r in parts],case,cfg)
                    row['late_halves'] = [pooled_hazard([dict(r['late_halves'][i],_side=r['box_over_scale']) for r in parts],case,cfg) for i in range(2)]
                output.append(row)
    floor = [r for r in output if r['case']=='floor']
    output.extend(dict(r,pair_class=pair) for pair in ('end/mid','mid/mid') for r in floor)
    return output


def run_all(cfg, workers, pilot=False, checkpoint=None):
    mapped = [mapped_parameters(cfg, case) for case in cfg['cases']]
    precision = {case['name']:case.get('precision_stages',cfg['precision_stages']) for case in cfg['cases']}
    stages = range(max(precision.values()))
    jobs = [(stage,run,cfg,m,pair,run['seed']+100*ci+pi+100000*stage)
            for stage in stages for ci,m in enumerate(mapped)
            for pi,pair in enumerate(cfg['pair_classes']) for run in cfg['runs']
            if stage < precision[m['case']] and (m['case']!='floor' or pair=='end/end')]
    if pilot:
        job = next(j for j in jobs if j[3]['case']=='N16' and j[4]=='end/end' and j[1]['name']=='base')
        return run_one(dict(job[1],replicas=128),*job[2:])
    rows = []
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
        futures = [(job[0],pool.submit(run_one,*job[1:])) for job in jobs]
        for index,(stage,future) in enumerate(futures):
            row = dict(future.result(),precision_stage=stage)
            rows.append(row)
            if checkpoint is not None:
                temporary = checkpoint / f'{index:03d}.pending'
                temporary.write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
                temporary.rename(checkpoint / f'{index:03d}.json')
    return result_from_stages(rows,cfg)


def result_from_stages(rows,cfg):
    mapped = [mapped_parameters(cfg,case) for case in cfg['cases']]
    pooled = aggregate_stages(rows,cfg,mapped)
    print('SAMPLER: actual OU trajectory covariances at the declared lags',flush=True)
    covariance = sampler_checks(cfg)
    print('CONTACT: independent bead Euler cross-check',flush=True)
    contacts = contact_checks(cfg)
    green = ewald_constant()
    assert abs(green-cfg['finite_box_green_constant']) < 1e-8
    # Use the same pooled-count arithmetic for each stage and the final estimate.
    stages_out = []
    for row in rows:
        case = next(m for m in mapped if m['case']==row['case'])
        row = dict(row)
        if 'long' in row:
            row['long'] = pooled_hazard([dict(row['long'],_side=row['box_over_scale'])],case,cfg)
            row['late_halves'] = [pooled_hazard([dict(h,_side=row['box_over_scale'])],case,cfg) for h in row['late_halves']]
        stages_out.append(row)
    return {'schema':3,'mapped':mapped,'simulation':pooled,'stage_simulation':stages_out,
            'precision_stages':cfg['precision_stages'],'case_precision_stages':{case['name']:case.get('precision_stages',cfg['precision_stages']) for case in cfg['cases']},'sampler_covariance':covariance,
            'bead_contact_crosscheck':contacts,'ewald_constant':green,
            'qualification':qualify(pooled,cfg,mapped)}


def markdown(result):
    lines = ['<!-- BEGIN REPRODUCED NUMBERS -->', '',
             'All +/- values are one Monte Carlo standard error. Individual cap checks and the joint numerical band are empirical finite-setting diagnostics, not rigorous continuum-error bounds.', '',
             '| Case | Pair | Reference / (D_unit length_unit) | Reference (m3 mol-1 s-1) | stat | step cap | contact cap | adaptive cap | box cap | plateau cap | Numerical log | Factor | Qualified |',
             '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|']
    for row in result['qualification']:
        budget = row['individual_qualification_log']
        vals = ' | '.join(f'{budget[key]:.5f}' for key in ('stat','step','contact','adaptive','box','plateau'))
        unit = row['rate_unit_m3_mol_s']
        lines.append(f"| {row['case']} | {row['pair_class']} | {row['reference_reduced']:.9g} +/- {row['reference_se_reduced']:.7g} | {row['reference_reduced']*unit:.9g} +/- {row['reference_se_reduced']*unit:.7g} | {vals} | {row['total_log_tolerance']:.5f} | {row['factor_tolerance']:.5f} | {'YES' if row['qualified'] else 'NO'} |")
    lines += ['', '| Case | Pair | Run | Events | k_box +/- SE | k_infinite +/- SE | Correction removed | Late half 1 | Late half 2 | Mixing times at start |',
              '|---|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for row in result['simulation']:
        if 'long' not in row:
            continue
        r = row['long']; h1,h2 = row['late_halves']
        lines.append(f"| {row['case']} | {row['pair_class']} | {row['name']} | {r['events']} | {r['k_box_reduced']:.8g} +/- {r['se_box_reduced']:.6g} | {r['k_infinite_reduced']:.8g} +/- {r['se_infinite_reduced']:.6g} | {r['box_removed_fraction']:.5f} | {h1['k_infinite_reduced']:.7g} +/- {h1['se_infinite_reduced']:.5g} | {h2['k_infinite_reduced']:.7g} +/- {h2['se_infinite_reduced']:.5g} | {r['mixing_times_at_start']:.5f} |")
    lines += ['', '| Case | Pair | Observed audit sum | Joint audit sampling margin | Primary statistics | Numerical band | Old sum of marginal bounds |',
              '|---|---|---:|---:|---:|---:|---:|']
    for row in result['qualification']:
        n = row['numerical_uncertainty']
        lines.append(f"| {row['case']} | {row['pair_class']} | {n['observed_audit_contrasts_log']:.7g} | {n['joint_audit_sampling_log']:.7g} | {n['statistics_log']:.7g} | {n['numerical_log_band']:.7g} | {row['legacy_sum_of_marginal_bounds_log']:.7g} |")
    lines += ['', 'Dedicated transient ensembles: finite-box conditional hazards, REPORT-ONLY. Times are in each case time unit; no Green correction or acceptance is applied.', '',
              '| Case | Pair | t1 | t2 | S(t2) +/- SE | k_bin reduced +/- SE | k_bin (m3 mol-1 s-1) +/- SE |',
              '|---|---|---:|---:|---:|---:|---:|']
    for row in result['simulation']:
        if 'transient' not in row:
            continue
        unit = next(case['rate_unit_m3_mol_s'] for case in result['mapped'] if case['case']==row['case'])
        transient = row['transient']
        for index, item in enumerate(transient['bins']):
            k,s = item['k_reduced'],item['se_reduced']
            lines.append(f"| {row['case']} | {row['pair_class']} | {item['t1']:g} | {item['t2']:g} | {transient['survival'][index+1]:.8g} +/- {transient['survival_se'][index+1]:.6g} | {k:.8g} +/- {s:.6g} | {k*unit:.8g} +/- {s*unit:.6g} |")
    lines += ['', '| Pair | Bead contacts | Mode contacts | Replicas | Combined SE of probabilities | Allowed difference | Pass |',
              '|---|---:|---:|---:|---:|---:|---|']
    for row in result['bead_contact_crosscheck']:
        lines.append(f"| {row['pair_class']} | {row['bead_contacts']} | {row['mode_contacts']} | {row['replicas']} | {row['combined_se']:.7g} | {row['bound']:.7g} | {row['pass']} |")
    max_z = max(abs(r['observed_covariance']-r['theory'])/r['se'] for r in result['sampler_covariance'])
    lines += ['', f"Sampled covariance: {len(result['sampler_covariance'])} checks; maximum absolute discrepancy / SE = {max_z:.8g} (limit {6:g}).",
              f"Independent Ewald constant: {result['ewald_constant']:.12g}.", '', '<!-- END REPRODUCED NUMBERS -->']
    return '\n'.join(lines)+'\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--workers', type=int, choices=range(1,9), default=8)
    parser.add_argument('--pilot',action='store_true')
    args = parser.parse_args()
    if hasattr(os, 'sched_getaffinity'):
        os.sched_setaffinity(0, set(sorted(os.sched_getaffinity(0))[:8]))
    begin = time.monotonic()
    cfg, parameter_hash = read_parameters()
    source_before = fingerprint(Path(__file__).read_bytes())
    checkpoint = None
    if not args.pilot:
        args.output.mkdir(parents=True,exist_ok=False)
        checkpoint = args.output/'ensembles'
        checkpoint.mkdir()
        (args.output/'parsed_inputs.json').write_text(json.dumps({'parameters_sha256':parameter_hash,'program_sha256':source_before,'parameters':cfg},indent=2,allow_nan=False)+'\n')
    result = run_all(cfg,args.workers,args.pilot,checkpoint)
    if args.pilot:
        print(json.dumps(result, indent=2), flush=True)
        return
    assert fingerprint(Path(__file__).read_bytes()) == source_before, 'oracle changed during run'
    result['parameters_sha256'] = parameter_hash
    result['program_sha256'] = source_before
    result['generation_provenance'] = [{'precision_stage':stage, 'parameters_sha256':parameter_hash, 'program_sha256':source_before, 'driver_program_sha256':source_before, 'cases':[case['name'] for case in cfg['cases'] if stage < case.get('precision_stages',cfg['precision_stages'])]} for stage in range(max(case.get('precision_stages',cfg['precision_stages']) for case in cfg['cases']))]
    result['versions'] = {'python':platform.python_version(), 'numpy':np.__version__, 'scipy':scipy.__version__,
                          'rng':'numpy.PCG64', 'platform':platform.system(), 'machine':platform.machine()}
    (args.output/'results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    (args.output/'results.md').write_text(markdown(result))
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    runtime = {'wall_seconds':time.monotonic()-begin, 'child_cpu_seconds':usage.ru_utime+usage.ru_stime, 'parent_cpu_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime, 'cpu_affinity':sorted(os.sched_getaffinity(0)),
               'child_max_rss_MiB':usage.ru_maxrss/1024,
               'parent_max_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
               'workers':args.workers, 'threads_per_worker':1}
    (args.output/'runtime.json').write_text(json.dumps(runtime, indent=2)+'\n')
    print('REFERENCE COMPUTED:', json.dumps(runtime), flush=True)


if __name__ == '__main__':
    main()
