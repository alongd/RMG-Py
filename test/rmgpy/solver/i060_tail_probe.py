#!/usr/bin/env python3
"""I-060 probe: the hybrid handshake across the tail-exhaustion boundary.

Usage: python i060_tail_probe.py <out.json> [states_to_replay.json]

Emits ONE json document to the named file. Run it against two different builds (via
PYTHONPATH, from a cwd that is not a worktree) and diff the documents; nothing
here interprets its own result, so the same probe serves both arms.

Two measurements:

  "field"  -- the residual evaluated on a grid of tail states at fixed mu0 and
              fixed PDI, sweeping the tail mean DP across the handshake guard's
              threshold xs. Reports the handshake flux (dn_dt at the explicit
              DP=xs species, which NOTHING else in this fixture feeds) and the
              moment derivatives.

  "traj"   -- a real trajectory of the same fixture, integrated with LSODA on
              the solver's own residual, from a fully realizable tail
              (mean DP >> xs) down through the guard's threshold. Reports the
              cone quantities at every step.

The fixture is a single pool with an explicit DP=xs oligomer, k_unzip > 0 (the
handshake's per-chain frequency) and k_scission > 0 (what drives the mean DP
down -- unzip alone is stationary at the tail boundary and never crosses it).
"""
import json
import sys

import numpy as np

from rmgpy.molecule import Molecule
from rmgpy.species import Species
import rmgpy.solver.polymer as polymer_mod
from rmgpy.solver.polymer import HybridPolymerSystem, PolymerPoolConfig

XS = 3
K_UNZIP = 0.1
K_SCISSION = 0.02
IDX_INERT, IDX_PXS, IDX_MU0, IDX_MU1, IDX_MU2, IDX_MON = range(6)


def _spc(smiles, label):
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.label = label
    return s


def build_system(mu0, mu1, mu2, k_scission=K_SCISSION, k_unzip=K_UNZIP):
    inert = _spc("N#N", "N2")
    pxs = _spc("CCC", "P3")          # explicit DP=xs oligomer (condensed)
    m0 = _spc("CO", "Mu0")
    m1 = _spc("C=O", "Mu1")
    m2 = _spc("C#N", "Mu2")
    mon = _spc("C", "M")             # released monomer (gas)
    core = [inert, pxs, m0, m1, m2, mon]
    mask = np.array([True, False, False, False, False, True], dtype=bool)
    pool = PolymerPoolConfig(
        label="poly",
        xs=XS,
        explicit_dp_to_species_index={XS: IDX_PXS},
        mu_indices=(IDX_MU0, IDX_MU1, IDX_MU2),
        monomer_poly_index=IDX_MON,
        k_scission=k_scission,
        k_unzip=k_unzip,
        tail_kinetics=None,
    )
    rs = HybridPolymerSystem(
        T=800.0, P=1.0e5,
        initial_mole_fractions={inert: 1.0},
        V_poly=1.0,
        polymer_pools=[pool],
        mass_transfer=[],
        gas_species_mask=mask,
        constant_gas_volume=False,
        initial_polymer_moments={"poly": (mu0, mu1, mu2)},
        initial_explicit_species={"poly": {XS: 0.0}},
        termination=[],
    )
    rs.initialize_model(core, [], [], [])
    return rs


def _rhs(rs):
    z = np.zeros_like(rs.y)

    def f(t, y):
        return np.asarray(rs.residual(t, np.asarray(y, dtype=float), z)[0],
                          dtype=float)
    return f


def field_scan(pdi=1.5, mu0=1.0, means=None):
    """Residual on a grid of tail states, sweeping mean DP across xs."""
    rs = build_system(1.0, 20.0, 600.0)
    f = _rhs(rs)
    if means is None:
        means = [XS + 4.0, XS + 3.0, XS + 2.0, XS + 1.5, XS + 1.0, XS + 0.5,
                 XS + 1e-6, XS, XS - 1e-6, XS - 0.5, XS - 1.0, XS - 1.5,
                 XS - 2.0, 1.5, 1.0]
    rows = []
    for mean in means:
        mu1 = mu0 * mean
        mu2 = pdi * mu1 * mu1 / mu0
        y = rs.y.copy()
        y[IDX_MU0], y[IDX_MU1], y[IDX_MU2] = mu0, mu1, mu2
        y[IDX_PXS] = 0.0
        d = f(0.0, y)
        rows.append(dict(
            mean=mean, mu0=mu0, mu1=mu1, mu2=mu2,
            handshake=float(d[IDX_PXS]).hex(),
            dmu0=float(d[IDX_MU0]).hex(),
            dmu1=float(d[IDX_MU1]).hex(),
            dmu2=float(d[IDX_MU2]).hex(),
            gas=float(d[IDX_MON]).hex(),
            # drift of the tail's own support invariant mu1 >= (xs+1)*mu0
            de_tail=float(d[IDX_MU1] - (XS + 1) * d[IDX_MU0]).hex(),
        ))
    return rows


def q_scan(seed=20260904, n_rand=4000, k_unzip=K_UNZIP):
    """Round 2: the OTHER half of the cone, Q := mu0*mu2 - mu1^2 >= 0.

    A handshake event removes one chain at exactly n = xs, so with moments
    taken about xs (m_k = SUM (n-xs)^k c_n) it moves m_0 by -F and leaves m_1
    and m_2 alone, giving dQ/dt = -F*m_2 <= 0 for every realizable state. The
    budget that keeps Q >= 0 is c_xs <= Q/m_2 (Cauchy-Schwarz on the measure
    with the n = xs atom removed), and it is TIGHT.

    Reports, per state: the attempted chain removal N = F/k against that
    budget, and the ratio F*m_2/Q, which the budget caps at exactly k. A
    hand-picked PDI x mean grid plus a seeded random sweep; low PDI is what
    is adversarial for Q (Q = mu0^2 * Var, so the budget vanishes as the
    distribution narrows), which is a different corner from the one that is
    adversarial for mu1 - mu0.
    """
    rs = build_system(1.0, 20.0, 600.0, k_scission=0.0, k_unzip=k_unzip)
    f = _rhs(rs)

    def one(mu0, mean, pdi, tag):
        mu1 = mu0 * mean
        mu2 = pdi * mu1 * mu1 / mu0
        y = rs.y.copy()
        y[IDX_MU0], y[IDX_MU1], y[IDX_MU2] = mu0, mu1, mu2
        y[IDX_PXS] = 0.0
        d = f(0.0, y)
        F = float(d[IDX_PXS])
        Q = mu0 * mu2 - mu1 * mu1
        m2 = mu2 - 2.0 * XS * mu1 + XS * XS * mu0
        dQ = float(mu2 * d[IDX_MU0] + mu0 * d[IDX_MU2] - 2.0 * mu1 * d[IDX_MU1])
        return dict(tag=tag, mean=mean, pdi=pdi, mu0=mu0, mu1=mu1, mu2=mu2,
                    F=F, Q=Q, m2=m2, dQ_total=dQ,
                    N=F / k_unzip,
                    B_Q=(Q / m2) if m2 > 0.0 else float("inf"),
                    B_excess=((mu1 - mu0) / (XS - 1)) if XS > 1 else float("inf"))

    rows = []
    for pdi in (1.000000102, 1.000002, 1.001, 1.01, 1.05, 1.2, 1.5, 3.0):
        for mean in (1.2, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, XS + 1.5, 6.0, 12.0, 40.0):
            rows.append(one(1.0, mean, pdi, "grid"))
    rng = np.random.RandomState(seed)
    for _ in range(n_rand):
        mean = float(np.exp(rng.uniform(np.log(1.0001), np.log(60.0))))
        pdi = 1.0 + float(np.exp(rng.uniform(np.log(1e-7), np.log(5.0))))
        rows.append(one(1.0, mean, pdi, "rand"))
    return rows


def q_trajectory(t_end=2.0e-4, n_out=201):
    """The state q_scan's search names as the worst realizable case for Q,
    integrated. mean DP = 4.503939 is xs + 1.5 -- ABOVE the cutoff, so the
    removed boolean was TRUE here and this trajectory is untouched by the
    round-1 change. k_scission = 0 so the handshake is the only channel that
    can drain Q, and nothing else is in the way of reading the result.
    """
    from scipy.integrate import solve_ivp
    mu0, mean, pdi = 1.0, 4.503939, 1.000000102
    mu1 = mu0 * mean
    mu2 = pdi * mu1 * mu1 / mu0
    rs = build_system(mu0, mu1, mu2, k_scission=0.0, k_unzip=K_UNZIP)
    f = _rhs(rs)
    sol = solve_ivp(f, (0.0, t_end), rs.y.copy(), method="LSODA",
                    t_eval=np.linspace(0.0, t_end, n_out),
                    rtol=1e-12, atol=1e-18, max_step=t_end / 400.0)
    rows = []
    for i, t in enumerate(sol.t):
        y = sol.y[:, i]
        a, b, c = float(y[IDX_MU0]), float(y[IDX_MU1]), float(y[IDX_MU2])
        d = f(t, y)
        rows.append(dict(t=float(t), mu0=a, mu1=b, mu2=c,
                         mean=b / a if a > 0 else float("nan"),
                         Q=a * c - b * b, e_1=b - a,
                         handshake=float(d[IDX_PXS]),
                         dQ=float(c * d[IDX_MU0] + a * d[IDX_MU2]
                                  - 2.0 * b * d[IDX_MU1])))
    return dict(status=int(sol.status), message=str(sol.message),
                mu=[mu0, mu1, mu2], rows=rows)


def _row(f, mu0, mu1, mu2, rs, tag):
    y = rs.y.copy()
    y[IDX_MU0], y[IDX_MU1], y[IDX_MU2] = mu0, mu1, mu2
    y[IDX_PXS] = 0.0
    d = f(0.0, y)
    return dict(
        tag=tag, mu0=mu0, mu1=mu1, mu2=mu2,
        mean=mu1 / mu0,
        handshake=float(d[IDX_PXS]).hex(),
        dmu0=float(d[IDX_MU0]).hex(),
        dmu1=float(d[IDX_MU1]).hex(),
        dmu2=float(d[IDX_MU2]).hex(),
        gas=float(d[IDX_MON]).hex(),
        de_tail=float(d[IDX_MU1] - (XS + 1) * d[IDX_MU0]).hex(),
    )


def random_states(seed=20260904, n_cone=120, n_exh=120):
    """Two seeded state sets evaluated by BOTH arms, so every comparison is
    state-matched rather than trajectory-matched.

      'cone' -- states strictly inside the tail's own realizable region,
                mean DP > xs + 1e-9 with PDI > 1, i.e. exactly where the
                removed boolean was TRUE. Nothing here may change by one
                ULP.
      'exh'  -- the exhausting regime, mean DP in (1, xs], where the boolean
                switched the outlet off.
    """
    rs = build_system(1.0, 20.0, 600.0)
    f = _rhs(rs)
    rng = np.random.RandomState(seed)
    rows = []
    for _ in range(n_cone):
        mu0 = float(10.0 ** rng.uniform(-6.0, 2.0))
        mean = float(XS + 1e-9 + 10.0 ** rng.uniform(-6.0, 3.0))
        pdi = float(1.0 + 10.0 ** rng.uniform(-4.0, 0.7))
        mu1 = mu0 * mean
        rows.append(_row(f, mu0, mu1, pdi * mu1 * mu1 / mu0, rs, "cone"))
    for _ in range(n_exh):
        mu0 = float(10.0 ** rng.uniform(-6.0, 2.0))
        mean = float(rng.uniform(1.0, XS))
        pdi = float(1.0 + 10.0 ** rng.uniform(-4.0, 0.7))
        mu1 = mu0 * mean
        rows.append(_row(f, mu0, mu1, pdi * mu1 * mu1 / mu0, rs, "exh"))
    return rows


def replay(states):
    """Evaluate the residual on a list of full state vectors handed in from
    outside -- how the UNFIXED build is measured on the FIXED build's own
    trajectory, so the two-arm comparison is per-step and state-matched."""
    rs = build_system(1.0, 20.0, 600.0)
    f = _rhs(rs)
    rows = []
    for st in states:
        y = np.asarray(st, dtype=float)
        d = f(0.0, y)
        mu0, mu1, mu2 = float(y[IDX_MU0]), float(y[IDX_MU1]), float(y[IDX_MU2])
        rows.append(dict(
            mu0=mu0, mu1=mu1, mu2=mu2,
            mean=(mu1 / mu0) if mu0 > 0 else float("nan"),
            handshake=float(d[IDX_PXS]).hex(),
            dmu0=float(d[IDX_MU0]).hex(),
            dmu1=float(d[IDX_MU1]).hex(),
            dmu2=float(d[IDX_MU2]).hex(),
            gas=float(d[IDX_MON]).hex(),
            de_tail=float(d[IDX_MU1] - (XS + 1) * d[IDX_MU0]).hex(),
        ))
    return rows


def trajectory(t_end=600.0, n_out=601, k_scission=K_SCISSION,
               k_unzip=K_UNZIP, mu=(1.0, 20.0, 600.0)):
    from scipy.integrate import solve_ivp
    rs = build_system(*mu, k_scission=k_scission, k_unzip=k_unzip)
    f = _rhs(rs)
    y0 = rs.y.copy()
    ts = np.linspace(0.0, t_end, n_out)
    sol = solve_ivp(f, (0.0, t_end), y0, method="LSODA", t_eval=ts,
                    rtol=1e-10, atol=1e-14, max_step=t_end / 400.0)
    rows = []
    for i, t in enumerate(sol.t):
        y = sol.y[:, i]
        mu0, mu1, mu2 = float(y[IDX_MU0]), float(y[IDX_MU1]), float(y[IDX_MU2])
        d = f(t, y)
        rows.append(dict(
            t=float(t), mu0=mu0, mu1=mu1, mu2=mu2,
            mean=(mu1 / mu0) if mu0 > 0 else float("nan"),
            pxs=float(y[IDX_PXS]), gas=float(y[IDX_MON]),
            handshake=float(d[IDX_PXS]),
            dmu0=float(d[IDX_MU0]), dmu1=float(d[IDX_MU1]),
            e_tail=mu1 - (XS + 1) * mu0,   # tail cone   mu1 >= (xs+1) mu0
            e_xs=mu1 - XS * mu0,           # brief's candidate
            e_1=mu1 - mu0,                 # generic moment cone
            var=mu0 * mu2 - mu1 * mu1,     # Cauchy-Schwarz
        ))
    # every 10th accepted output state, verbatim, for the other arm to replay
    states = [[float(v) for v in sol.y[:, i]] for i in range(0, len(sol.t), 10)]
    return dict(status=int(sol.status), message=str(sol.message), rows=rows,
                states=states, t_states=[float(sol.t[i]) for i in
                                         range(0, len(sol.t), 10)])


def main():
    out = dict(
        build=polymer_mod.__file__,
        xs=XS, k_unzip=K_UNZIP, k_scission=K_SCISSION,
        tail_conc_min=float(polymer_mod.TAIL_CONC_MIN),
        field=field_scan(),
        # PDI == 1 exactly: _gamma_params_from_mu012 refuses (pdi <= 1+1e-6),
        # so p_cond comes from the in-block monodisperse triangle fallback
        # instead of the gamma leg. Second instance of the same guard shape.
        field_mono=field_scan(pdi=1.0),
        random=random_states(),
        traj=trajectory(),
        q_field=q_scan(),
        q_traj=q_trajectory(),
    )
    if len(sys.argv) > 2:
        with open(sys.argv[2]) as fh:
            out["replay"] = replay(json.load(fh))
    # To a FILE, never stdout: initialize_model prints a solver diagnostic
    # banner to stdout, which would corrupt the document.
    with open(sys.argv[1], "w") as fh:
        json.dump(out, fh, indent=1)
        fh.write("\n")


if __name__ == "__main__":
    main()
