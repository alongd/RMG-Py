"""Offline thermodynamic diagnostics; no RMG product model is changed.

Command: PYTHONPATH=$PWD /home/alon/anaconda3/envs/rmg_env/bin/python
test/rmgpy/kmc/fixtures/i042_probe/models.py

PC-SAFT equations: Gross & Sadowski, doi:10.1021/ie0003887.
Universal constants and independent component regressions: FeOs 0.2.0,
https://github.com/feos-org/feos-pcsaft/tree/main/src/eos
UNIFAC constants: https://www.ddbst.com/published-parameters-unifac.html
Equations are implemented here from their mathematical definitions.
"""

from __future__ import annotations

import math
import numpy as np

R = 8.314472
NA = 6.02214179e23
MW = 104.152  # g/mol, C8H8; named literature convention, not a fitted density
RHO = 1050.0  # kg/m3, dispatch volume convention, held constant with T
C0 = 1000.0  # mol/m3
DP = 100000.0 / MW  # illustrative monodisperse carrier, not measured MWD
CTOT = RHO * 1000.0 / MW  # conserved styrene equivalents per m3

# m per (g/mol), sigma in Angstrom, epsilon/k in K; association absent.
# A: Lopez-Dominguez et al. 2023, Table 2, doi:10.1002/cjce.24908.
# B: Aerts, Residual monomer reduction in polymer latex products by
# extraction with supercritical carbon dioxide, TU/e thesis 2012, Ch.2 Table 2,
# https://pure.tue.nl/ws/files/3406427/723144.pdf (printed p.15).
PC_SETS = {
    "PC-SAFT A": {"monomer_m": 0.023 * MW, "monomer_sigma": 3.78,
                  "monomer_epsilon": 297.66, "polymer_m_per_mass": 0.019,
                  "polymer_sigma": 4.1071, "polymer_epsilon": 267.0, "kij": 0.005},
    "PC-SAFT B": {"monomer_m": 3.08, "monomer_sigma": 3.712,
                  "monomer_epsilon": 295.98, "polymer_m_per_mass": 0.019,
                  "polymer_sigma": 4.107, "polymer_epsilon": 267.0, "kij": 0.0},
}

A = np.array([
    [0.91056314451539, 0.63612814494991, 2.68613478913903, -26.5473624914884,
     97.7592087835073, -159.591540865600, 91.2977740839123],
    [-0.30840169182720, 0.18605311591713, -2.50300472586548, 21.4197936296668,
     -65.2558853303492, 83.3186804808856, -33.7469229297323],
    [-0.09061483509767, 0.45278428063920, 0.59627007280101, -1.72418291311787,
     -4.13021125311661, 13.7766318697211, -8.67284703679646],
])
B = np.array([
    [0.72409469413165, 2.23827918609380, -4.00258494846342, -21.00357681484648,
     26.8556413626615, 206.5513384066188, -355.60235612207947],
    [-0.57554980753450, 0.69950955214436, 3.89256733895307, -17.21547164777212,
     192.6722644652495, -161.8264616487648, -165.2076934555607],
    [0.09768831158356, -0.25575749816100, -9.15585615297321, 20.64207597439724,
     -38.80443005206285, 93.6267740770146, -29.66690558514725],
])


def helmholtz_density(t, c, m, sigma, epsilon, kij, components=False):
    """A_res/(RT V) in mol/m3; c is a vector of molecular concentrations."""
    c, m, sigma, epsilon = map(np.asarray, (c, m, sigma, epsilon))
    # Work in molecule/Angstrom3 to apply the published segment diameters.
    rho = c * NA * 1e-30
    d = sigma * (1 - 0.12 * np.exp(-3 * epsilon / t))
    z = np.array([np.pi / 6 * np.sum(rho * m * d**k) for k in range(4)])
    eta = z[3]
    if not 0 < eta.real < 1:
        raise ValueError("PC-SAFT packing fraction outside (0,1)")
    x = rho / np.sum(rho)
    mb = np.sum(x * m)
    # BMCSL hard spheres, plus Wertheim chain connectivity.
    hs = 6 / np.pi * (3*z[1]*z[2]/(1-eta) + z[2]**3/(eta*(1-eta)**2)
                      + (z[2]**3/eta**2-z[0])*np.log1p(-eta))
    g = 1/(1-eta) + 1.5*d*z[2]/(1-eta)**2 + 0.5*d**2*z[2]**2/(1-eta)**3
    chain = -np.sum(rho * (m-1) * np.log(g))
    weights = np.array([1, (mb-1)/mb, (mb-1)*(mb-2)/mb**2])
    i1 = np.sum((weights @ A) * eta**np.arange(7))
    i2 = np.sum((weights @ B) * eta**np.arange(7))
    c1 = 1/(1 + mb*(8*eta-2*eta**2)/(1-eta)**4
            + (1-mb)*(20*eta-27*eta**2+12*eta**3-2*eta**4)/((1-eta)*(2-eta))**2)
    eij = np.sqrt(epsilon[:, None]*epsilon[None, :])*(1-np.asarray(kij))/t
    sij = (sigma[:, None]+sigma[None, :])/2
    pair = (rho*m)[:, None]*(rho*m)[None, :]*sij**3
    disp = -2*np.pi*i1*np.sum(pair*eij) - np.pi*mb*c1*i2*np.sum(pair*eij**2)
    scale = 1/(NA*1e-30)
    values = np.array([hs, chain, disp])*scale
    return values if components else np.sum(values)


def chemical_potentials(t, c, m, sigma, epsilon, kij):
    """Residual mu/RT via a complex-step concentration derivative at fixed V."""
    result = []
    for i in range(len(c)):
        trial = np.array(c, dtype=complex)
        trial[i] += 1e-20j
        result.append(helmholtz_density(t, trial, m, sigma, epsilon, kij).imag/1e-20)
    return np.array(result)


def pc_state(label, masses, monomer=C0, density=RHO, dp=DP):
    """Carrier + monomer + trace Rn + trace Rn+1; radicals use PS mass scaling."""
    p = PC_SETS[label]
    total = density*1000/MW
    c = np.array([(total-monomer)/dp, monomer, 0.0, 0.0])
    if c[0] <= 0 or monomer <= 0:
        raise ValueError("carrier and monomer must be positive")
    m = np.array([p["polymer_m_per_mass"]*dp*MW, p["monomer_m"],
                  *(p["polymer_m_per_mass"]*mass for mass in masses)])
    sig = np.array([p["polymer_sigma"], p["monomer_sigma"], *[p["polymer_sigma"]]*2])
    eps = np.array([p["polymer_epsilon"], p["monomer_epsilon"], *[p["polymer_epsilon"]]*2])
    kij = np.zeros((4, 4))
    kij[1, [0, 2, 3]] = p["kij"]
    kij[[0, 2, 3], 1] = p["kij"]
    return c, m, sig, eps, kij


def pc_values(t, state):
    c, m, sig, eps, kij = state
    mu = chemical_potentials(t, *state)
    f = helmholtz_density(t, *state)
    pressure = R*t*(sum(c)+np.dot(c, mu)-f)
    return {"G": float(R*t*(mu[3]-mu[2]-mu[1])),
            "mu_J_mol": (R*t*mu).tolist(), "P_Pa": float(pressure)}


def derivative(function, x, step):
    """Fourth-order centered derivative."""
    return (function(x-2*step)-8*function(x-step)+8*function(x+step)-function(x+2*step))/(12*step)


def pc_transfer(t, state):
    result = pc_values(t, state)
    gt = derivative(lambda u: pc_values(u, state)["G"], t, 0.05)
    pt = derivative(lambda u: pc_values(u, state)["P_Pa"], t, 0.05)
    def scaled(log_scale):
        return (state[0]*np.exp(log_scale), *state[1:])
    gl = derivative(lambda u: pc_values(t, scaled(u))["G"], 0, 1e-4)
    pl = derivative(lambda u: pc_values(t, scaled(u))["P_Pa"], 0, 1e-4)
    s_path = -gt
    # These are derivatives of equal-concentration transfer, not fixed-pressure
    # fugacity coefficients. The gas standard remains fixed at c0.
    s_pressure = -gt + gl*pt/pl
    result.update(S_path=s_path, H_path=result["G"]+t*s_path,
                  S_pressure=s_pressure, H_pressure=result["G"]+t*s_pressure,
                  dP_dlogrho_Pa=pl)
    return result


# DDBST original UNIFAC subgroups: CH3, CH2, CH, CH2=CH, ACH, AC, ACCH2, ACCH.
SUBGROUPS = [1, 2, 3, 5, 9, 10, 12, 13]
UR = np.array([0.9011, 0.6744, 0.4469, 1.3454, 0.5313, 0.3652, 1.0396, 0.8121])
UQ = np.array([0.848, 0.540, 0.228, 1.176, 0.400, 0.120, 0.660, 0.348])
MAIN = np.array([1, 1, 1, 2, 3, 3, 4, 4])-1
INTERACTION = np.array([[0, 86.02, 61.13, 76.50], [-35.36, 0, 38.81, 74.15],
                        [-11.12, 3.446, 0, 167.0], [-69.7, -113.6, -146.8, 0]])


def unifac_log_gamma(t, groups, x):
    """Original UNIFAC mole-fraction gamma; trace fractions may be zero."""
    groups, x = np.asarray(groups), np.asarray(x)
    r, q = groups @ UR, groups @ UQ
    rm, qm = x @ r, x @ q
    ell = 5*(r-q)-(r-1)  # coordination number z=10
    comb = np.log(r/rm)+5*q*np.log(q*rm/(r*qm))+ell-r/rm*(x @ ell)
    psi = np.exp(-INTERACTION[MAIN[:, None], MAIN[None, :]]/t)
    def log_group_gamma(counts):
        theta = counts*UQ / np.sum(counts*UQ)
        denominator = theta @ psi
        return UQ*(1-np.log(denominator)-psi @ (theta/denominator))
    mixture = log_group_gamma(x @ groups)
    residual = np.array([np.dot(g, mixture-log_group_gamma(g)) for g in groups])
    return comb+residual


def unifac_state(monomer=C0, dp=DP):
    # End-saturated actual I039 proxies: one ACCH2 benzyl end, one CH3 end.
    # Bath is DP copies of CH2 + ACCH + 5 ACH, no finite-chain end correction.
    groups = np.array([[0, dp, 0, 0, 5*dp, 0, 0, dp],
                       [0, 0, 0, 1, 5, 1, 0, 0],
                       [1, 1, 0, 0, 10, 0, 1, 1],
                       [1, 2, 0, 0, 15, 0, 1, 2]])
    c = np.array([(CTOT-monomer)/dp, monomer, 0, 0])
    return groups, c/np.sum(c), np.sum(c)


def unifac_values(t, state):
    groups, x, molecular_c = state
    lg = unifac_log_gamma(t, groups, x)
    dg = R*t*(lg[3]-lg[2]-lg[1])
    # nu=-1. Change from pure-liquid mole standard to a fixed 1 M standard.
    concentration_term = R*t*math.log(molecular_c/C0)
    return {"G_excess": float(dg), "G_mixing_1M": float(dg+concentration_term),
            "monomer_activity": float(x[1]*np.exp(lg[1])), "log_gamma": lg.tolist()}


def fh_activity(t, monomer, chi_a=0.0, chi_b=0.0, dp=DP):
    """FH pure-monomer activity under equal repeat/monomer lattice volumes.

    chi=chi_a+chi_b/T. The default is an athermal control, not measured styrene/PS chi.
    """
    phi = monomer/CTOT
    return phi*math.exp((1-1/dp)*(1-phi)+(chi_a+chi_b/t)*(1-phi)**2)


def self_check():
    # Public FeOs component tests, propane at 250 K and 1 molecule/1000 A3.
    c = np.array([1e27/NA])
    parts = helmholtz_density(250., c, [2.001829], [3.618353], [208.1101], [[0]], True)/c[0]
    expected = [0.410610492598808, -0.12402626171926148, -1.0622531100351962]
    np.testing.assert_allclose(parts, expected, atol=2e-12, rtol=0)
    for i in range(4):
        groups, _, _ = unifac_state()
        x = np.zeros(4); x[i] = 1
        assert abs(unifac_log_gamma(700, groups, x)[i]) < 1e-12
    print("I042 FeOs three component regressions and UNIFAC pure-component limits reproduced")


if __name__ == "__main__":
    self_check()
