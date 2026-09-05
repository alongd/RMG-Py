"""
I-103 measurement probe (measurement only; changes nothing).

Question: is the whole-molecule group-additivity thermo that Polymer bare-assigns
(`self.thermo = proxy.thermo`) per-chain (extensive, grows with DP) or per-site
(intensive, flat in DP)?

Method: build the pool's OWN capped chain at DP in {2,3,4,6} via the class's own
`Polymer._capped_chain_species(dp)` -- the same stitching recipe the proxy uses --
so the ONLY thing that varies across the four cases is the repeat-unit count. Then
obtain thermo SYNCHRONOUSLY via `ThermoDatabase.get_thermo_data(species)` (NOT the
async thermoengine.submit route that stalled the previous attempt).
"""
import os

from rmgpy import settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.polymer import Polymer


def load_db():
    db = RMGDatabase()
    db.load_thermo(
        os.path.join(settings["database.directory"], "thermo"),
        thermo_libraries=["primaryThermoLibrary"],
        surface=False,
    )
    return db


def main():
    db = load_db()
    tdb = db.thermo

    # Backend MRO (Verifier item 1): the class whose get_thermo_data we call.
    print("=== Backend MRO (ThermoDatabase.get_thermo_data owner) ===")
    for c in type(tdb).__mro__:
        print("   ", c.__module__ + "." + c.__name__)

    # PE pool: ethylene diradical monomer, methyl/H end caps. A plain aliphatic
    # backbone keeps group additivity unambiguous and avoids aromatic resonance
    # confounds. The monomer and both end caps are FIXED across all DP.
    # Mn/Mw/initial_mass are pool-distribution moments required by the
    # constructor; they do NOT enter _capped_chain_species(dp), which builds
    # from end_groups + monomer + dp only. Fixed identically for all DP.
    pol = Polymer(label="PE", monomer="[CH2][CH2]", end_groups=["[CH3]", "[H]"],
                  Mn=5000.0, Mw=6000.0, initial_mass=1.0)

    rows = []
    for dp in (2, 3, 4, 6):
        spc = pol._capped_chain_species(dp)
        if spc is None:
            raise RuntimeError(f"DP={dp}: capped chain not constructible")
        mol = spc.molecule[0]
        formula = mol.get_formula()
        n_heavy = sum(1 for a in mol.atoms if a.is_non_hydrogen())

        thermo = tdb.get_thermo_data(spc)
        if thermo is None:
            raise RuntimeError(f"DP={dp}: get_thermo_data returned None")

        h298 = thermo.get_enthalpy(298.0) / 1000.0   # J/mol -> kJ/mol
        s298 = thermo.get_entropy(298.0)             # J/mol/K
        cp1000 = thermo.get_heat_capacity(1000.0)    # J/mol/K
        rows.append((dp, formula, n_heavy, h298, s298, cp1000))

        print(f"\n=== DP={dp}  formula={formula}  heavy_atoms={n_heavy} ===")
        print(f"    thermo class : {type(thermo).__name__}")
        print(f"    comment      : {thermo.comment[:200]}")
        print(f"    H298         : {h298:10.3f} kJ/mol")
        print(f"    S298         : {s298:10.3f} J/mol/K")
        print(f"    Cp(1000K)    : {cp1000:10.3f} J/mol/K")

    print("\n=== TABLE ===")
    print(f"{'DP':>3} {'formula':>10} {'heavy':>6} {'H298[kJ/mol]':>14} {'S298[J/mol/K]':>15} {'Cp1000[J/mol/K]':>16}")
    for dp, formula, nh, h, s, cp in rows:
        print(f"{dp:>3} {formula:>10} {nh:>6} {h:14.3f} {s:15.3f} {cp:16.3f}")

    # Per-DP increments (a first difference; equal increments => linear/extensive)
    print("\n=== First differences (per added repeat unit) ===")
    print(f"{'DP->DP':>8} {'dH298':>12} {'dS298':>12} {'dCp1000':>12} {'dHeavy':>8}")
    for i in range(1, len(rows)):
        dp0, _, nh0, h0, s0, cp0 = rows[i - 1]
        dp1, _, nh1, h1, s1, cp1 = rows[i]
        span = dp1 - dp0
        print(f"{dp0:>3}->{dp1:<3} "
              f"{(h1 - h0) / span:12.3f} {(s1 - s0) / span:12.3f} "
              f"{(cp1 - cp0) / span:12.3f} {(nh1 - nh0) / span:8.2f}")


if __name__ == "__main__":
    main()
