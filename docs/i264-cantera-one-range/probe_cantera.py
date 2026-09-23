import cantera as ct
import numpy as np

print("=== cantera", ct.__version__, "===")

# Does the Python object API have NasaPoly1 (single-range)?
print("has NasaPoly1 attr:", hasattr(ct, "NasaPoly1"))
print("has NasaPoly2 attr:", hasattr(ct, "NasaPoly2"))

# A real single-range NASA7 coefficient set for monatomic Ar (Cp = 5/2 R constant).
# a0 = 2.5, a1..a4 = 0, a5 = enthalpy const, a6 = entropy const.
coeffs = [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.37967]

def try_yaml(label, tranges, data):
    y = f"""
phases:
- name: gas
  thermo: ideal-gas
  elements: [Ar]
  species: [Ar]
  state: {{T: 300.0, P: 1 atm}}
species:
- name: Ar
  composition: {{Ar: 1}}
  thermo:
    model: NASA7
    temperature-ranges: {tranges}
    data: {data}
"""
    try:
        sol = ct.Solution(yaml=y)
        sp = sol.species('Ar')
        # report thermo back
        cp = []
        for T in (300.0, 1000.0, 3000.0):
            sol.TP = T, ct.one_atm
            cp.append((T, sol.cp_mole/ct.gas_constant, sol.enthalpy_mole, sol.entropy_mole))
        print(f"[OK]   {label}: loaded. thermo type = {type(sp.thermo).__name__}")
        for T, cpr, h, s in cp:
            print(f"          T={T:6.0f}  Cp/R={cpr:.4f}  H={h:.3f}  S={s:.4f}")
    except Exception as e:
        print(f"[FAIL] {label}: {type(e).__name__}: {str(e).splitlines()[0]}")

# 1) genuine one-range: 2 breakpoints, 1 coeff set
try_yaml("one-range (2 breakpoints, 1 data)", [200.0, 6000.0], [coeffs])
# 2) duplicated across two ranges: 3 breakpoints, 2 identical coeff sets
try_yaml("two-range duplicated (3 breakpoints, 2 data)", [200.0, 1000.0, 6000.0], [coeffs, coeffs])
