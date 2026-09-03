#!/usr/bin/env python3
"""Emit a polymer_pools.json sidecar for a pool carrying k_unzip > 0.

Runs under rmg_env against THIS worktree. Writes the artifact to argv[1] and
asserts locally that the unzip channel still carries A > 0; the liveness
verdict itself is CKMG's and is taken in a separate ck_env process
(i055_sidecar_audit.py), because the two live in different conda envs.
"""

import json
import sys

from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.polymer import Polymer, build_polymer_moments_artifact

K_UNZIP = 100.0


def _spc(smiles, label, index=-1):
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.label = label
    s.index = index
    return s


def _mu_dummy(label):
    s = Species(label=label, reactive=False)
    s.molecule = [Molecule().from_smiles("[Ne]")]
    s.is_moment_dummy = True
    s.index = -1
    return s


def main():
    out_path = sys.argv[1]
    pool = Polymer(
        label="PE", monomer="[CH2][CH2]", end_groups=["[H]", "[H]"],
        cutoff=3, Mn=1500.0, Mw=1800.0, initial_mass=1.0,
        k_scission=0.0, k_unzip=K_UNZIP,
    )
    core = [
        _spc("CC", "PE", index=2),
        _mu_dummy("PE_mu0"), _mu_dummy("PE_mu1"), _mu_dummy("PE_mu2"),
        _spc("[CH3]", "G", index=7),
    ]
    core[0].is_polymer_proxy = True

    artifact = build_polymer_moments_artifact(
        [pool], core_species=core, core_reactions=[],
        configured_pool_labels=["PE"], condensed_species=core[:4],
    )
    with open(out_path, "w") as f:
        json.dump(artifact, f, indent=1)

    block = artifact["pools"][0]
    a = float(block["channels"]["unzip"]["A"])
    print(f"  emitted {out_path}")
    print(f"  schema_version      = {artifact['schema_version']}")
    print(f"  pool                = {block['label']!r}")
    print(f"  channels.unzip.A    = {a!r}   (configured k_unzip = {K_UNZIP!r})")
    if a != K_UNZIP:
        print(f"  FAIL: unzip A {a!r} != configured k_unzip {K_UNZIP!r}")
        return 1
    if not a > 0.0:
        print("  FAIL: unzip channel A is not > 0 -- the channel would audit dead")
        return 1
    print("  OK: unzip channel emitted with A > 0")
    return 0


if __name__ == "__main__":
    sys.exit(main())
