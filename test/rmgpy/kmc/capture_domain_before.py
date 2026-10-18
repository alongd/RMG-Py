"""Freeze ``domain.check`` on the molecules of the existing domain tests
(certificatesTest: accepted references, refusals, the multiplicity
convention, the styryl end) -> ``data/domain_check_before.json``.  Run ONCE
on the S2a source; database-free.  ``declaredPolicyTest`` compares the live
``check`` against it (``universal_check`` must leave ``check`` unchanged)."""

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def domain_molecules() -> dict:
    """name -> Molecule (atom ids from 0, RMG atom order), the inputs of
    certificatesTest's domain tests."""
    import channelsTest as ct
    mol = ct.mol
    radical, closed, h_atom = ct.dp3_molecules()
    out = {
        "DP3 radical": radical, "DP3 closed": closed, "H": h_atom,
        "DP8 head radical": mol("[CH2]C(c1ccccc1)" + ct.ps(7), 0),
        "DP4": mol(ct.ps(4), 0), "1-phenylethyl": mol("C[CH]c1ccccc1", 0),
        "DP8 closed": mol(ct.ps(8), 0), "CH3": mol("[CH3]", 0), "C2H6": mol("CC", 0),
        "CH4": mol("C", 0), "H2": mol("[H][H]", 0),
        "vinylidene end": mol("C=C(c1ccccc1)CC(c1ccccc1)C", 0),
        "styryl end": mol("c1ccccc1C=CCC(c1ccccc1)C", 0),
        "Kekule DP2": mol("CC(C1=CC=CC=C1)CC(C1=CC=CC=C1)", 0),
        "two far benzylic radicals": mol("C[C](c1ccccc1)" + ct.ps(2) + "C[C](c1ccccc1)C", 0),
        "secondary backbone radical": mol("CC(c1ccccc1)[CH]C(c1ccccc1)C", 0),
        "styryl with a second phenyl": mol("c1ccccc1C=C(c1ccccc1)CC", 0),
        "styrene": mol("C=Cc1ccccc1", 0),
    }
    refused = {
        "branched": "CC(C)(c1ccccc1)CC",
        "aryl radical on a chain": "CC(c1ccc[c]c1)CC(c1ccccc1)C",
        "phenyl radical": "[c]1ccccc1",
        "oxygen": "CC(O)C",
        "internal unsaturation": "CC(c1ccccc1)=CC(c1ccccc1)C",
        "radicals two bonds apart": "C[C](c1ccccc1)C[CH]C",
        "three radicals": "[CH2]C(c1ccccc1)" + ct.ps(3) + "C[C](c1ccccc1)" + ct.ps(3) +
                          "C[CH]c1ccccc1",
        "allylic radical": "C=C(c1ccccc1)[CH]C(c1ccccc1)C",
        "naphthyl": "CC(c1ccc2ccccc2c1)C",
        "carbene": "[CH2]",
        "two phenyls on one carbon": "CC(c1ccccc1)(c1ccccc1)C",
    }
    out.update({f"refused: {k}": mol(v, 0) for k, v in refused.items()})
    singlet = mol("C[C](c1ccccc1)" + ct.ps(2) + "C[C](c1ccccc1)C", 0)
    singlet.multiplicity = 1
    out["singlet biradical"] = singlet
    return out


def verdicts() -> dict:
    from rmgpy.kmc import domain as dom
    return {name: json.loads(json.dumps(dom.check(m).as_dict(), sort_keys=True, default=str))
            for name, m in domain_molecules().items()}


if __name__ == "__main__":
    path = HERE / "data" / "domain_check_before.json"
    path.write_text(json.dumps(verdicts(), sort_keys=True, indent=1) + "\n")
    print(path)
