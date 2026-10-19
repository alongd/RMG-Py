"""Records of a declared polymer certify under the declaration's own species
rules (``emit.certification_rules``) and are refused under the v1 PS rules.
Needs the RMG database (``rmg_moves_harness.DATABASE``)."""

from pathlib import Path

import pytest

import channelsTest as ct
import rmg_moves_harness as harness
from rmgpy.kmc import certificates as cm
from rmgpy.kmc import emit
from rmgpy.kmc.run import chemistry_sha256, declaration_sha256, read_deck

pytestmark = pytest.mark.skipif(not harness.DATABASE.is_dir(),
                                reason=f"needs the RMG database at {harness.DATABASE}")

DECKS = Path(__file__).parent / "fixtures" / "decks"


@pytest.fixture(scope="module")
def fam():
    return ct.database().kinetics.families


def _declared_artifact(spec):
    return {"provenance": {"proxy_source": {
        "kind": "declaration", "chemistry_sha256": chemistry_sha256(spec),
        "declaration_sha256": declaration_sha256(spec)}}, "records": []}


def _pair(fam, smiles):
    from certificatesTest import _record_pair
    owner, partner, _cs, _ch = _record_pair(
        fam, "R_Recombination", [ct.mol(smiles[0], 1), ct.mol(smiles[1], 100)])
    return [owner, partner]


@pytest.mark.parametrize("deck, smiles", [
    ("pp_deck.py", ("C[C](C)CC(C)C", "[CH3]")),           # a tertiary PP radical
    ("novolac_deck.py", ("Cc1ccccc1[O]", "[CH3]")),       # a phenoxy radical
])
def test_records_certify_under_declared_rules_and_are_refused_under_ps_rules(fam, deck, smiles):
    spec = read_deck(DECKS / deck).spec
    records = _pair(fam, smiles)
    rules = emit.certification_rules(_declared_artifact(spec), spec)
    declared, _store, _ = cm.certify_records(records, fam, cm.Certifier(fam, **rules))
    assert [r.status for r in declared] != ["refused", "refused"]
    assert all(r.status != "refused" for r in declared), [r.status_reason for r in declared]
    legacy, _store, _ = cm.certify_records(records, fam, cm.Certifier(fam))
    assert all(r.status == "refused" for r in legacy)
    assert all("refused:domain" in r.status_reason for r in legacy), \
        [r.status_reason for r in legacy]
