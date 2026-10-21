"""Window extension of a linear statistical copolymer (schema v2): past an
artificial port the next residue is not fixed by the declaration, so the
recertified fixed point extends a sequence chain over EVERY residue the
sequence support allows there (one branch per residue word) and accepts a
window only when every branch gives the same window and certificate
(selection, support, outcome).  Branches that disagree are refused by name
(``refused:extension_branches``) with the breaking branch; no residue stands
for the others.

The v1 homopolymer fixed points (a PE and a PP declared chain) are pinned
byte-identical.

Needs the RMG database (``rmg_moves_harness.DATABASE``)."""

import hashlib
import json

import pytest

import channelsTest as ct
import declared_parity_support as sup
import rmg_moves_harness as harness
from polymer_spec_v2_support import DECKS
from rmgpy.kmc import boundaries as b
from rmgpy.kmc import certificates as cm
from rmgpy.kmc import declared_chain as dc
from rmgpy.kmc import template_emission as te
from rmgpy.kmc import windows as win
from rmgpy.kmc.declared_chain import EditOp
from rmgpy.kmc.declared_sequence import SequenceChainSpec, sequence_provider
from rmgpy.kmc.run import read_deck

pytestmark = pytest.mark.skipif(not harness.DATABASE.is_dir(),
                                reason=f"needs the RMG database at {harness.DATABASE}")

# sha256 of json.dumps(value, sort_keys=True, default=repr) at 2d76fceb5 (before
# the branching): /home/alon/runs/pk-i236/probe/p3_v1_pins.py
V1_WINDOW_PINS = {
    "pe_fixed_points": "ddd03c84ed851d78684a7dfa084df85600fea7432bb5d3c64e3c9748a24a339b",
    "pe_member_row": "938a109f74287efd591ab1482a09bad3541d919b9b5d0bfe22546e0c15ad278c",
    "pp_fixed_points": "8999e663286abd624103be7c67b74febc1e01b4ff80269249b173aaddb32889a",
    "pp_member_row": "d76e95b34e0e5e9c14f0df6611c37079db72b8551976dff9d037d4e1868ca888",
}


@pytest.fixture(scope="module")
def fam():
    return ct.database().kinetics.families


@pytest.fixture(scope="module")
def certifier(fam):
    return cm.Certifier(fam)


@pytest.fixture(scope="module")
def ep():
    return sequence_provider(read_deck(DECKS / "ep_v2_deck.py").spec)


def _radical_chain(residues, unit, key):
    return SequenceChainSpec(residues, edits=(EditOp("remove_h", unit, key),
                                              EditOp("set_radical", unit, key, 1)))


def _member(certifier, family, provider, chain, suffix):
    specs = (win.DeclaredSpec(provider.name, chain),)
    cases = [c for c in win._cover_cases(certifier, family, specs, "forward")
             if c.label.endswith(suffix)]
    assert cases, suffix
    return cases[0]


def _fixed_point(certifier, family, case):
    observer = win.MemberObserver(certifier, family, case.direction, case.wanted)
    return win.fixed_point(case.witness, observer)


def test_a_window_away_from_the_ports_holds_on_every_residue_beyond_them(fam, certifier, ep):
    """An H shift in the middle of EPEPE (the radical on the first P's
    methine): the window stays inside the chain, and every check (each port
    alone and both jointly, at +1 and +2 residues) holds on every residue word
    the support allows there, 2 ** (residues added) branches each.  The fixed
    point cuts a local template."""
    family = fam["intra_H_migration"]
    case = _member(certifier, family, ep, _radical_chain(("E", "P", "E", "P", "E"), 1, "p2"),
                   "R2H_S;C_rad_out_Cs2;Cs_H_out_2H|[CH2]C(CCC)CCCC(C)CC|m0")
    fp = _fixed_point(certifier, family, case)
    assert fp.status == "fixed_point", (fp.status, fp.reason)
    assert all(depth == 0 for _port, depth in fp.depths)
    checks = fp.checks[-1]
    assert {(len(c["ports"]), c["step"]) for c in checks} == {(1, 1), (1, 2), (2, 1), (2, 2)}
    for check in checks:
        assert check["branches"] == 2 ** (len(check["ports"]) * check["step"]), check
        assert check["same_window"] and check["same_certificate"], check
    cut = te.cut_template(fp.witness, fp.window.atoms, family=family.label,
                          pathway=case.pathway, direction=case.direction,
                          labels=te.labels_of(case.moves[0]))
    assert not cut.whole_species and cut.ports


def test_branches_that_read_different_residues_are_refused_by_name(fam, certifier, ep):
    """An H shift near the head of EEEEE reads past the head port: extended by
    E and by P the windows differ, so the window depends on the residue the
    declaration leaves open.  Refused by name with the breaking branch; never
    the E window standing for both."""
    family = fam["intra_H_migration"]
    case = _member(certifier, family, ep, _radical_chain(("E",) * 5, 2, "e1"),
                   "|CCC[CH]CCCCCC|m0")
    fp = _fixed_point(certifier, family, case)
    assert fp.status == "refused:extension_branches", (fp.status, fp.reason)
    assert "0:head+E" in fp.reason and "0:head+P" in fp.reason, fp.reason
    assert fp.window is None


def test_an_extension_step_adds_one_named_residue_and_verifies(fam, certifier, ep):
    """One step at a sequence port adds the named residue (an allowed one)
    by its declared junction: the step records it, the witness carries the
    word per port, and the map verifies (replayed from the original witness).
    Without a residue a copolymer port still refuses: the declaration does
    not fix it."""
    case = _member(certifier, fam["intra_H_migration"], ep,
                   _radical_chain(("E", "P", "E", "P", "E"), 1, "p2"),
                   "R2H_S;C_rad_out_Cs2;Cs_H_out_2H|[CH2]C(CCC)CCCC(C)CC|m0")
    witness, tail = case.witness, (0, "tail")
    assert b.extension_residues(witness, tail) == ("E", "P")
    grown, record = b.extend(witness, tail, residue="P")
    assert [step.residue for step in record.steps] == ["P"]
    assert dict(grown.extension_words) == {tail: ("P",)}
    longer, _record = b.extend(grown, tail, residue="E")
    assert dict(longer.extension_words) == {tail: ("P", "E")}
    added = longer.heavy_count() - witness.heavy_count()
    assert added == 5        # P (3 carbons) then E (2)
    assert grown != b.extend(witness, tail, residue="E")[0]
    with pytest.raises(b.ExtensionRefused, match="refused:sequence_extension"):
        b.extend(witness, tail)
    with pytest.raises(b.ExtensionRefused, match="not a declared residue"):
        b.extend(witness, tail, residue="X")


def test_v1_homopolymer_fixed_points_are_byte_identical(fam, certifier):
    """PE (an H shift reading one repeat past the tail) and PP (one member at
    depths 3 and 2) fixed points and full member rows (contexts, covers) as at
    2d76fceb5: one repeat per step, no branch."""
    family = fam["intra_H_migration"]
    pe = dc.DeclaredChainProvider(sup.pe_declaration())
    pp = dc.DeclaredChainProvider(read_deck(DECKS / "pp_deck.py").spec)
    got = {}
    for name, provider, chain in (
            ("pe", pe, dc.ChainSpec(9, edits=(EditOp("remove_h", 4, "c1"),
                                              EditOp("set_radical", 4, "c1", 1)))),
            ("pp", pp, dc.ChainSpec(4, edits=(EditOp("remove_h", 1, "c2"),
                                              EditOp("set_radical", 1, "c2", 1))))):
        specs = (win.DeclaredSpec(provider.name, chain),)
        cases = win._cover_cases(certifier, family, specs, "forward")
        fps = []
        for case in cases[:4]:
            fp = _fixed_point(certifier, family, case)
            fps.append([case.label, fp.status, fp.reason, repr(fp.depths), fp.retries,
                        sorted(fp.window.atoms) if fp.reached else None,
                        repr(fp.observation.view(fp.window.atoms)) if fp.reached else None,
                        repr(fp.checks)])
        row = win.run_member(certifier, family, cases[0], specs=specs)
        row.pop("seconds", None)
        got[f"{name}_fixed_points"] = _digest(fps)
        got[f"{name}_member_row"] = _digest(row)
    assert got == V1_WINDOW_PINS


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=repr).encode()).hexdigest()
