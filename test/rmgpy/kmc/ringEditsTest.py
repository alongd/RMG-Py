"""Ring edits on declared chains (I-264).

A declared chain admits ring-touching chemistry (``ChemistryPolicy.ring_edits
== "admit"``); the legacy v1 chemistry does not.  What is admitted is decided
by the window: closure holds every ring of at most ``RING_LIMIT`` atoms
through a window atom whole, and a window that still cuts a ring (a bond
leaving it lies on a ring) is refused by name, never truncated.  At emission,
a record whose channel RMG reads under another node of the family's tree on
the record's own graphs (the node depends on the drawing of an aromatic ring)
is refused by name, not an error.

Fast part: the policies and closure on SMILES witnesses (no database).  The
novolac part compiles one short-molecule proxy of the v2 fixture for one
family and emits it (RMG database at ``harness.DATABASE``)."""

import json
from pathlib import Path

import pytest

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import rmg_moves as rm
from rmgpy.kmc import windows as win

DECKS = Path(__file__).parent / "fixtures" / "decks"
NOVOLAC = DECKS / "novolac_v2_deck.py"


# --------------------------------------------------------------------------
# the policies (no database)
# --------------------------------------------------------------------------
def test_declared_chains_admit_ring_edits_and_the_legacy_chemistry_does_not():
    assert win.DECLARED_INTERIM.ring_edits == "admit"
    assert win.declared_interim({"C", "H", "O"}).ring_edits == "admit"
    assert win.LEGACY_V1.ring_edits == "exclude"
    # the policy names (and so every certificate's domain name) are unchanged
    assert win.declared_interim({"O"}).name == "declared_interim_v1[C,H,O]"
    assert win.DECLARED_INTERIM.domain.name == "declared_interim_v1"


# --------------------------------------------------------------------------
# closure on rings (no database)
# --------------------------------------------------------------------------
def smiles_witness(smiles, centre=0):
    """The molecule of ``smiles`` (atom ids in SMILES order) as a witness
    whose product side puts a radical on ``centre``."""
    from rmgpy.molecule import Molecule
    m = Molecule(smiles=smiles)
    for k, atom in enumerate(m.atoms):
        atom.id = k
    rows = rm.molecule_rows(m)
    product = dict(rows)
    product[centre] = rows[centre][:3] + (1,) + rows[centre][4:]
    return b.make_witness(tuple(sorted(rows.items())), tuple(sorted(product.items())), ())


def heavy(window, witness):
    rows = witness.reactant_rows()
    return {a for a in window.atoms if rows[a][0] != "H"}


def test_closure_holds_every_small_ring_through_a_window_atom_whole():
    """Decalin, the support on a ring atom next to the fusion: both rings
    pass through the fusion atoms, so both are closed (one smallest ring per
    atom closes one of them and cuts the other)."""
    w = smiles_witness("C1CCC2CCCCC2C1", centre=1)
    window = win.close(w, win.Support(atoms=frozenset({1}), reacting=frozenset({1})))
    assert window.refusal is None
    assert heavy(window, w) == set(range(10))


def test_a_window_that_cuts_a_ring_is_refused_by_name():
    """A six-membered ring fused to a ten-membered one (above ``RING_LIMIT``):
    the six-ring is closed whole, its fusion atoms' smallest ring is the
    six-ring, and the ten-ring leaves the window through two ring bonds."""
    w = smiles_witness("C12CCCCC1CCCCCCCC2", centre=2)
    support = win.Support(atoms=frozenset({2}), reacting=frozenset({2}))
    window = win.close(w, support)
    assert window.refusal is not None and window.refusal.startswith("refused:ring_cut")
    assert "10-membered" in window.refusal
    assert window.atoms == support.atoms        # refused, never truncated


def test_disjoint_rings_close_as_one_smallest_ring_each():
    """Two phenyls on a CH2 (the declared domain's rings): the window of a
    ring atom is that ring and its hydrogens, as before; no ring bond leaves
    it."""
    w = smiles_witness("c1ccc(Cc2ccccc2)cc1", centre=0)
    window = win.close(w, win.Support(atoms=frozenset({0}), reacting=frozenset({0})))
    assert window.refusal is None
    assert heavy(window, w) == {0, 1, 2, 3, 11, 12}
    assert win.boundary_bond_problems(w, window) == []


# --------------------------------------------------------------------------
# novolac v2: one short-molecule proxy, one family, compiled and emitted
# --------------------------------------------------------------------------
def _database_or_skip():
    import rmg_moves_harness as harness
    if not harness.DATABASE.is_dir():
        pytest.skip(f"needs the RMG database at {harness.DATABASE}")
    return str(harness.DATABASE)


def novolac_subset(tmp_path_factory, family, proxy):
    """Compile ``family`` on the novolac v2 fixture's proxy ``proxy`` alone
    (both depths), emit every record against the roster of that proxy, and
    assemble: (rows, catalog)."""
    from rmgpy.kmc import declared_proxies as dp
    from rmgpy.kmc import emit
    from rmgpy.kmc.compile import CompileRequest, compile_artifact
    from rmgpy.kmc.run import read_deck
    database = _database_or_skip()
    out = tmp_path_factory.mktemp(family.replace(",", "_"))
    original = dp.sequence_proxy_set

    def subset(spec, units, families):
        proxies, record = original(spec, units, families)
        kept = tuple(p for p in proxies if p.site_type.split("@")[0] == proxy)
        record = dict(record, pairs={d: [] for d in record["pairs"]})
        return kept, record

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(dp, "sequence_proxy_set", subset)
        path, _id = compile_artifact(CompileRequest(
            database_path=database, output_dir=str(out), families=(family,),
            family_universe=None, database_sha=None, barrier_e0_fixed_b=None, proxy_units=1,
            deck=str(NOVOLAC)))
    artifact = json.loads(Path(path).read_text())
    roster = emit.declared_reference_systems(read_deck(str(NOVOLAC)).spec, 1,
                                             pairs=emit.declared_roster_pairs(artifact))
    emit.write_roster(roster, out / "roster.json")
    doc = json.loads((out / "roster.json").read_text())
    doc["systems"] = {n: v for n, v in doc["systems"].items()
                      if n.split(":", 1)[1].split("@")[0] == proxy}
    (out / "roster_proxy.json").write_text(json.dumps(doc, sort_keys=True,
                                                      separators=(",", ":")) + "\n")
    emit.main([str(path), str(out / "roster_proxy.json"), "--database", database,
               "--out", str(out / "rows.json")])
    emit.main(["assemble", str(out / "rows.json"), "--out", str(out / "catalog.json")])
    decoder, text, k, rows = json.JSONDecoder(), (out / "rows.json").read_text(), 0, []
    while k < len(text):
        obj, k = decoder.raw_decode(text, k)
        while k < len(text) and text[k].isspace():
            k += 1
        if obj.get("kind") == "row":
            rows.append(obj)
    return artifact, rows, json.loads((out / "catalog.json").read_text())


@pytest.fixture(scope="module")
def shift(tmp_path_factory):
    return novolac_subset(tmp_path_factory, "1,2_shiftC", "short_radical:1")


@pytest.fixture(scope="module")
def exocyclic(tmp_path_factory):
    return novolac_subset(tmp_path_factory, "Intra_R_Add_Exocyclic", "short_radical:2")


def test_novolac_ring_edits_emit_templates(shift):
    """Every 1,2_shiftC mapping of the novolac short molecule touches a ring
    atom (each was ``refused:unmapped`` under the excluding policy): the
    records emit, each from fixed-point members, and the catalog assembles
    complete."""
    artifact, rows, catalog = shift
    assert rows and len(rows) == len(artifact["records"])
    assert not [r for r in rows if r.get("refusal") == "refused:unmapped"]
    emitted = [r for r in rows if r["status"] == "emitted"]
    assert emitted
    assert all(m["status"] == "fixed_point" for r in emitted for m in r["member_rows"])
    assert catalog["complete"] and catalog["gates"]["failures"] == []
    templates = [t for p in catalog["payloads"] for t in p["templates"]]
    assert templates and {t["family"] for t in templates} == {"1,2_shiftC"}


def test_a_drawing_dependent_node_is_refused_by_name(exocyclic):
    """Intra_R_Add_Exocyclic closes a three-membered ring onto an aromatic
    ring; its reversed (ring-opening) record is read by RMG under another
    node on the record's own graphs, so no count of its channel can be made:
    refused by name, never an emission error."""
    _artifact, rows, catalog = exocyclic
    assert not [r for r in rows if r["status"] == "error"]
    named = [r for r in rows if r.get("refusal") == "refused:node_attribution"]
    assert named and all(r["orientation"] == "reversed" for r in named)
    reasons = [r["payload"]["refusal"]["reason"] for r in named]
    assert all("node" in why for why in reasons)
    assert not [f for f in catalog["gates"]["failures"] if f.startswith("error ")]
