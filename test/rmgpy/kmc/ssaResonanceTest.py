"""The SSA multiplier's same-species test is resonance-aware and spin-aware.

A compiled bimolecular record carries ``ssa_multiplier`` 2 (and the
``unordered-identical-pair`` convention) exactly when its two reactants are
one species: one resonance set at one spin multiplicity.  RMG's generated
reactions hold one resonance form per reactant Species, so two forms of one
species (crotyl head and tail) must still compare equal, and two spin states
of one skeleton must not.

Light: the unit cases need no database; the production cases load the
R_Recombination family only (``RMG_DATABASE_PATH`` or
/home/alon/runs/phase2b/database).
"""

import os
from pathlib import Path
from types import SimpleNamespace

import pytest

DATABASE = Path(os.environ.get("RMG_DATABASE_PATH", "/home/alon/runs/phase2b/database"))
IDENTICAL = "unordered-identical-pair N(N-1)/2"
DISTINCT = "distinct-pair N_A*N_B"
CROTYL = ("C=C[CH]C", "[CH2]C=CC")


def _species(smiles, multiplicity=None):
    from rmgpy.molecule import Molecule
    from rmgpy.species import Species
    mol = Molecule().from_smiles(smiles)
    if multiplicity is not None:
        mol.multiplicity = multiplicity
    return Species(molecule=[mol])


def _verdicts(a, b):
    from rmgpy.kmc.compiler import _pair_convention
    from rmgpy.kmc.event_record import ssa_multiplier_for
    return ssa_multiplier_for([a, b]), _pair_convention([a, b])


# --------------------------------------------------------------------------
# the same-species test
# --------------------------------------------------------------------------
@pytest.mark.parametrize("first,second", [CROTYL, CROTYL[::-1]])
def test_crotyl_forms_are_one_species(first, second):
    assert _verdicts(_species(first), _species(second)) == (2.0, IDENTICAL)


@pytest.mark.parametrize("first,second", [
    ("[CH2]c1ccccc1", "C=C1C=CC=C[CH]1"),
    ("C=C1C=CC=C[CH]1", "[CH2]c1ccccc1"),
    ("[CH2]C1=CC=CC=C1", "C=C1C=C[CH]C=C1"),
])
def test_benzyl_forms_are_one_species(first, second):
    """Benzyl drawn exocyclic (aromatic or Kekule) and quinoid (ortho, para)."""
    assert _verdicts(_species(first), _species(second)) == (2.0, IDENTICAL)


@pytest.mark.parametrize("first,second", [("C=C[CH]C", "[CH2]C=C"), ("[CH2]C=C", "[CH2]C=CC")])
def test_crotyl_and_allyl_are_two_species(first, second):
    assert _verdicts(_species(first), _species(second)) == (1.0, DISTINCT)


@pytest.mark.parametrize("triplet,singlet", [
    ("[CH2]CCCC[CH2]", "[CH2]CCCC[CH2]"),     # one biradical skeleton, two spins
    ("[CH2]C=CC[CH2]", "C=C[CH]C[CH2]"),      # resonance-related forms, two spins
])
@pytest.mark.parametrize("order", ["triplet first", "singlet first"])
def test_two_spin_states_are_two_species(triplet, singlet, order):
    """Resonance generation resets a singlet biradical to the triplet, so the
    multiplicity is compared on the original forms."""
    a, b = _species(triplet), _species(singlet, multiplicity=1)
    assert a.molecule[0].multiplicity == 3 and b.molecule[0].multiplicity == 1
    if order == "singlet first":
        a, b = b, a
    assert _verdicts(a, b) == (1.0, DISTINCT)


def test_resonance_related_biradical_forms_at_one_spin_are_one_species():
    """The control for the spin case: the same two forms, both triplet."""
    assert _verdicts(_species("[CH2]C=CC[CH2]"), _species("C=C[CH]C[CH2]")) == \
        (2.0, IDENTICAL)


def test_the_same_species_test_leaves_the_species_alone():
    """Resonance structures are generated on copies: a record's fields are
    captured from the reaction's Species after the test."""
    a, b = _species(CROTYL[0]), _species(CROTYL[1])
    before = [(len(s.molecule), s.molecule[0].to_adjacency_list()) for s in (a, b)]
    _verdicts(a, b)
    assert [(len(s.molecule), s.molecule[0].to_adjacency_list()) for s in (a, b)] == before


def test_bare_molecules_and_one_object():
    from rmgpy.kmc.event_record import same_species
    from rmgpy.molecule import Molecule
    a, b = (Molecule().from_smiles(s) for s in CROTYL)
    assert same_species(a, b) and same_species(b, a)
    assert same_species(a, a)
    assert not same_species(a, Molecule().from_smiles("[CH2]C=C"))


@pytest.mark.parametrize("first,second", [
    ("[CH2]CC", "C[CH]C"),            # n-propyl / isopropyl
    ("[CH2]CC=C", "[CH2]C=CC"),       # but-3-en-1-yl / crotyl (head form)
    ("[CH2]CC=C", "C=C[CH]C"),        # but-3-en-1-yl / crotyl (tail form)
    ("[CH]=CCC", "[CH2]C=CC"),        # but-1-en-1-yl / crotyl
])
@pytest.mark.parametrize("order", ["as listed", "swapped"])
def test_constitutional_isomers_are_two_species(first, second, order):
    """Same formula, same spin, not resonance-related: never one species."""
    a, b = _species(first), _species(second)
    if order == "swapped":
        a, b = b, a
    assert a.molecule[0].get_formula() == b.molecule[0].get_formula()
    assert a.molecule[0].multiplicity == b.molecule[0].multiplicity
    assert _verdicts(a, b) == (1.0, DISTINCT)


def test_an_unexpandable_form_is_compared_as_given(monkeypatch):
    """When RMG's resonance generation raises on one participant, that
    participant's given form still enters the comparison: the crotyl head
    form, expanded, meets the tail form as given."""
    from rmgpy.molecule import resonance
    real = resonance.generate_resonance_structures
    calls = []

    def flaky(mol, *args, **kwargs):
        calls.append(mol.to_smiles())
        if mol.is_isomorphic(_species(CROTYL[0]).molecule[0]):
            raise RuntimeError("unexpandable form")
        return real(mol, *args, **kwargs)

    monkeypatch.setattr(resonance, "generate_resonance_structures", flaky)
    assert _verdicts(_species(CROTYL[1]), _species(CROTYL[0])) == (2.0, IDENTICAL)
    assert calls, "the helper must expand forms through rmgpy.molecule.resonance"


def test_limit_charged_resonance_partners_are_not_merged():
    """Stated limit: RMG returns a charged structure without expanding its
    resonance forms (``resonance.generate_resonance_structures``), so two
    forms of one charged radical are compared as given and stay two species
    (ssa 1).  This pins today's behaviour; it is not the intended one."""
    assert _verdicts(_species("[CH2]C=C[NH3+]"), _species("C=C[CH][NH3+]")) == \
        (1.0, DISTINCT)


def test_template_rows_ssa_agrees_on_crotyl():
    """``template_emission.ssa_of_rows`` keys species through
    ``channels.SpeciesKeys`` (resonance-aware; rows carry no spin): the two
    crotyl forms are one species there too."""
    from rmgpy.kmc import rmg_moves as rm
    from rmgpy.kmc import template_emission as te
    from rmgpy.molecule import Molecule

    def rows(*smiles):
        out, start = {}, 0
        for s in smiles:
            mol = Molecule().from_smiles(s)
            for i, atom in enumerate(mol.atoms):
                atom.id = start + i
            start += len(mol.atoms)
            out.update(rm.molecule_rows(mol))
        return out

    counter = te.RmgCounter({})
    assert te.ssa_of_rows(rows(*CROTYL), counter) == 2.0
    assert te.ssa_of_rows(rows(CROTYL[0], "[CH2]C=C"), counter) == 1.0


# --------------------------------------------------------------------------
# the production path: RMG's crotyl + crotyl reactions through the
# compiler's record builder (EventSetCompiler._record, which sets
# ssa_multiplier and reactant_pair_convention)
# --------------------------------------------------------------------------
@pytest.fixture(scope="module")
def database():
    if not (DATABASE / "input").is_dir():
        pytest.skip(f"needs the RMG database at {DATABASE}")
    from rmgpy.data.rmg import RMGDatabase
    db = RMGDatabase()
    db.load_kinetics(str(DATABASE / "input" / "kinetics"), reaction_libraries=[],
                     seed_mechanisms=None, kinetics_families=["R_Recombination"],
                     kinetics_depositories=[])
    return db


def _record(reaction, orientation="as_generated"):
    """``EventSetCompiler._record`` with the rate given (as the linked-pair
    path gives it): the record fields under test need no kinetics, thermo or
    compile.  The rate is RMG's degeneracy times a unit rule (k_rule = 1)."""
    from rmgpy.kmc import compiler as cp
    stub = cp.EventSetCompiler.__new__(cp.EventSetCompiler)
    stub.reference_thermo_provider = SimpleNamespace(provenance={"test": "stub"})
    grid = [300.0, 1000.0]
    proxy = cp.SiteProxy(site_type="end", reactants=list(reaction.reactants))
    return stub._record(proxy, reaction, {"test": True}, orientation,
                        rate_override=({"T": grid, "k": [float(reaction.degeneracy)] * 2,
                                        "interpolation": "linear-ln-k",
                                        "extrapolation": "refuse"},
                                       {"kind": "RMG family estimate", "available": True,
                                        "units": "m^3/(mol*s)"}))


def _smiles(species):
    return species.molecule[0].to_smiles()


@pytest.fixture(scope="module")
def crotyl_reactions(database):
    from rmgpy.species import Species
    reactants = [Species(smiles=s) for s in CROTYL]
    for s in reactants:
        s.generate_resonance_structures(keep_isomorphic=False)
    return database.kinetics.generate_reactions_from_families(
        reactants, only_families=["R_Recombination"], resonance=True)


def test_production_head_to_tail_record_is_an_identical_pair(crotyl_reactions):
    """RMG's head-to-tail crotyl coupling (degeneracy 1: the A + A half is
    already in it) between two one-form Species is compiled with ssa 2."""
    by_product = {_smiles(r.products[0]): r for r in crotyl_reactions}
    r = by_product["C=CC(C)CC=CC"]
    assert sorted(_smiles(x) for x in r.reactants) == sorted(CROTYL)
    assert all(len(x.molecule) == 1 for x in r.reactants)
    assert r.degeneracy == 1.0
    rec = _record(r)
    assert rec.ssa_multiplier == 2.0
    assert rec.reactant_pair_convention == IDENTICAL
    # the reaction's Species are captured unchanged
    assert all(len(x.molecule) == 1 for x in r.reactants)


def test_production_same_form_records_stay_identical_pairs(crotyl_reactions):
    for r in crotyl_reactions:
        rec = _record(r)
        assert rec.ssa_multiplier == 2.0, _smiles(r.products[0])
        assert rec.reactant_pair_convention == IDENTICAL


def test_production_reverse_record_of_a_homolysis_into_two_forms(crotyl_reactions):
    """The reverse direction and the product side: the homolysis view of the
    head-to-tail coupling (``compiler._reverse_view``; its products are the
    two one-form crotyl Species) is unimolecular, and the reverse record the
    linked-pair path builds from it (``_reverse_view`` again, orientation
    ``reversed``) is an identical pair."""
    from rmgpy.kmc import compiler as cp
    by_product = {_smiles(r.products[0]): r for r in crotyl_reactions}
    homolysis = cp._reverse_view(by_product["C=CC(C)CC=CC"])
    assert sorted(_smiles(p) for p in homolysis.products) == sorted(CROTYL)
    assert _record(homolysis).ssa_multiplier == 1.0
    rec = _record(cp._reverse_view(homolysis), "reversed")
    assert (rec.ssa_multiplier, rec.reactant_pair_convention) == (2.0, IDENTICAL)


def test_production_reactions_generated_from_the_dimer(database):
    """RMG's reactions generated from the dimer come back in the family's
    forward direction with every resonance form on each radical; their
    crotyl + crotyl record is an identical pair as well."""
    from rmgpy.species import Species
    rxns = database.kinetics.generate_reactions_from_families(
        [Species(smiles="C=CC(C)CC=CC")], only_families=["R_Recombination"], resonance=True)
    pairs = [r for r in rxns if len(r.reactants) == 2
             and all({m.to_smiles() for m in x.molecule} == set(CROTYL) for x in r.reactants)]
    assert len(pairs) == 1, [[_smiles(x) for x in r.reactants] for r in rxns]
    rec = _record(pairs[0])
    assert (rec.ssa_multiplier, rec.reactant_pair_convention) == (2.0, IDENTICAL)


def test_production_crotyl_and_allyl_record_is_a_distinct_pair(database):
    from rmgpy.species import Species
    reactants = [Species(smiles="C=C[CH]C"), Species(smiles="[CH2]C=C")]
    rxns = database.kinetics.generate_reactions_from_families(
        reactants, only_families=["R_Recombination"], resonance=True)
    pairs = [r for r in rxns if len(r.reactants) == 2]
    assert pairs
    for r in pairs:
        rec = _record(r)
        assert (rec.ssa_multiplier, rec.reactant_pair_convention) == (1.0, DISTINCT), \
            [_smiles(x) for x in r.reactants]


# --------------------------------------------------------------------------
# the engines that read ssa_multiplier / reactant_pair_convention, on one
# hand-built state of four crotyl molecules (two in each localized form) and
# one of four allyl molecules.  Units: k_rule = 1, N_A * V = 1, so RMG's rate
# law is sum over channels of d * N (N - 1).
# --------------------------------------------------------------------------
def _graph_strand(strand_id, graph, names):
    from rmgpy.kmc.state import AtomRef, Strand
    formula = {}
    for node in graph.values():
        formula[node["element"]] = formula.get(node["element"], 0) + 1
    return Strand(
        strand_id, 1, formula=formula,
        radical_count=sum(node["radical"] for node in graph.values()),
        atom_refs={i: AtomRef(names[i], strand_id, 0, "feature") for i in graph},
        atom_graph={names[i]: {**{k: v for k, v in node.items() if k not in ("edges", "label")},
                               "edges": {names[j]: o for j, o in node["edges"].items()}}
                    for i, node in graph.items()})


def _legacy_state(graphs, seed):
    """One strand per molecule; atom UUIDs drawn from ``seed`` (a renaming)."""
    import random
    from rmgpy.kmc.state import KMCState
    rng = random.Random(seed)
    return KMCState([_graph_strand(f"m{m}", g, {i: f"{rng.random():.12f}" for i in g})
                     for m, g in enumerate(graphs)])


def _legacy(records, graphs, seed):
    """Per record (eligible candidates, propensity) of legacy ``ssa.py``."""
    from rmgpy.kmc import ssa
    state = _legacy_state(graphs, seed)
    index = ssa.SiteIndex(state, records)
    return [(len(ssa.eligible_candidates(state, index, r)),
             ssa.record_propensity(state, index, r, 300.0, 1.0 / ssa.N_A)) for r in records]


@pytest.fixture(scope="module")
def crotyl_records(crotyl_reactions):
    from rmgpy.kmc.state import _parse_adjacency
    out = {}
    for r in crotyl_reactions:
        rec = _record(r)
        out[_smiles(r.products[0])] = rec
    forms = {_smiles(x): _parse_adjacency(g)
             for r in crotyl_reactions for x, g in zip(r.reactants, out[_smiles(r.products[0])]
                                                       .reactant_graphs)}
    return out, forms


@pytest.fixture(scope="module")
def allyl_records(database):
    from rmgpy.kmc.state import _parse_adjacency
    from rmgpy.species import Species
    reactants = [Species(smiles="[CH2]C=C"), Species(smiles="[CH2]C=C")]
    rxns = database.kinetics.generate_reactions_from_families(
        reactants, only_families=["R_Recombination"], resonance=True)
    records = [_record(r) for r in rxns if len(r.reactants) == 2]
    assert len(records) == 1
    return records, _parse_adjacency(records[0].reactant_graphs[0])


SEEDS = range(8)
HEAD_TAIL = "C=CC(C)CC=CC"


@pytest.mark.parametrize("seed", SEEDS)
def test_legacy_head_to_tail_count_is_invariant_under_uuid_renaming(crotyl_records, seed):
    """The identical-pair dedup drops a candidate only when its reverse is a
    candidate too.  Head-to-tail between two localized forms has no reverse
    candidate: all 2 x 2 pairs stay, whatever the atom UUIDs."""
    records, forms = crotyl_records
    graphs = [forms[CROTYL[1]]] * 2 + [forms[CROTYL[0]]] * 2
    (n, a), = _legacy([records[HEAD_TAIL]], graphs, seed)
    assert (n, a) == (4, pytest.approx(8.0, rel=1e-12))


@pytest.mark.parametrize("seed", SEEDS)
def test_legacy_same_form_pair_equals_the_rate_law(allyl_records, seed):
    """Allyl + allyl (one form): N = 4 gives 4*3/2 = 6 unordered pairs, ssa 2,
    k = d = 2: propensity 2 * 2 * 6 = 24 = d * N (N - 1), as before the fix."""
    records, graph = allyl_records
    assert records[0].degeneracy == 2.0 and records[0].ssa_multiplier == 2.0
    assert _legacy(records, [graph] * 4, seed) == [(6, pytest.approx(24.0, rel=1e-12))]


def _melt_total(pairs, state_smiles):
    """The melt window matcher: each channel's template cut from its own
    two-molecule witness (whole species), weight ssa * k / T(witness), times
    its transitions on the state."""
    from rmgpy.kmc import rmg_moves as rm
    from rmgpy.kmc import template_emission as te
    from rmgpy.kmc import window_matcher as wm
    from rmgpy.molecule import Molecule

    def rows_of(smiles):
        rows, radicals, base = {}, [], 0
        for s in smiles:
            mol = Molecule().from_smiles(s)
            heavy = [x for x in mol.atoms if not x.is_hydrogen()]
            for k, atom in enumerate(heavy + [x for x in mol.atoms if x.is_hydrogen()]):
                atom.id = base + k
            radicals.append(next(x.id for x in heavy if x.radical_electrons))
            rows.update(rm.molecule_rows(mol))
            base += len(mol.atoms)
        return rows, radicals

    class Witness:
        participants = ()

        def __init__(self, r, p):
            self.r, self.p = r, p

        def reactant_rows(self):
            return self.r

        def product_rows(self):
            return self.p

    def count(template, smiles):
        result = wm.match([template], te.melt_state(rows_of(smiles)[0]))
        stats = result.classes[template.name]
        assert stats.pair_count in (None, len(result.transitions))
        return len(result.transitions)

    total = 0.0
    for name, (pair, d, ssa) in pairs.items():
        r, (a, b) = rows_of(pair)
        p = {x: list(v) for x, v in r.items()}
        for x, y in ((a, b), (b, a)):
            edges = dict(p[x][5])
            edges[y] = 1.0
            p[x][5] = tuple(sorted(edges.items()))
            p[x][3] -= 1
        p = {x: tuple(v) for x, v in p.items()}
        template = te.cut_template(Witness(r, p), set(r), family="", pathway=(),
                                   direction="forward").matcher(1.0, name)
        total += ssa * d / count(template, pair) * count(template, state_smiles)
    return total


def test_melt_matcher_equals_the_rate_law(crotyl_records, allyl_records):
    """Hand values (N = 4): crotyl d = 1 (head-tail), 0.5, 0.5 gives
    (1 + 0.5 + 0.5) * 12 = 24; allyl d = 2 gives 24.  The melt matcher's
    windows see every resonance form, so it reaches them, with ssa 2 on the
    head-to-tail channel (with ssa 1 it would give 18)."""
    records, _forms = crotyl_records
    by = {k: (r.degeneracy, r.ssa_multiplier) for k, r in records.items()}
    assert by == {HEAD_TAIL: (1.0, 2.0), "C=CC(C)C(C)C=C": (0.5, 2.0),
                  "CC=CCCC=CC": (0.5, 2.0)}
    state = [CROTYL[1]] * 2 + [CROTYL[0]] * 2
    crotyl = {"head-tail": ((CROTYL[1], CROTYL[0]), 1.0, 2.0),
              "head-head": ((CROTYL[1], CROTYL[1]), 0.5, 2.0),
              "tail-tail": ((CROTYL[0], CROTYL[0]), 0.5, 2.0)}
    assert _melt_total(crotyl, state) == pytest.approx(24.0, rel=1e-12)
    assert _melt_total({"allyl": (("[CH2]C=C", "[CH2]C=C"), 2.0, 2.0)},
                       ["[CH2]C=C"] * 4) == pytest.approx(24.0, rel=1e-12)


def test_limit_legacy_ssa_counts_localized_forms_only(crotyl_records):
    """Stated limit of legacy ``ssa.py``: its matcher compares a record's
    reactant graph with the state's localized bonds exactly, so a crotyl
    molecule plays only the role of the form it is stored in.  On two head
    forms and two tail forms it counts head-tail 2 * 2 = 4 (ssa 2: 8),
    head-head 1 pair and tail-tail 1 pair (ssa 2, d 0.5: 1 each): 10, against
    the rate law's 24 (which the melt matcher reaches).  No pair-count
    convention can close this; it needs resonance-aware matching."""
    records, forms = crotyl_records
    graphs = [forms[CROTYL[1]]] * 2 + [forms[CROTYL[0]]] * 2
    for seed in SEEDS:
        out = _legacy(list(records.values()), graphs, seed)
        assert sum(a for _n, a in out) == pytest.approx(10.0, rel=1e-12)
