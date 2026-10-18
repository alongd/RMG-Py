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
    compile."""
    from rmgpy.kmc import compiler as cp
    stub = cp.EventSetCompiler.__new__(cp.EventSetCompiler)
    stub.reference_thermo_provider = SimpleNamespace(provenance={"test": "stub"})
    grid = [300.0, 1000.0]
    proxy = cp.SiteProxy(site_type="end", reactants=list(reaction.reactants))
    return stub._record(proxy, reaction, {"test": True}, orientation,
                        rate_override=({"T": grid, "k": [1.0, 1.0],
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
