"""Round-75 findings, measured end to end. One check per finding, each ending at
``cantera.Solution``.

Why ``cantera.Solution`` and never deck text, and never ``ck2yaml.convert_mech`` alone:
``convert_mech`` writes YAML without running ``Kinetics::checkDuplicates``, so the measured
signature of this whole defect class is *convert_mech ACCEPTED; Solution REJECTED*. Three of
the five findings below produce exactly that signature.

Why ``transport_model=None``: these fixtures carry no transport data, and the phases these
writers emit declare a transport model, so ``GasTransport::getTransportData`` aborts before
the kinetics are ever built. That rejection has nothing to do with duplicates, and reporting
it as one is the same defect this ticket is about -- a failure attributed to the stage that
was asked rather than the stage that threw. Kinetics, and therefore the duplicate check,
still run.

Run it against the build as committed to see the findings RED; run it after the repairs to
see them green. Exit status is 0 only when every check passes.
"""
import os
import sys
import tempfile

import cantera as ct

from rmgpy.species import Species
from rmgpy.molecule import Molecule
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.reaction import Reaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.kinetics import (Arrhenius, MultiArrhenius, MultiPDepArrhenius, PDepArrhenius,
                            ThirdBody)
from rmgpy.chemkin import mark_duplicate_reaction, chemkin_duplicate_flags, save_chemkin_file

TMP = tempfile.mkdtemp(prefix="r75-")
FAILURES = []


def _species(label, index, smiles):
    spc = Species(label=label, molecule=[Molecule(smiles=smiles)])
    spc.index = index
    coeffs = [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, -11.7246]
    spc.thermo = NASA(
        polynomials=[NASAPolynomial(coeffs=coeffs, Tmin=(200, "K"), Tmax=(1000, "K")),
                     NASAPolynomial(coeffs=coeffs, Tmin=(1000, "K"), Tmax=(6000, "K"))],
        Tmin=(200, "K"), Tmax=(6000, "K"))
    return spc


H = _species("H", 1, "[H]")
OH = _species("OH", 2, "[OH]")
H2 = _species("H2", 3, "[H][H]")
O = _species("O", 4, "[O]")
AR = _species("Ar", 5, "[Ar]")
HE = _species("He", 6, "[He]")
QUAD = [H, OH, H2, O]


def arr(a=1.0e12):
    return Arrhenius(A=(a, "cm^3/(mol*s)"), n=0.0, Ea=(0.0, "kcal/mol"))


def production_mark(reactions):
    """The marking loop ``rmgpy.rmg.model`` runs as the model grows, in its shape."""
    checked = []
    for rxn in reactions:
        mark_duplicate_reaction(rxn, checked)
        checked.append(rxn)


def load(path):
    try:
        ct.Solution(path, transport_model=None)
        return None
    except Exception as exc:
        lines = [ln.strip() for ln in str(exc).splitlines()
                 if ln.strip() and set(ln.strip()) != {"*"}]
        return " | ".join(lines[:3])


def check(label, reason):
    """``reason`` is None when the mechanism loaded, else the rejection."""
    if reason is None:
        print("    PASS  {0}".format(label))
    else:
        print("    FAIL  {0}\n          {1}".format(label, reason))
        FAILURES.append(label)


def writer1(species, reactions, name):
    from rmgpy.yaml_cantera1 import write_cantera
    path = os.path.join(TMP, name)
    elements = {atom.element for s in species for atom in s.molecule[0].atoms}
    write_cantera(species, reactions, elements_in_use=elements, path=path)
    return path


def writer2(species, reactions, name):
    from rmgpy.rmg.model import ReactionModel
    from rmgpy.yaml_cantera2 import save_cantera_model
    path = os.path.join(TMP, name)
    save_cantera_model(ReactionModel(species=list(species), reactions=list(reactions)), path)
    return path


def chemkin_then_cantera(species, reactions, name):
    """The production Chemkin path: write the deck, convert it, then LOAD it."""
    from cantera import ck2yaml
    deck = os.path.join(TMP, name + ".inp")
    out = os.path.join(TMP, name + ".yaml")
    save_chemkin_file(deck, species, reactions, verbose=False)
    ck2yaml.convert_mech(deck, out_name=out, quiet=True, permissive=True)
    return out


print("=" * 78)
print("F1 (HIGH 1) -- CanteraWriter1 must recompute the answer over its own list")
print("=" * 78)
core = LibraryReaction(reactants=[H, OH], products=[H2, O], library="LibA",
                       kinetics=arr(1e12), reversible=False)
edge = LibraryReaction(reactants=[H, OH], products=[H2, O], library="LibB",
                       kinetics=arr(2e12), reversible=False)
production_mark([core, edge])
assert [core.duplicate, edge.duplicate] == [True, True], (
    "precondition: production marking must mark both, or F1 never reaches the defect")
print("  production flags over core+edge:", [core.duplicate, edge.duplicate])
print("  group answer for the CORE-ONLY list:", chemkin_duplicate_flags([core]))
check("Writer1, core-only export after production marking",
      load(writer1(QUAD, [core], "f1_w1.yaml")))
check("Writer2, core-only export after production marking",
      load(writer2(QUAD, [core], "f1_w2.yaml")))

print()
print("=" * 78)
print("F2 (HIGH 2) -- a cross-class pair through the PAIRWISE production path")
print("=" * 78)
lib = LibraryReaction(reactants=[H, OH], products=[H2, O], library="LibA",
                      kinetics=arr(1e12), reversible=False)
tpl = TemplateReaction(reactants=[H, OH], products=[H2, O], family="FamA",
                       kinetics=arr(2e12), reversible=False)
production_mark([lib, tpl])
print("  production (pairwise) flags:", [lib.duplicate, tpl.duplicate])
print("  group answer over the same list:", chemkin_duplicate_flags([lib, tpl]))
if [lib.duplicate, tpl.duplicate] != [True, True]:
    print("    FAIL  pairwise marking still splits a cross-class pair the group key joins")
    FAILURES.append("F2 pairwise cross-class")
else:
    print("    PASS  pairwise marking agrees with the group key")
check("Writer1, cross-class pair", load(writer1(QUAD, [lib, tpl], "f2_w1.yaml")))
check("Writer2, cross-class pair", load(writer2(QUAD, [lib, tpl], "f2_w2.yaml")))
check("Chemkin -> ck2yaml -> Solution, cross-class pair",
      load(chemkin_then_cantera(QUAD, [lib, tpl], "f2_ck")))

print()
print("=" * 78)
print("F3 (HIGH 3) -- a Multi wrapper holding exactly ONE entry")
print("=" * 78)
for label, kin in (
        ("MultiArrhenius([one])", MultiArrhenius(arrhenius=[arr(1e12)])),
        ("MultiPDepArrhenius([one])", MultiPDepArrhenius(arrhenius=[PDepArrhenius(
            pressures=([0.1, 10.0], "bar"), arrhenius=[arr(1e12), arr(2e12)])]))):
    lone = LibraryReaction(reactants=[H, OH], products=[H2, O], library="LibA",
                           kinetics=kin, reversible=False)
    print("  --", label, "-- leaves:", len(kin.arrhenius))
    check("{0}: Chemkin -> ck2yaml -> Solution".format(label),
          load(chemkin_then_cantera(QUAD, [lone], "f3_ck")))
    check("{0}: Writer2".format(label), load(writer2(QUAD, [lone], "f3_w2.yaml")))
    check("{0}: Writer1".format(label), load(writer1(QUAD, [lone], "f3_w1.yaml")))

# A two-entry wrapper must still be marked -- the repair must not swing the other way.
pair = LibraryReaction(reactants=[H, OH], products=[H2, O], library="LibA",
                       kinetics=MultiArrhenius(arrhenius=[arr(1e12), arr(2e12)]),
                       reversible=False)
print("  -- control: MultiArrhenius([two]) must STILL be marked")
check("MultiArrhenius([two]): Chemkin -> ck2yaml -> Solution",
      load(chemkin_then_cantera(QUAD, [pair], "f3_ck2")))
check("MultiArrhenius([two]): Writer2", load(writer2(QUAD, [pair], "f3_w2b.yaml")))
check("MultiArrhenius([two]): Writer1", load(writer1(QUAD, [pair], "f3_w1b.yaml")))

print()
print("=" * 78)
print("F4 (MEDIUM) -- a collider the equation builder cannot render")
print("=" * 78)


def plog():
    return PDepArrhenius(pressures=([0.1, 10.0], "bar"), arrhenius=[arr(1e12), arr(2e12)])


a = LibraryReaction(reactants=[H, OH], products=[H2, O], library="LibA",
                    kinetics=plog(), reversible=False)
a.specific_collider = AR
b = LibraryReaction(reactants=[H, OH], products=[H2, O], library="LibA",
                    kinetics=plog(), reversible=False)
b.specific_collider = HE
print("  group answer:", chemkin_duplicate_flags([a, b]))
for name, writer in (("Writer2", writer2), ("Writer1", writer1)):
    try:
        path = writer([H, OH, H2, O, AR, HE], [a, b], "f4_{0}.yaml".format(name))
    except Exception as exc:
        print("    PASS  {0} refuses a PLOG collider it cannot render: {1}".format(
            name, type(exc).__name__))
        continue
    reason = load(path)
    print("    FAIL  {0} emitted a mechanism instead of refusing ({1})".format(
        name, reason or "and it loaded, with the two colliders conflated"))
    FAILURES.append("F4 {0}".format(name))

# The control: a ThirdBody collider IS renderable -- Writer2 writes it explicitly on both
# sides, which Cantera accepts and which Chemkin cannot express at all. The refusal above
# must not take this with it.
tb_a = Reaction(index=1, reactants=[H, H], products=[H2], reversible=False,
                specific_collider=AR,
                kinetics=ThirdBody(arrheniusLow=Arrhenius(
                    A=(1e14, "cm^6/(mol^2*s)"), n=0.0, Ea=(0.0, "kcal/mol"))))
tb_b = Reaction(index=2, reactants=[H, H], products=[H2], reversible=False,
                specific_collider=HE,
                kinetics=ThirdBody(arrheniusLow=Arrhenius(
                    A=(2e14, "cm^6/(mol^2*s)"), n=0.0, Ea=(0.0, "kcal/mol"))))
print("  -- control: ThirdBody colliders ARE renderable and must stay distinct")
check("Writer2, two ThirdBody reactions with different colliders",
      load(writer2([H, H2, AR, HE], [tb_a, tb_b], "f4_tb.yaml")))

print()
print("=" * 78)
if FAILURES:
    print("FAILURES ({0}): {1}".format(len(FAILURES), "; ".join(FAILURES)))
    sys.exit(1)
print("all checks passed")
