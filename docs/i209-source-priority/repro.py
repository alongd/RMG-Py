"""
i209 reproduction: "the sourced value wins" is an accident of ordering AND is
not even permanent, because edge pruning erases the incumbent.

Run in rmg_env from the repo root against the committed testing_database:
    python docs/i209-source-priority/repro.py

Three demonstrations, all through the real model path (make_new_species dedup +
check_for_existing_reaction), none needing a full rmg.py run:

  CASE 1  Order decides the winner. The predicate keeps whichever reaction was
          registered FIRST, ignoring kinetics and provenance. Library first ->
          sourced value survives; a library-keyed estimate first -> estimate
          survives. (A USER-SUPPLIED seed mechanism produces exactly such
          library-keyed reactions; note this is NOT the auto-generated --restart
          path, which reloads family estimates as family-keyed TemplateReactions
          that the seed sweep skips -- see report.md, Correction 1.)

  CASE 2  Pruning erases a library incumbent (the headline). A sourced library
          reaction registered first is deleted from the global registry when one
          of its edge species is pruned, after which a re-proposed reaction faces
          NO incumbent -- so a later family estimate installs as the operative
          rate.

  CASE 3  Library-vs-library collisions are real, confirmed through the identity
          path the model actually uses (Species.__eq__ is reference equality, so
          two libraries only collide after make_new_species unifies their
          species).
"""
import os

import rmgpy.data.rmg
from rmgpy import settings
from rmgpy.data.base import ForbiddenStructures
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.kinetics import Arrhenius
from rmgpy.rmg.model import CoreEdgeReactionModel, are_identical_species_references
from rmgpy.species import Species

path = os.path.join(settings["test_data.directory"], "testing_database")
kdb = KineticsDatabase()
kdb.load(os.path.join(path, "kinetics"),
         families=["H_Abstraction", "Disproportionation"],
         libraries=["GRI-Mech3.0", "ethane-oxidation"])
db = RMGDatabase()
db.kinetics = kdb
db.forbidden_structures = ForbiddenStructures()
rmgpy.data.rmg.database = db

SOURCED = Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(5.0, "kcal/mol"))   # 1.0e6 SI
ESTIMATE = Arrhenius(A=(9.9e9, "cm^3/(mol*s)"), n=0.5, Ea=(30.0, "kcal/mol"))  # 9900 SI


def species():
    a = Species().from_smiles("[H]")
    b = Species().from_smiles("C=C[CH2]C")
    c = Species().from_smiles("C=C=CC")
    d = Species().from_smiles("[H][H]")
    b.generate_resonance_structures()
    return a, b, c, d


def library(rs, ps, lib, kin):
    return LibraryReaction(reactants=list(rs), products=list(ps), library=lib,
                           kinetics=kin, reversible=True)


def template(rs, ps, kin):
    return TemplateReaction(reactants=list(rs), products=list(ps),
                            family="H_Abstraction", template=["Csd", "H"], kinetics=kin)


# ---- CASE 1: order decides the winner ----
def model_with(*spc):
    m = CoreEdgeReactionModel()
    for s in spc:
        m.add_species_to_core(s)
    return m


a, b, c, d = species()
m = model_with(a, b, c, d)
lib = library([a, b], [c, d], "GRI-Mech3.0", SOURCED)
lib.reactants.sort(); lib.products.sort()
m.add_reaction_to_core(lib); m.register_reaction(lib)
fam = template([a, b], [c, d], ESTIMATE)
fam.reactants.sort(); fam.products.sort()
found, kept = m.check_for_existing_reaction(fam)
print("CASE1a library-first: are_identical=%s found=%s kept_A=%s (SOURCED=%s)"
      % (are_identical_species_references(fam, lib), found, kept.kinetics.A.value_si, SOURCED.A.value_si))

a, b, c, d = species()
m = model_with(a, b, c, d)
est = library([a, b], [c, d], "ethane-oxidation", ESTIMATE)   # a library-keyed estimate, registered first
est.reactants.sort(); est.products.sort()
m.add_reaction_to_core(est); m.register_reaction(est)
srcd = library([a, b], [c, d], "GRI-Mech3.0", SOURCED)
srcd.reactants.sort(); srcd.products.sort()
found, kept = m.check_for_existing_reaction(srcd)
print("CASE1b estimate-first: found=%s kept_A=%s (ESTIMATE=%s, SOURCED discarded)"
      % (found, kept.kinetics.A.value_si, ESTIMATE.A.value_si))


# ---- CASE 2: pruning erases a library incumbent (headline) ----
m = CoreEdgeReactionModel()
a = m.make_new_species(Species().from_smiles("[H]"), generate_thermo=False)[0]
b = m.make_new_species(Species().from_smiles("C=C[CH2]C"), generate_thermo=False)[0]
c = m.make_new_species(Species().from_smiles("C=C=CC"), generate_thermo=False)[0]
d = m.make_new_species(Species().from_smiles("[H][H]"), generate_thermo=False)[0]
for s in (a, b, c, d):
    m.add_species_to_edge(s)
srcd = library([a, b], [c, d], "GRI-Mech3.0", SOURCED)
incumbent, _ = m.make_new_reaction(srcd, generate_thermo=False, generate_kinetics=False)


def key_of(rxn):
    for f in m.reaction_dict:
        for k1 in m.reaction_dict[f]:
            for k2 in m.reaction_dict[f][k1]:
                if rxn in m.reaction_dict[f][k1][k2]:
                    return f
    return None


e = template([a, b], [c, d], ESTIMATE); e.reactants.sort(); e.products.sort()
print("CASE2 before prune: registry_key=%s found=%s" % (key_of(incumbent), m.check_for_existing_reaction(e)[0]))
m.remove_species_from_edge([], b)
e2 = template([a, b], [c, d], ESTIMATE); e2.reactants.sort(); e2.products.sort()
print("CASE2 after prune:  registry_key=%s found=%s  (incumbent erased; a re-proposed estimate would face no collision)"
      % (key_of(incumbent), m.check_for_existing_reaction(e2)[0]))


# ---- CASE 3: library-vs-library collisions, identity-confirmed ----
m = CoreEdgeReactionModel()
for rxn in kdb.libraries["GRI-Mech3.0"].get_library_reactions():
    m.make_new_reaction(rxn, generate_thermo=False, generate_kinetics=False)
collisions = 0
for rxn in kdb.libraries["ethane-oxidation"].get_library_reactions():
    _, is_new = m.make_new_reaction(rxn, generate_thermo=False, generate_kinetics=False)
    collisions += (not is_new)
print("CASE3 confirmed library-vs-library collisions (GRI first, ethane-oxidation offered): %d of 18"
      % collisions)
