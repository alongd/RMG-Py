"""Prototype reproduction for i209: does 'sourced value wins' survive ordering?"""
import os
from rmgpy import settings
from rmgpy.data.base import ForbiddenStructures
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.data.rmg import RMGDatabase
from rmgpy.kinetics import Arrhenius
from rmgpy.molecule import Molecule
from rmgpy.rmg.main import RMG
from rmgpy.rmg.model import (
    CoreEdgeReactionModel,
    are_identical_species_references,
    get_family_library_object,
)
from rmgpy.species import Species

path = os.path.join(settings["test_data.directory"], "testing_database")
rmg = RMG()
rmg.database = RMGDatabase()
rmg.database.load_kinetics(
    os.path.join(path, "kinetics"),
    kinetics_families=["H_Abstraction"],
    reaction_libraries=["GRI-Mech3.0", "ethane-oxidation"],
)
for fam in rmg.database.kinetics.families.values():
    fam.forbidden = ForbiddenStructures()
rmg.database.forbidden_structures = ForbiddenStructures()

print("families:", list(rmg.database.kinetics.families.keys()))
print("libraries:", list(rmg.database.kinetics.libraries.keys()))


def mk_species():
    a = Species().from_smiles("[H]")
    b = Species().from_smiles("C=C[CH2]C")
    c = Species().from_smiles("C=C=CC")
    d = Species().from_smiles("[H][H]")
    a.label, b.label, c.label, d.label = "[H]", "C=C[CH2]C", "C=C=CC", "[H][H]"
    b.generate_resonance_structures()
    return a, b, c, d


SOURCED = Arrhenius(A=(1.0e12, "cm^3/(mol*s)"), n=0.0, Ea=(5.0, "kcal/mol"))
ESTIMATE = Arrhenius(A=(9.9e9, "cm^3/(mol*s)"), n=0.5, Ea=(30.0, "kcal/mol"))


def fresh_model(*species):
    cerm = CoreEdgeReactionModel()
    for s in species:
        cerm.add_species_to_core(s)
    return cerm


# ---- CASE 1: normal ordering (library incumbent, family newcomer) ----
a, b, c, d = mk_species()
cerm = fresh_model(a, b, c, d)
lib = LibraryReaction(reactants=[a, b], products=[c, d], library="GRI-Mech3.0",
                      kinetics=SOURCED, reversible=True)
lib.reactants.sort(); lib.products.sort()
cerm.add_reaction_to_core(lib)
cerm.register_reaction(lib)

fam = TemplateReaction(reactants=[a, b], products=[c, d], family="H_Abstraction",
                       template=["Csd", "H"], kinetics=ESTIMATE)
fam.reactants.sort(); fam.products.sort()
print("\nCASE1 are_identical:", are_identical_species_references(fam, lib))
found, kept = cerm.check_for_existing_reaction(fam)
print("CASE1 found:", found, "kept is library:", kept is lib,
      "kept kinetics A:", kept.kinetics.A.value_si if kept else None)

# ---- CASE 2: inverted ordering (estimate incumbent, sourced newcomer) ----
a, b, c, d = mk_species()
cerm = fresh_model(a, b, c, d)
# seed carrying a family-ESTIMATED rate, saved as a library reaction (restart-from-seed)
seed = LibraryReaction(reactants=[a, b], products=[c, d], library="ethane-oxidation",
                       kinetics=ESTIMATE, reversible=True)
seed.reactants.sort(); seed.products.sort()
cerm.add_reaction_to_core(seed)
cerm.register_reaction(seed)

libnew = LibraryReaction(reactants=[a, b], products=[c, d], library="GRI-Mech3.0",
                         kinetics=SOURCED, reversible=True)
libnew.reactants.sort(); libnew.products.sort()
print("\nCASE2 are_identical:", are_identical_species_references(libnew, seed))
found, kept = cerm.check_for_existing_reaction(libnew)
print("CASE2 found:", found, "kept is seed(estimate):", kept is seed,
      "kept kinetics A:", kept.kinetics.A.value_si if kept else None,
      "SOURCED A:", SOURCED.A.value_si, "ESTIMATE A:", ESTIMATE.A.value_si)
