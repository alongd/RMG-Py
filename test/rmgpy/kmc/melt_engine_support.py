"""Declarations, toy templates and catalogs for the S4 melt-engine tests.

Not collected by pytest (no ``Test`` suffix).  Two polymers are declared:

* PS (``ps_spec``): the meltStateTest declaration plus ``transport_model.met``.
  Its numbers are the values met.py carries for the reference arm
  ``A0_REF_H_CROSS_NEc`` (WLF package H anchored at 473 K, crossover at
  N_e = 14800 g/mol / M0, <R^2>/M = 0.434 A^2 mol/g, contact radius
  SIGMA_CONTACT); they enter the engine only through this declaration.
* PE (``pe_spec``): a hand-declared linear polyethylene, -CH2-CH2- repeat
  unit (b = 2), H caps, with TOY transport numbers (declared for the test,
  not physical, not fitted) and a hand-written toy catalog
  (``pe_templates``): homolysis, end-radical recombination, H abstraction,
  end/end disproportionation, mid-radical beta scission, methyl + end
  recombination (a declared small species), and two branch-creating classes (mid + end
  recombination, intra- and intermolecular) that must accrue as leak.
"""

import json
from pathlib import Path

from rmgpy.kmc import met
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc.window_matcher import AtomEdit, BondEdit, Template, WindowAtom
from rmgpy.polymer_spec import LinearHomopolymerSpec

from meltStateTest import ps_mapping

FIXTURES = Path(__file__).resolve().parent / "fixtures"


# --------------------------------------------------------------------------
# PS
# --------------------------------------------------------------------------
def ps_met():
    """PS transport declaration: met.py's A0_REF_H_CROSS_NEc reference arm."""
    return {"schema": me.MET_SCHEMA,
            "segment_diffusivity": {"law": "wlf", "D_ref_m2_s": met._D_FERRY, "T_ref_K": 473.0,
                                    "C1": 13.39, "C2_K": 45.28, "T_wlf_K": 378.15},
            "chain_scaling": {"law": "crossover", "entanglement_units": 14800.0 / met.PS_M0},
            "mean_square_end_to_end_per_molar_mass_m2_mol_per_g": met.PS_C_R2,
            "contact_radius_m": met.SIGMA_CONTACT, "spin_factor": 1.0}


def ps_spec(counts, mobile=None, transport=True, met_block=None):
    mapping = ps_mapping(counts, mobile=mobile if mobile is not None else {})
    if transport:
        mapping["transport_model"] = {"met": met_block or ps_met()}
    return LinearHomopolymerSpec.from_mapping(mapping)


def payloads(*names):
    out = []
    for name in names:
        out += json.loads((FIXTURES / name).read_text())["payloads"]
    return out


PS_CATALOG = ("i151_ps_dp8_sampled_catalog.json", "i151_ps_dp4_hab_sampled_catalog.json")


def ps_catalog():
    return me.PayloadCatalog(payloads(*PS_CATALOG))


# --------------------------------------------------------------------------
# PE (toy)
# --------------------------------------------------------------------------
def _atom(key, element="C", hydrogens=0, radical=0):
    return {"key": key, "element": element, "hydrogens": hydrogens,
            "electronic_state": {"radical_electrons": radical}}


def pe_met(d_ref=1.0e-17):
    return {"schema": me.MET_SCHEMA,
            "segment_diffusivity": {"law": "arrhenius", "D_ref_m2_s": d_ref, "T_ref_K": 700.0,
                                    "Ea_J_mol": 0.0},
            "chain_scaling": {"law": "rouse"},
            "mean_square_end_to_end_per_molar_mass_m2_mol_per_g": 1.0e-20,
            "contact_radius_m": 5.0e-10, "spin_factor": 1.0}


def pe_mapping(counts, *, transport=True, d_ref=1.0e-17, small=None, mobile=None):
    repeat = {
        "atoms": [_atom("c1", hydrogens=2), _atom("c2", hydrogens=2)],
        "bonds": [{"atoms": ["c1", "c2"], "order": 1}],
        "ports": {"head": {"atom": "c1", "bond_order": 1}, "tail": {"atom": "c2", "bond_order": 1}},
        "electronic_state": {"multiplicity": 1},
    }

    def cap():
        return {"atoms": [_atom("h", "H")], "bonds": [],
                "ports": {"cap": {"atom": "h", "bond_order": 1}},
                "electronic_state": {"multiplicity": 1}}

    mapping = {
        "schema_version": 1, "label": "toy linear polymer (CH2-CH2)", "repeat_unit": repeat,
        "end_groups": {"head": cap(), "tail": cap()},
        "small_species": dict(small or {}),
        "initial_ensemble": {"dp_counts": dict(counts)},
        "initial_mobile_counts": dict(mobile or {}),
        "density_basis": "whole_initial_mixture" if mobile else "dry_polymer",
        "initial_mass_kg": 1.0, "density_kg_m3": 850.0,
        "locality_policy": {"radius": 2, "bounds": {
            "elements": ["C", "H"], "isotopes": [-1], "formal_charges": [0],
            "radical_electrons": [0, 1], "lone_pairs": [0],
            "bond_orders": [1, 1.5, 2, 3], "spin_multiplicities": [1, 2, 3],
            "maximum_heavy_atom_degree": 4,
            "maximum_explicit_hydrogens_per_heavy_atom": 4,
            "maximum_total_hydrogens_on_free_small_species": 12,
            "maximum_family_arity": 2, "maximum_root_tuple_size": 4,
            "maximum_free_small_species_heavy_atoms": 8,
            "attachment_classes": [], "topology_metadata_schema": {}}},
        "family_manifest": {}, "library_manifest": {},
        "transport_model": {"met": pe_met(d_ref)} if transport else {}, "diagnostics": {},
    }
    return mapping


def pe_spec(counts, **kw):
    return LinearHomopolymerSpec.from_mapping(pe_mapping(counts, **kw))


METHYL = {"graph": {"atoms": [_atom("c", hydrogens=3, radical=1)], "bonds": [],
                    "electronic_state": {"multiplicity": 2}}, "disposition": "mobile"}

MID = ((1.0, "C", False), (1.0, "H", False), (1.0, "H", False))  # a CH2 port seen from inside
PE_WEIGHTS = {"homolysis": 0.2, "recombination": 1.0e-2, "h_abstraction": 5.0e-6,
              "disproportionation": 5.0e-5, "mid_end_recombination": 2.0e-2,
              "intra_mid_end": 0.3, "beta_scission": 5.0, "methyl_recombination": 1.0e-2}


def _ch2_radical(p, port=MID):
    """C(radical) closed + 2 H closed + one CH2 port: atoms [C, H, H, port]."""
    return [WindowAtom("C", p, radical=1), WindowAtom("H", p), WindowAtom("H", p),
            WindowAtom("C", p, "port", port)]


def pe_templates(temperature=700.0, weights=None, scale=1.0):
    """The toy PE catalog; weights scale linearly with ``scale`` (used to
    emulate a temperature dependence in the refresh test)."""
    w = dict(PE_WEIGHTS, **(weights or {}))
    w = {k: v * scale for k, v in w.items()}
    out = []
    # homolysis of a C-C bond between two CH2 whose outer neighbours are CH2
    atoms = [WindowAtom("C"), WindowAtom("C"), WindowAtom("H"), WindowAtom("H"),
             WindowAtom("H"), WindowAtom("H"), WindowAtom("C", 0, "port", MID),
             WindowAtom("C", 0, "port", MID)]
    bonds = [(0, 1, 1), (0, 2, 1), (0, 3, 1), (1, 4, 1), (1, 5, 1), (0, 6, 1), (1, 7, 1)]
    out.append(Template("homolysis", "toy_homolysis", atoms, bonds,
                        [BondEdit(0, 1, 0), AtomEdit(0, radical=1), AtomEdit(1, radical=1)],
                        weight=w["homolysis"]))
    # recombination of two end radicals
    atoms = _ch2_radical(0) + _ch2_radical(1)
    bonds = [(0, 1, 1), (0, 2, 1), (0, 3, 1), (4, 5, 1), (4, 6, 1), (4, 7, 1)]
    out.append(Template("recombination", "toy_recombination", atoms, bonds,
                        [BondEdit(0, 4, 1), AtomEdit(0, radical=0), AtomEdit(4, radical=0)],
                        weight=w["recombination"]))
    # H abstraction: end radical + mid CH2 (both neighbours CH2)
    atoms = _ch2_radical(0) + [WindowAtom("C", 1), WindowAtom("H", 1), WindowAtom("H", 1),
                               WindowAtom("C", 1, "port", MID), WindowAtom("C", 1, "port", MID)]
    bonds = [(0, 1, 1), (0, 2, 1), (0, 3, 1), (4, 5, 1), (4, 6, 1), (4, 7, 1), (4, 8, 1)]
    out.append(Template("h_abstraction", "toy_h_abstraction", atoms, bonds,
                        [BondEdit(4, 5, 0), BondEdit(0, 5, 1), AtomEdit(0, radical=0),
                         AtomEdit(4, radical=1)], weight=w["h_abstraction"]))
    # disproportionation: end radical takes a beta H of another end radical
    atoms = _ch2_radical(0) + [WindowAtom("C", 1, radical=1), WindowAtom("H", 1), WindowAtom("H", 1),
                               WindowAtom("C", 1), WindowAtom("H", 1), WindowAtom("H", 1),
                               WindowAtom("C", 1, "port", MID)]
    bonds = [(0, 1, 1), (0, 2, 1), (0, 3, 1), (4, 5, 1), (4, 6, 1), (4, 7, 1), (7, 8, 1),
             (7, 9, 1), (7, 10, 1)]
    out.append(Template("disproportionation", "toy_disproportionation", atoms, bonds,
                        [BondEdit(7, 8, 0), BondEdit(0, 8, 1), BondEdit(4, 7, 2),
                         AtomEdit(0, radical=0), AtomEdit(4, radical=0)],
                        weight=w["disproportionation"]))
    # mid radical + end radical recombination (creates a branch: leak)
    atoms = [WindowAtom("C", 0, radical=1), WindowAtom("H", 0), WindowAtom("C", 0, "port", MID),
             WindowAtom("C", 0, "port", MID)] + _ch2_radical(1)
    bonds = [(0, 1, 1), (0, 2, 1), (0, 3, 1), (4, 5, 1), (4, 6, 1), (4, 7, 1)]
    out.append(Template("mid_end_recombination", "toy_branching", atoms, bonds,
                        [BondEdit(0, 4, 1), AtomEdit(0, radical=0), AtomEdit(4, radical=0)],
                        weight=w["mid_end_recombination"]))
    # the same inside one molecule (unimolecular, disconnected window: leak)
    atoms = [WindowAtom("C", 0, radical=1), WindowAtom("H", 0), WindowAtom("C", 0, "port", MID),
             WindowAtom("C", 0, "port", MID)] + _ch2_radical(0)
    out.append(Template("intra_mid_end", "toy_branching", atoms, bonds,
                        [BondEdit(0, 4, 1), AtomEdit(0, radical=0), AtomEdit(4, radical=0)],
                        weight=w["intra_mid_end"]))
    # beta scission of a mid radical: C(.)H-CH2-CH2- -> CH=CH2 + .CH2-
    atoms = [WindowAtom("C", 0, radical=1), WindowAtom("H"), WindowAtom("C", 0, "port", MID),
             WindowAtom("C"), WindowAtom("H"), WindowAtom("H"), WindowAtom("C"), WindowAtom("H"),
             WindowAtom("H"), WindowAtom("C", 0, "port", MID)]
    bonds = [(0, 1, 1), (0, 2, 1), (0, 3, 1), (3, 4, 1), (3, 5, 1), (3, 6, 1), (6, 7, 1),
             (6, 8, 1), (6, 9, 1)]
    out.append(Template("beta_scission", "toy_beta_scission", atoms, bonds,
                        [BondEdit(3, 6, 0), BondEdit(0, 3, 2), AtomEdit(0, radical=0),
                         AtomEdit(6, radical=1)], weight=w["beta_scission"]))
    # declared small species: methyl radical + end radical
    atoms = [WindowAtom("C", 0, radical=1), WindowAtom("H", 0), WindowAtom("H", 0),
             WindowAtom("H", 0)] + _ch2_radical(1)
    bonds = [(0, 1, 1), (0, 2, 1), (0, 3, 1), (4, 5, 1), (4, 6, 1), (4, 7, 1)]
    out.append(Template("methyl_recombination", "toy_recombination", atoms, bonds,
                        [BondEdit(0, 4, 1), AtomEdit(0, radical=0), AtomEdit(4, radical=0)],
                        weight=w["methyl_recombination"]))
    return out


class ToyCatalog:
    """``templates(T)`` from a builder; records every requested temperature."""

    def __init__(self, builder=pe_templates, names=None):
        self.builder = builder
        self.names = names
        self.requested = []

    def templates(self, temperature):
        self.requested.append(float(temperature))
        out = self.builder(temperature)
        if self.names is not None:
            out = [t for t in out if t.name in self.names]
        return out


# --------------------------------------------------------------------------
# state helpers
# --------------------------------------------------------------------------
def chain_atoms(state, ordinal):
    return sorted(a for a, o in state.origin.items() if o.kind == "chain" and o.ordinal == ordinal)


def carbons(state, ordinal):
    return [a for a in chain_atoms(state, ordinal) if state.atoms[a].element == "C"]


def strip_h(state, carbon):
    """Break one C-H bond of ``carbon``: the carbon and the released H each
    carry one radical electron (atoms are never deleted)."""
    from rmgpy.kmc.melt_state import SetBond, SetRadical
    h = min(b for b in state.atoms[carbon].edges if state.atoms[b].element == "H")
    state.apply([SetBond(carbon, h, 1.0, 0.0), SetRadical(carbon, state.atoms[carbon].radical,
                                                          state.atoms[carbon].radical + 1),
                 SetRadical(h, 0, 1)])
    return h


def end_carbons(state, ordinal):
    """The two terminal carbons (CH3) of a PE chain, by backbone degree."""
    cs = carbons(state, ordinal)
    return [c for c in cs if sum(1 for b in state.atoms[c].edges if state.atoms[b].element == "C") == 1]


def witness_state(smiles):
    """A MeltState of the given molecules (explicit H, atom ids 0..N-1)."""
    from rmgpy.molecule import Molecule
    from rmgpy.kmc import rmg_moves as rm
    from rmgpy.kmc import template_emission as te
    rows = {}
    for s in smiles:
        mol = Molecule(smiles=s)
        offset = len(rows)
        for i, atom in enumerate(mol.atoms):
            atom.id = offset + i
        rows.update(rm.molecule_rows(mol))
    return te.melt_state(rows)


PE2_RADICAL = "CC(c1ccccc1)C[CH]c1ccccc1"   # the pe2_identical witness radical


def pe2_two_state(temperature=700.0):
    """(catalog, state, derived forward template payload) of the pe2
    two-state system: two identical benzylic radicals <-> their dimer."""
    cat = payloads("i151_pe2_catalog.json")
    state = witness_state([PE2_RADICAL, PE2_RADICAL])
    return me.PayloadCatalog(cat), state, cat
