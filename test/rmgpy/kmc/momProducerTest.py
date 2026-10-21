"""Red tests: the method-of-moments producer for the runner's ``--mom`` file
(``rmgpy.kmc.mom_producer``, I-169).  No RMG run: the MoM side is a
hand-built solver profile, pool table, species dictionary and annotated
mechanism with known answers.

Expected values are derived here without the producer: the moments of the
declared ensemble by hand, the ledger-booked condensed mass written out, the
piecewise-linear interpolant of samples of t^2 in closed form, and the
profile column of a long moment-dummy label by the MoM writer's rule.  The
producer must export a key only under kMC's definition of it: at t = 0 every
exported quantity equals the kMC observation of the same deck.
"""

import csv
import hashlib
import json
import re

import numpy as np
import pytest

import melt_engine_support as S
from meltRunTest import pe_catalog_file, write_deck, write_catalog, toy_payloads

from rmgpy import constants
from rmgpy.kmc import melt_engine as me
from rmgpy.kmc import melt_observables as mo
from rmgpy.kmc import mom_producer as mp
from rmgpy.kmc import run as kr
from rmgpy.molecule import Molecule
from rmgpy.species import Species
from rmgpy.polymer_spec import LinearHomopolymerSpec

CHECKOUT = "/nonexistent/mom-checkout"   # never run here; only recorded
PE_FAMILIES = ["toy_beta_scission", "toy_branching", "toy_disproportionation", "toy_h_abstraction",
               "toy_homolysis", "toy_recombination"]
PE_LIVE = ("homolysis", "recombination", "h_abstraction", "disproportionation", "beta_scission",
           "mid_end_recombination", "intra_mid_end")
LABEL = "toylinearp"                      # mom_label of "toy linear polymer (CH2-CH2)"
DAUGHTER = LABEL + "_rad_primary_end"     # a spawned pool with a long label
BASE = {10: 1, 20: 1}


def test_refused_sidecar_is_not_treated_as_missing(tmp_path):
    path = tmp_path / "polymer_pools.json"
    path.write_text(json.dumps({"refused": True, "reason": "core_topology_identity",
                                "species": ["bad(7)"]}))
    with pytest.raises(mp.ProducerError, match="refused:core_topology_identity.*bad\\(7\\)"):
        mp.read_pools(path)


def test_reconstruction_refuses_malformed_species_by_name():
    mixed = Species(molecule=[Molecule().from_smiles("C")])
    mixed.molecule.append(Molecule().from_smiles("[CH3]"))
    mixed.label = "mixed"
    doc = {"pools": [{"label": "pool", "monomer_mw_g_mol": 1.0}],
           "conventions": {"condensed_species": []}}
    with pytest.raises(mp.ProducerError, match="reconstructed phase: mixed.*malformed"):
        mp.core_phases(["time", "volume", "mixed"], doc, [mixed], "profile.csv")


def test_reconstruction_refuses_undecidable_species_by_name():
    class UncopyableStructure:
        multiplicity = 1

        def copy(self, deep=False):
            raise ValueError("structure cannot be copied")

    undecidable = Species(label="undecidable")
    undecidable.molecule = [UncopyableStructure()]
    doc = {"pools": [{"label": "pool", "monomer_mw_g_mol": 1.0}],
           "conventions": {"condensed_species": []}}
    with pytest.raises(mp.ProducerError, match="reconstructed phase: undecidable.*resonance"):
        mp.core_phases(["time", "volume", "undecidable"], doc, [undecidable], "profile.csv")


GOLDEN = '''polymer(
    label='toylinearp',
    monomer="""multiplicity 3
1 *2 C u1 p0 c0 {2,S} {3,S} {4,S}
2 *1 C u1 p0 c0 {1,S} {5,S} {6,S}
3    H u0 p0 c0 {1,S}
4    H u0 p0 c0 {1,S}
5    H u0 p0 c0 {2,S}
6    H u0 p0 c0 {2,S}
""",
    end_groups=['[H]', '[H]'],
    cutoff=4,
    moments=[2.3651094350938804, 35.47664152640821, 591.2773587734702],
    initial_mass=0.9952320765335197,
)
phase = polymer_phase(label='toylinearp_melt', species=['toylinearp'], solvent='toylinearp', density=(850.0, 'kg/m^3'), residence={'*': 'melt', 'H': 'melt', 'H2': 'melt', 'N2': 'gas'})
hybridPolymerReactor(
    temperature=(700.0, 'K'),
    pressure=(1.0, 'bar'),
    initialMoles={'toylinearp': 2.3651094350938804, 'N2': 1.0},
    polymerPhase=phase,
    terminationTime=(40.0, 's'),
    constant_gas_volume=False,
)
'''

ADJ = {
    "N2": "1 N u0 p1 c0 {2,T}\n2 N u0 p1 c0 {1,T}\n",
    "H": "multiplicity 2\n1 H u1 p0 c0\n",
    "H2": "1 H u0 p0 c0 {2,S}\n2 H u0 p0 c0 {1,S}\n",
    "dummy": "1 Ne u0 p4 c0\n",
}


def spec_of(counts, **kw):
    return LinearHomopolymerSpec.from_mapping(S.pe_mapping(counts, **kw))


def masses(counts=BASE):
    spec = spec_of(counts)
    return spec.repeat_unit.composition.molar_mass_g_mol(), spec.end_composition.molar_mass_g_mol()


def mw(smiles):
    return Molecule(smiles=smiles).get_molecular_weight() * 1000.0


def make_setup(tmp_path, counts=BASE, *, edit=None, **kw):
    """``edit``: a function that changes the declaration mapping in place."""
    mapping = S.pe_mapping(counts, **kw)
    if edit is not None:
        edit(mapping)
    deck = write_deck(tmp_path / "pe_deck.py", mapping)
    catalog = pe_catalog_file(tmp_path)
    (tmp_path / "db" / "input" / "kinetics" / "families").mkdir(parents=True, exist_ok=True)
    return {"deck": deck, "catalog": catalog, "db": str(tmp_path / "db"), "tmp": tmp_path,
            "spec": LinearHomopolymerSpec.from_mapping(mapping)}


@pytest.fixture
def setup(tmp_path):
    return make_setup(tmp_path)


def produce(setup, out, **kw):
    args = dict(checkout=CHECKOUT, temperature=700.0, end=40.0, snapshots=9, families=PE_FAMILIES,
                catalogs=[setup["catalog"]], out=out, database=setup["db"])
    args.update(kw)
    return mp.produce(setup["deck"], **args)


# --------------------------------------------------------------------------
# a hand-built MoM run
# --------------------------------------------------------------------------
def proxy_of(label=LABEL):
    return f"{label}(3)"


def family_rows(admitted, label=LABEL):
    """{family: (rows, admitted rows)} -> pool-coupled reaction rows, each
    with its own stoichiometry (proxy <=> proxy + X<k>)."""
    rows = []
    for family, (n, ok) in sorted(admitted.items()):
        for j in range(n):
            rows.append({"family": family, "admitted": j < ok, "reactants": [proxy_of(label)],
                         "products": [proxy_of(label), f"X{len(rows)}"], "reversible": True})
    return rows


def write_run(out, *, times, pools, gas=(), admitted=None, rows=None, label=LABEL, end_groups=("[H]", "[H]"),
              mass_override=None, duplicate=None, extra_columns=(), n2=lambda t: 1.0, gas_reactions=(),
              pool_extra=None, proxy_amount=lambda t: 7.0, ea_units="KCAL/MOLE", enthalpy=None, bookkeeping=None,
              condensed=None, spawned=None, proxy_index=3, proxy_first=False, solver_mask=None,
              pool_order=None, mask_reactor=1, mask_profile=None, mask_structures=None, mask_record=None,
              species_comments=None):
    """A MoM output directory under ``out/rmg``.

    ``pools``: [(label, f(t) -> (mu0, mu1, mu2), defect_g_mol)]; ``gas``:
    [(identifier, adjacency, f(t) -> mol)]; ``admitted``: {family: (rows,
    admitted rows)} of pool-coupled reactions (default one admitted
    toy_homolysis row); ``rows``: the reaction rows themselves (dicts with
    family, admitted, reactants, products, reversible and optionally
    cantera_null / expanded: no Cantera index / a MultiArrhenius the Chemkin
    writer expands without its family annotation; A: the rate coefficient's
    A, default 1; no_chemkin: the row has no Chemkin entry; rate: the
    Chemkin main line's (A, n, Ea); aux: Chemkin lines after it (PLOG, LOW,
    ...); terms: the main lines of an expanded row; equation: the Chemkin
    equation (colliders); sidecar_kinetics: the sidecar's kinetics record,
    None as the serializer writes for non-Arrhenius kinetics);
    ``gas_reactions``: Chemkin-only reactions (reactants, products,
    reversible, A).  The Chemkin file is complete for RMG's own reader
    (SPECIES block, every reaction species in the dictionary, a THERM block:
    every species Cp = 5/2 R with H/R = -745.375 K, or ``enthalpy``'s value
    {identifier: H/R in K}); its Ea unit is ``ea_units``.  The sidecar's
    conventions.condensed_species is MoM's (every pool's mu-dummies by full
    label and the main proxy) plus ``condensed``; ``bookkeeping``: {pool:
    labels} replacing a pool's bookkeeping species; ``spawned``: the
    conventions.spawned_pools.  The SPECIES block lists the core as RMG's
    Chemkin writer does: index-less species first, then by index; the
    profile header in MoM's core order (``proxy_first``: the main proxy,
    ``label(proxy_index)``, before the gas columns).  ``solver_mask``: the
    sidecar records the solver's gas mask (conventions.solver_gas_mask) over
    the profile's species columns, each gas unless the dict {column: gas}
    says otherwise; a list of (name, gas) is written as given.
    The record is bound as MoM binds it: reactor ``mask_reactor``, profile
    ``mask_profile`` (default the profile written, solver/simulation_1_8.csv),
    each species' structure its dictionary entry's adjacency list and
    multiplicity unless ``mask_structures`` {name: structure} says otherwise;
    ``mask_record``: the conventions.solver_gas_mask value written as given.
    ``pool_order``: the sidecar's pool rows in this label order.
    ``species_comments``: {identifier: comment} on its SPECIES line, as RMG's
    verbose writer gives the species' own label(index)."""
    rmg = out / "rmg"
    for sub in ("solver", "chemkin"):
        (rmg / sub).mkdir(parents=True, exist_ok=True)
    a, _e = masses()
    moment_ids = {p: [mp.moment_column(f"{p}_mu{k}") for k in range(3)] for p, _f, _d in pools}
    proxy = f"{label}({proxy_index})"
    header = ["Time (s)", "Volume (m^3)", "N2"] + [c for p, _f, _d in pools for c in moment_ids[p]] + \
        ([proxy] if proxy_first else []) + [g[0] for g in gas] + ([] if proxy_first else [proxy]) + list(extra_columns)
    table = []
    for t in times:
        row = [t, 1.0, n2(t)]
        for _p, f, _d in pools:
            row += list(f(t))
        gas_amounts = [g[2](t) for g in gas]
        row += (([proxy_amount(t)] + gas_amounts) if proxy_first else (gas_amounts + [proxy_amount(t)])) \
            + [0.0] * len(extra_columns)
        table.append(row)
    if duplicate is not None:
        k, delta = duplicate
        bad = list(table[k])
        bad[3] += delta
        table.insert(k + 1, bad)
    with open(rmg / "solver" / "simulation_1_8.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(table)
    with open(rmg / "solver" / "simulation_1_7.csv", "w", newline="") as handle:   # an earlier core: ignored
        csv.writer(handle).writerow(["Time (s)", "Volume (m^3)"])
    if rows is None:
        rows = family_rows({"toy_homolysis": (1, 1)} if admitted is None else admitted, label)
    dictionary = [f"N2\n{ADJ['N2']}"]
    for p, _f, _d in pools:
        dictionary += [f"{c}\n{ADJ['dummy']}" for c in moment_ids[p]]
    dictionary += [f"{g[0]}\n{g[1]}" for g in gas]
    dictionary.append(f"{proxy}\n" + Molecule(smiles="CCCCCC").to_adjacency_list())
    named = {entry.split("\n", 1)[0] for entry in dictionary}
    every = [x for r in list(rows) + list(gas_reactions) for x in list(r["reactants"]) + list(r["products"])]
    every += [x for r in rows for x in r.get("colliders", ())]
    # RMG's Chemkin identifiers fit 16 characters: a longer indexed label is written as formula(index)
    ck = {x: (re.sub(r"^.*\((\d+)\)$", r"C2H2(\1)", x) if len(x) > 16 else x) for x in every}
    for x in dict.fromkeys(every):
        if ck[x] not in named:              # a structure for RMG's reader; no profile column
            dictionary.append(f"{ck[x]}\n" + Molecule(smiles="C#C").to_adjacency_list())
            named.add(ck[x])
    (rmg / "chemkin" / "species_dictionary.txt").write_text("\n".join(dictionary) + "\n")
    pool_rows = [{"label": p, "monomer_mw_g_mol": a if mass_override is None else mass_override,
                  "chain_mass_defect_g_mol": d, "end_groups": list(end_groups),
                  "bookkeeping_species": (bookkeeping or {}).get(p, moment_ids[p] + ([proxy] if p == label else [])),
                  **(pool_extra or {})}
                 for p, _f, d in pools]
    reactions, annotated, k, c = [], [f"REACTIONS    {ea_units}   MOLES", ""], 0, 0
    pool_labels = {p for p, _f, _d in pools}

    def is_proxy(x):
        return x.rsplit("(", 1)[0] in pool_labels and x.endswith(")")

    for r in rows:
        arrow = "<=>" if r["reversible"] else "=>"
        equation = " + ".join(r["reactants"]) + f" {arrow} " + " + ".join(r["products"])
        if r.get("cantera_null"):
            row = {"id": f"d{r['family']}:{equation}:0", "cantera": None}
        else:
            row = {"id": f"r{k}", "cantera": {"index": k, "equation": equation}}
            k += 1
        big_a = r.get("A", 1.0)
        kinetics = r.get("sidecar_kinetics", {"A": big_a, "n": 0.0, "Ea": 0.0})
        row.update({"kinetics": None if kinetics is None else {**kinetics, "reversible": r["reversible"]},
                    "reactants": list(r["reactants"]), "products": list(r["products"]),
                    "proxy_reactants": [x for x in r["reactants"] if is_proxy(x)],
                    "proxy_products": [x for x in r["products"] if is_proxy(x)],
                    "refused": not r["admitted"], "archetype": r.get("archetype", "same_pool/1"),
                    **r.get("extra", {})})
        if not r["admitted"]:
            row["refused_reason"] = "conduit-deferred"
        reactions.append(row)
        line = r.get("equation") or ("+".join(ck[x] for x in r["reactants"]) + arrow
                                     + "+".join(ck[x] for x in r["products"]))
        terms = r.get("terms") or [r.get("rate", (big_a, 0.0, 0.0))] * (2 if r.get("expanded") else 1)
        for term in ([] if r.get("no_chemkin") else terms):
            c += 1
            annotated.append(f"! Reaction index: Chemkin #{c}; RMG #{len(reactions)}")
            if not r.get("expanded"):
                annotated.append(f"! Template reaction: {r['family']}")
            annotated += ([f"{line:<51} {term[0]!r} {term[1]!r} {term[2]!r}"] + list(r.get("aux", ()))
                          + (["DUPLICATE"] if r.get("expanded") else []) + [""])
    for j, g in enumerate(gas_reactions):
        c += 1
        line = ("+".join(ck[x] for x in g["reactants"]) + ("<=>" if g["reversible"] else "=>")
                + "+".join(ck[x] for x in g["products"]))
        annotated += [f"! Reaction index: Chemkin #{c}; RMG #{len(reactions) + 1 + j}",
                      "! Template reaction: toy_h_abstraction", f"{line:<51} {g['A']} 0.0 0.0", ""]
    if pool_order is not None:
        pool_rows = [next(r for r in pool_rows if r["label"] == x) for x in pool_order]
    if isinstance(solver_mask, dict):
        solver_mask = [(h, solver_mask.get(h, True)) for h in header[2:]]
    adjacency = {entry.split("\n", 1)[0]: entry.split("\n", 1)[1] for entry in dictionary}

    def structure(name):
        if name in (mask_structures or {}):
            return mask_structures[name]
        if name not in adjacency:
            return None
        return {"adjacency_list": adjacency[name],
                "multiplicity": Molecule().from_adjacency_list(adjacency[name]).multiplicity}
    recorded = {}
    if solver_mask is not None:
        recorded = {"solver_gas_mask": {"source": "solver", "reactors": [{
            "reactor": mask_reactor, "profile": mask_profile or "solver/simulation_1_8.csv",
            "core": [{"species": x, "structure": structure(x), "gas": g} for x, g in solver_mask]}]}}
    if mask_record is not None:
        recorded = {"solver_gas_mask": mask_record}
    (rmg / "chemkin" / "polymer_pools.json").write_text(json.dumps(
        {"schema_version": "2.0", "pools": pool_rows, "reactions": reactions,
         "conventions": {"stale_topology": False, **recorded,
                         "condensed_species": sorted([proxy] + [f"{p}_mu{k}" for p, _f, _d in pools for k in range(3)]
                                                     + list(condensed or ())),
                         **({"spawned_pools": list(spawned)} if spawned else {})}}))
    order = sorted(named, key=lambda x: (mp.species_key(x)[0] == "index", mp.species_key(x)[1:]))
    species_block = ["SPECIES"] + [f"    {x}" + (f"    ! {species_comments[x]}" if x in (species_comments or {}) else "")
                                   for x in order] + ["END", ""]
    therm = ["THERM ALL", "   300.000  1000.000  5000.000", ""]
    for x in sorted(named):
        h = f"{(enthalpy or {}).get(x, -745.375):15.8E}"
        therm += [f"{x:<24}Ne  1               G   200.000  6000.000 1000.00      1",
                  " 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2",
                  f"{h} 3.35532000E+00 2.50000000E+00 0.00000000E+00 0.00000000E+00    3",
                  f" 0.00000000E+00 0.00000000E+00{h} 3.35532000E+00                   4", ""]
    therm += ["END", ""]
    (rmg / "chemkin" / "chem_annotated.inp").write_text("\n".join(species_block + therm + annotated + ["END", ""]))


def mom_entries(out):
    """The entries of a produced --mom file (header form)."""
    doc = json.loads((out / "mom.json").read_text())
    assert set(doc) == {"header", "entries"}
    return doc["entries"]


WHOLE = ("whole_mu0_mol", "whole_mu1_mol", "whole_mu2_mol", "whole_Mn_g_mol", "whole_Mw_g_mol")


def declared(spec):
    return [float(x) for x in spec.macroscopic_moments_mol]


ZERO_GAS = [("H(1)", ADJ["H"], lambda t: 0.0), ("H2(2)", ADJ["H2"], lambda t: 0.0)]


def drift(mu):
    """Pool moments whose chain count grows (scission-like) at fixed mu1: the
    profile moves, the booked mass (ledger) stays 1000 g."""
    return lambda t: (mu[0] * (1 + 1e-3 * t), mu[1], mu[2])


def constant_run(out, spec, times=(0.0, 10.0, 20.0, 40.0, 45.0), flat=False, **kw):
    """A run whose mass stays put; ``flat``: nothing in the profile moves."""
    mu = declared(spec)
    write_run(out, times=list(times), pools=[(LABEL, (lambda t: mu) if flat else drift(mu), 0.0)],
              gas=ZERO_GAS, **kw)


PROFILE_T = [0.0, 10.0, 20.0, 40.0, 45.0]       # the solver's last step passes the end time


def interpolant(f, t):
    """The piecewise-linear interpolant of f sampled at PROFILE_T."""
    for lo, hi in zip(PROFILE_T, PROFILE_T[1:]):
        if lo <= t <= hi:
            return f(lo) + (f(hi) - f(lo)) * (t - lo) / (hi - lo)
    raise AssertionError(t)


def moving_run(out, spec, **kw):
    """Main pool loses mass to H2 and chains to a spawned daughter pool;
    total mass closes exactly (linear in the columns, so also between rows)."""
    mu = declared(spec)
    a, _e = masses()
    m_h2 = mw("[H][H]")
    d, h = 1e-3, 2e-3
    r = 5 * d + m_h2 * h / a

    def main(t):
        return mu[0] * (1 + 1e-3 * t * t), mu[1] - r * t * t, mu[2] * (1 + 1e-3 * t * t)

    def daughter(t):
        return d * t * t, 5 * d * t * t, 40 * d * t * t

    write_run(out, times=PROFILE_T, pools=[(LABEL, main, 0.0), (DAUGHTER, daughter, 0.0)],
              gas=[("H(1)", ADJ["H"], lambda t: 0.0), ("H2(2)", ADJ["H2"], lambda t: h * t * t)], **kw)
    return main, daughter, h, m_h2


# --------------------------------------------------------------------------
# deck -> input file
# --------------------------------------------------------------------------
def test_deck_gives_the_golden_mom_input(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    text = (out / "rmg" / "input.py").read_text()
    start = text.index("polymer(\n")
    assert text[start:text.index("model(")] == GOLDEN
    a, e = masses()
    pair = (10 * a + e) + (20 * a + e)
    mu0 = 2 * 1000.0 / pair
    assert [2.3651094350938804, 35.47664152640821, 591.2773587734702] == pytest.approx(
        [mu0, 15 * mu0, 250 * mu0], rel=1e-14)
    assert 0.9952320765335197 == pytest.approx(a * 15 * mu0 / 1000.0, rel=1e-14)
    assert f"kineticsFamilies={PE_FAMILIES!r}," in text
    assert "reactionLibraries=[]," in text and "seedMechanisms=[]," in text
    for kernel in ("k_scission", "k_unzip", "k_homolysis", "radical_qssa_unzip", "side_group_homolysis",
                   "explicit_dp=True"):
        assert kernel not in text
    for name, smiles in (("H", "[H]"), ("H2", "[H][H]")):
        assert f"species(label={name!r}, reactive=True, structure=SMILES({smiles!r}))" in text
    assert "species(label='N2', reactive=False, structure=SMILES('N#N'))" in text
    assert "saveSimulationProfiles=True" in text
    assert (out / "rmg" / "rmgrc").read_text() == f"database.directory : {setup['db']}/input\n"


def test_conventions_are_written_and_name_every_setting(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    doc = json.loads((out / "conventions.json").read_text())
    names = {c["name"]: c for c in doc["conventions"]}
    for name in ("label", "legacy_translation", "initial_moments", "mass_booking", "mass_closure",
                 "tail_keys", "volume", "declared_species", "bath_gas", "pressure", "families",
                 "family_check", "database", "kernels_off", "no_cage", "model_settings", "profile",
                 "time_interpolation", "species_constraints", "phase_kinetics", "inactive_chemistry",
                 "initial_inventory"):
        assert name in names and names[name]["statement"], name
    assert names["time_interpolation"]["value"] == "linear in time"
    assert names["label"]["value"] == LABEL
    assert set(doc["quantity_sources"]) | set(doc["not_produced"]) == mp.runner_quantity_keys()


# --------------------------------------------------------------------------
# item 1: a key is exported only under kMC's definition of it
# --------------------------------------------------------------------------
def test_tail_keys_are_not_produced():
    for key in ("mu0_mol", "mu1_mol", "mu2_mol", "Mn_g_mol", "Mw_g_mol"):
        assert key in mp.NOT_PRODUCED and key not in mp.QUANTITY_SOURCES
    assert set(mp.QUANTITY_SOURCES) == {"polymer_condensed_mass_g", "polymer_conversion",
                                        "polymer_derived_molecules_mol", "species_mol", *WHOLE}
    assert set(mp.QUANTITY_SOURCES).isdisjoint(mp.NOT_PRODUCED)
    assert set(mp.QUANTITY_SOURCES) | set(mp.NOT_PRODUCED) == mp.runner_quantity_keys()


@pytest.mark.parametrize("counts", [{2: 1, 8: 1}, BASE])
def test_every_exported_quantity_equals_the_kmc_value_of_the_same_deck_at_t0(tmp_path, counts):
    """Spar's fixture (DP2 + DP8 with xs = 4: the kMC tail holds one chain, the
    MoM pool both): no exported row can disagree at t = 0 by definition."""
    s = make_setup(tmp_path, counts)
    out = tmp_path / "out"
    produce(s, out)
    constant_run(out, s["spec"])
    produce(s, out, mode="convert")
    quantities = mom_entries(out)[0]["quantities"]
    engine = me.MeltEngine(s["spec"], S.ToyCatalog(names=PE_LIVE), 700.0, np.random.default_rng(1))
    obs = mo.MeltObserver(s["spec"], xs=4, identity=engine.identity).observe(engine)
    rows = [r for r in mo.side_by_side(obs, quantities) if r["mom"] is not None]
    assert {r["quantity"] for r in rows} >= {"polymer_condensed_mass_g", "polymer_derived_molecules_mol",
                                             "polymer_conversion", "species_mol.H", "species_mol.H2", *WHOLE}
    assert all(r["agreement"] == "equal" for r in rows), [r for r in rows if r["agreement"] != "equal"]


# --------------------------------------------------------------------------
# item 2: the condensed mass is MoM's ledger booking
# --------------------------------------------------------------------------
def test_scission_does_not_invent_end_group_mass(tmp_path):
    """PE DP6 -> two DP3 (chain count doubles, repeat units kept): the mass
    stays 1000 g (adding e per chain would give 1011.8 g)."""
    s = make_setup(tmp_path, {6: 1})
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    mu = declared(s["spec"])

    def split(t):
        return mu[0] * (1 + t / 10.0), mu[1], mu[2] * (1 - 0.5 * t / 10.0)

    write_run(out, times=[0.0, 10.0], pools=[(LABEL, split, 0.0)], gas=ZERO_GAS)
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    entries = mom_entries(out)
    a, e = masses()
    assert 1000.0 * (1 + e / (6 * a + e)) == pytest.approx(1011.8, abs=0.1)     # the invented value
    for entry in entries:
        assert entry["quantities"]["polymer_condensed_mass_g"] == pytest.approx(1000.0, rel=1e-12)
        assert entry["quantities"]["polymer_conversion"] == pytest.approx(0.0, abs=1e-12)
    assert entries[1]["quantities"]["polymer_derived_molecules_mol"] == pytest.approx(2 * mu[0], rel=1e-12)


def test_profile_converts_to_the_mom_file_at_the_observation_times(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    main, daughter, h, m_h2 = moving_run(out, setup["spec"])
    produce(setup, out, mode="convert")
    entries = mom_entries(out)
    times = kr.observation_times(end=40.0, snapshots=9)
    assert [x["time_s"] for x in entries] == times
    for entry in entries:
        t, q = entry["time_s"], entry["quantities"]
        mu0 = interpolant(lambda s: main(s)[0], t) + interpolant(lambda s: daughter(s)[0], t)
        mass = 1000.0 - m_h2 * interpolant(lambda s: h * s * s, t)
        if t == 25.0:   # closed form between the 20 s and 40 s samples of t^2: (20 + 40) 25 - 800 = 700
            assert interpolant(lambda s: s * s, t) == pytest.approx(700.0, rel=1e-15)
        assert q["polymer_derived_molecules_mol"] == pytest.approx(mu0, rel=1e-12)
        assert q["polymer_condensed_mass_g"] == pytest.approx(mass, rel=1e-12)
        assert q["polymer_conversion"] == pytest.approx(1.0 - mass / 1000.0, rel=1e-9, abs=1e-15)
        assert q["species_mol"]["H"] == 0.0
        assert q["species_mol"]["H2"] == pytest.approx(interpolant(lambda s: h * s * s, t), rel=1e-12)
        assert not set(q) & set(mp.NOT_PRODUCED)
    manifest = json.loads((out / "producer.json").read_text())
    assert manifest["read"]["profile"].endswith("simulation_1_8.csv")
    assert manifest["read"]["pool_labels"] == [LABEL, DAUGHTER]


def test_chain_mass_defect_enters_the_booked_mass(setup):
    """A feature pool whose chains lost H2 (defect = M(H2) per chain): chains
    move to it while H2 is released; the booked mass is 1000 - M(H2) y t^2."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    m_h2 = mw("[H][H]")
    y = 1e-3

    def main(t):
        return mu[0] - y * t * t, mu[1] - 15 * y * t * t, mu[2] - 250 * y * t * t

    def feature(t):
        return y * t * t, 15 * y * t * t, 250 * y * t * t

    write_run(out, times=PROFILE_T, pools=[(LABEL, main, 0.0), (LABEL + "_mod", feature, m_h2)],
              gas=[("H(1)", ADJ["H"], lambda t: 0.0), ("H2(2)", ADJ["H2"], lambda t: y * t * t)])
    produce(setup, out, mode="convert")
    for entry in mom_entries(out):
        t = entry["time_s"]
        assert entry["quantities"]["polymer_condensed_mass_g"] == pytest.approx(
            1000.0 - m_h2 * interpolant(lambda s: y * s * s, t), rel=1e-12)


def test_a_booking_that_does_not_close_exports_no_mass(setup):
    """A capped DP3 chain released to the gas phase: MoM's ledger books its end
    groups as condensed while the gas species carries them too, so the booked
    total grows by e per released chain: mass and conversion are not exported
    where the booking does not close."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    x = 1e-4                        # the pool's variance stays positive to 45 s

    def main(t):
        return mu[0] - x * t * t, mu[1] - 3 * x * t * t, mu[2] - 9 * x * t * t

    hexane = Molecule(smiles="CCCCCC").to_adjacency_list()
    write_run(out, times=PROFILE_T, pools=[(LABEL, main, 0.0)],
              gas=ZERO_GAS + [("C6H14(4)", hexane, lambda t: x * t * t)])
    produce(setup, out, mode="convert")
    entries = mom_entries(out)
    assert "polymer_condensed_mass_g" in entries[0]["quantities"]
    for entry in entries[1:]:
        assert "polymer_condensed_mass_g" not in entry["quantities"]
        assert "polymer_conversion" not in entry["quantities"]
        assert "polymer_derived_molecules_mol" in entry["quantities"]
    diagnostics = json.loads((out / "mom_diagnostics.json").read_text())
    _a, e = masses()
    assert diagnostics[-1]["mass_closure_residual_g"] == pytest.approx(e * x * 40.0 ** 2, rel=1e-6)


# --------------------------------------------------------------------------
# item 3: inactive chemistry
# --------------------------------------------------------------------------
def test_inactive_chemistry_is_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], admitted={"toy_homolysis": (17, 0), "R_Recombination": (2, 2)})
    with pytest.raises(mp.ProducerError, match="inactive_chemistry"):
        produce(setup, out, mode="convert")
    assert not (out / "mom.json").exists()
    manifest = json.loads((out / "producer.json").read_text())
    assert manifest["status"] == "refused:inactive_chemistry"
    census = manifest["census"]
    assert census["families"]["toy_homolysis"] == {"rows": 17, "admitted": 0}
    assert census["families"]["R_Recombination"] == {"rows": 2, "admitted": 2}   # not a requested family
    assert census["admitted_in_requested_families"] == 0


def test_no_pool_reaction_at_all_is_inactive(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], admitted={})
    with pytest.raises(mp.ProducerError, match="inactive_chemistry"):
        produce(setup, out, mode="convert")


def test_inactive_chemistry_allowed_is_flagged(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], admitted={"toy_homolysis": (3, 0)})
    produce(setup, out, mode="convert", allow_inactive=True)
    assert (out / "mom.json").exists()
    manifest = json.loads((out / "producer.json").read_text())
    conv = json.loads((out / "conventions.json").read_text())
    for doc in (manifest, conv):
        assert doc["flags"]["inactive_chemistry"]["admitted_in_requested_families"] == 0
    assert manifest["status"] == "converted"


def test_admitted_chemistry_converts_unflagged(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], admitted={"toy_homolysis": (3, 1)})
    manifest = produce(setup, out, mode="convert")
    assert manifest["census"]["admitted_in_requested_families"] == 1
    assert "inactive_chemistry" not in manifest["flags"]


# --------------------------------------------------------------------------
# item 4: phase kinetics
# --------------------------------------------------------------------------
def test_phase_kinetics_ratio_is_recorded_and_flagged(setup):
    out = setup["tmp"] / "out"
    manifest = produce(setup, out)
    conv = json.loads((out / "conventions.json").read_text())
    names = {c["name"]: c["value"] for c in conv["conventions"]}
    v_kmc = 1.0 / 850.0                                 # 1 kg dry polymer at 850 kg/m^3
    # H and H2 are melt-resident: MoM's concentration is n / V_poly, its pool
    # mass a mu1 (no end groups) over the density -- not 1, so still flagged
    facts = json.loads((out / "producer.json").read_text())["facts"]
    v_poly = facts["a_g_mol"] * facts["moments_mol"][1] / 1000.0 / 850.0
    ratio = v_kmc / v_poly
    assert abs(ratio - 1.0) > 1.0e-6
    kinetics = names["phase_kinetics"]
    assert set(kinetics["ratio_by_species"]) == {"H", "H2"}
    for value in kinetics["ratio_by_species"].values():
        assert value == pytest.approx(ratio, rel=1e-12)
    for doc in (manifest, conv):
        assert doc["flags"]["kinetics_not_comparable"]["ratio_by_species"]["H"] == pytest.approx(ratio, rel=1e-12)


# --------------------------------------------------------------------------
# item 5: species naming
# --------------------------------------------------------------------------
def test_a_formula_named_species_column_resolves(tmp_path):
    """RMG writes monomer_product as C2H4(4) (label + index > 16 characters):
    it is found by its structure in the species dictionary."""
    s = make_setup(tmp_path, BASE, monomer=S.ETHYLENE_GRAPH)
    out = tmp_path / "out"
    produce(s, out)
    mu = declared(s["spec"])
    ethylene = Molecule(smiles="C=C").to_adjacency_list()
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: mu, 0.0)],
              gas=ZERO_GAS + [("C2H4(4)", ethylene, lambda t: 1e-3 * t)])
    produce(s, out, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    assert q["species_mol"]["monomer_product"] == pytest.approx(0.04, rel=1e-12)
    assert set(q["species_mol"]) == {"monomer_product", "H", "H2"}


def test_a_species_label_colliding_with_the_pool_label_is_refused(tmp_path):
    s = make_setup(tmp_path, BASE, small={LABEL: S.METHYL})
    with pytest.raises(mp.ProducerError, match="collide"):
        produce(s, tmp_path / "out")
    assert not (tmp_path / "out" / "rmg" / "input.py").exists()


# --------------------------------------------------------------------------
# items 6, 7 and the profile's own checks
# --------------------------------------------------------------------------
def test_a_wrong_initial_inventory_is_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: [0.5 * x for x in mu], 0.0)], gas=ZERO_GAS)
    with pytest.raises(mp.ProducerError, match="initial inventory"):
        produce(setup, out, mode="convert")


def test_conflicting_duplicate_times_are_refused_anywhere(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: mu, 0.0)], gas=ZERO_GAS,
              duplicate=(4, 1e-3))                     # at 45 s, past every requested time
    with pytest.raises(mp.ProducerError, match="duplicate"):
        produce(setup, out, mode="convert")


def test_a_repeat_mass_that_differs_from_the_deck_is_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    a, _e = masses()
    constant_run(out, setup["spec"], mass_override=a * (1 + 1e-6))
    with pytest.raises(mp.ProducerError, match="repeat-unit mass"):
        produce(setup, out, mode="convert")


def test_a_first_row_not_at_zero_is_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], times=(-1.0, 0.0, 10.0, 20.0, 40.0, 45.0))
    with pytest.raises(mp.ProducerError, match="first row"):
        produce(setup, out, mode="convert")


def test_a_time_on_a_profile_row_takes_the_row_exactly(setup):
    out = setup["tmp"] / "out"
    produce(setup, out, end=20.0, snapshots=3)
    main, daughter, _h, _m = moving_run(out, setup["spec"])
    produce(setup, out, end=20.0, snapshots=3, mode="convert")
    entries = mom_entries(out)
    assert [x["time_s"] for x in entries] == [0.0, 10.0, 20.0]
    for entry in entries:
        t = entry["time_s"]
        assert entry["quantities"]["polymer_derived_molecules_mol"] == main(t)[0] + daughter(t)[0]


def test_a_time_outside_the_profile_is_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out, end=50.0, snapshots=3)
    constant_run(out, setup["spec"])
    with pytest.raises(mp.ProducerError, match="outside the profile"):
        produce(setup, out, end=50.0, snapshots=3, mode="convert")


def test_moment_columns_of_an_unknown_pool_are_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], extra_columns=("other_mu0",))
    with pytest.raises(mp.ProducerError, match="belong to no pool"):
        produce(setup, out, mode="convert")


def _add_pool_row(out, **row):
    """Append a pool row to the run's pool table (the MoM writer lists every
    registered pool, also one with no core species: its phase and bookkeeping
    species are then empty and the solver holds no state for it)."""
    path = out / "rmg" / "chemkin" / "polymer_pools.json"
    doc = json.loads(path.read_text())
    a, _e = masses()
    doc["pools"].append({"label": "stateless_tail", "monomer_mw_g_mol": a, "end_groups": ["[H]", "[H]"],
                         "phase_species": [], "bookkeeping_species": [], "moments": [0.0, 0.0, 0.0], **row})
    path.write_text(json.dumps(doc))


def test_a_pool_with_no_core_species_has_no_columns_and_converts(setup):
    base = setup["tmp"] / "base"
    produce(setup, base)
    constant_run(base, setup["spec"])
    produce(setup, base, mode="convert")
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"])
    _add_pool_row(out)
    produce(setup, out, mode="convert")
    assert mom_entries(out) == mom_entries(base)
    read = json.loads((out / "producer.json").read_text())["read"]
    assert read["pool_labels"] == [LABEL]
    assert read["pools_without_state"] == ["stateless_tail"]


def test_a_pool_with_no_core_species_but_moments_is_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"])
    _add_pool_row(out, moments=[1e-3, 3e-2, 1.0])
    with pytest.raises(mp.ProducerError, match="no core species but moments"):
        produce(setup, out, mode="convert")


def test_a_long_pool_label_maps_its_hashed_columns_through_the_species_comments(setup):
    """RMG's verbose Chemkin writer names a long moment dummy by its 16-character
    identifier and gives the full label in the comment, and RMG's reader takes
    the comment as the label: the profile column is the identifier."""
    base = setup["tmp"] / "base"
    produce(setup, base)
    moving_run(base, setup["spec"])
    produce(setup, base, mode="convert")
    out = setup["tmp"] / "out"
    produce(setup, out)
    comments = {mp.moment_column(f"{DAUGHTER}_mu{k}"): f"{DAUGHTER}_mu{k}" for k in range(3)}
    assert all(len(full) > 16 and column != full for column, full in comments.items())
    moving_run(out, setup["spec"], species_comments=comments)
    produce(setup, out, mode="convert")
    assert mom_entries(out) == mom_entries(base)


def test_a_pool_with_core_species_and_no_columns_is_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"])
    _add_pool_row(out, phase_species=["stateless_tail_mu0"], bookkeeping_species=["stateless_tail_mu0"])
    with pytest.raises(mp.ProducerError, match="no column .* for pool 'stateless_tail'"):
        produce(setup, out, mode="convert")


def test_moment_column_follows_the_mom_writer_rule():
    assert mp.moment_column("PS_mu0") == "PS_mu0"
    assert mp.moment_column("toylinearp_mu2") == "toylinearp_mu2"            # <= 16: verbatim
    long = DAUGHTER + "_mu1"
    assert mp.moment_column(long) == "toyli-" + hashlib.md5(long.encode()).hexdigest()[:6] + "_mu1"
    assert len(mp.moment_column(long)) == 16


# --------------------------------------------------------------------------
# quantity keys
# --------------------------------------------------------------------------
def test_quantity_keys_are_the_runners():
    keys = mp.runner_quantity_keys()
    assert keys == {"mu0_mol", "mu1_mol", "mu2_mol", "end_inventory_mol", "explicit_dp_mol",
                    "polymer_condensed_mass_g", "Mn_g_mol", "Mw_g_mol", "polymer_conversion", "released_mol",
                    "species_mol", "polymer_derived_molecules_mol", *WHOLE}
    mp.check_quantity_keys({"polymer_derived_molecules_mol": 1.0, "species_mol": {"H": 0.0}}, ["H", "H2"])
    with pytest.raises(mp.ProducerError, match="not keys of"):
        mp.check_quantity_keys({"polymer_derived_molecules_mol": 1.0, "Mn": 1.0}, [])
    with pytest.raises(mp.ProducerError, match="not produced"):
        mp.check_quantity_keys({"mu0_mol": 1.0}, [])
    with pytest.raises(mp.ProducerError, match="not declared"):
        mp.check_quantity_keys({"species_mol": {"CH4": 1.0}}, ["H", "H2"])


# --------------------------------------------------------------------------
# families
# --------------------------------------------------------------------------
def test_a_family_mismatch_is_refused_unless_allowed_and_then_recorded(setup):
    with pytest.raises(mp.ProducerError, match="differ from the catalogs' provenance families"):
        produce(setup, setup["tmp"] / "refused", families=["R_Recombination"])
    assert not (setup["tmp"] / "refused" / "rmg" / "input.py").exists()
    out = setup["tmp"] / "allowed"
    manifest = produce(setup, out, families=["R_Recombination"], allow_family_mismatch=True)
    record = manifest["family_check"]
    assert record["match"] is False and record["allowed_mismatch"] is True
    assert record["only_in_families"] == ["R_Recombination"] and record["only_in_catalogs"] == PE_FAMILIES
    conv = {c["name"]: c["value"] for c in json.loads((out / "conventions.json").read_text())["conventions"]}
    assert conv["family_check"] == record
    assert "kineticsFamilies=['R_Recombination']," in (out / "rmg" / "input.py").read_text()
    match = produce(setup, setup["tmp"] / "match", families=list(reversed(PE_FAMILIES)))
    assert match["family_check"]["match"] is True and "allowed_mismatch" not in match["family_check"]


def test_two_catalogs_give_the_union_of_their_families(setup):
    tmp = setup["tmp"]
    first = write_catalog(tmp / "a.json", toy_payloads([t for t in S.pe_templates() if t.name == "homolysis"]),
                          tag="a")
    second = write_catalog(tmp / "b.json", toy_payloads([t for t in S.pe_templates()
                                                         if t.name == "recombination"]), tag="b")
    manifest = produce(setup, tmp / "two", catalogs=[first, second],
                       families=["toy_homolysis", "toy_recombination"])
    assert manifest["family_check"]["catalog_families"] == ["toy_homolysis", "toy_recombination"]


# --------------------------------------------------------------------------
# round trip and the command line
# --------------------------------------------------------------------------
def test_a_produced_mom_file_is_accepted_by_the_runner(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    moving_run(out, setup["spec"])
    produce(setup, out, mode="convert")
    times = kr.observation_times(end=40.0, snapshots=9)
    mom, record = kr.load_mom(out / "mom.json", times, deck_sha256=kr.read_deck(setup["deck"]).sha256,
                              temperature_K=700.0)
    assert record["times"] == len(times) and all(q is not None for q in mom)


def test_command_line_writes_and_refuses(setup, capsys):
    out = setup["tmp"] / "cli"
    argv = [str(setup["deck"]), "--mom-checkout", CHECKOUT, "--temperature", "700", "--end", "40",
            "--snapshots", "9", "--families", *PE_FAMILIES, "--catalog", str(setup["catalog"]),
            "--out", str(out), "--database", setup["db"], "--write-only"]
    assert mp.main(argv) == 0
    assert (out / "rmg" / "input.py").is_file() and (out / "conventions.json").is_file()
    assert mp.main(argv) == 2                                  # a non-empty output directory
    assert "not empty" in capsys.readouterr().err
    bad = argv[:argv.index("--families") + 1] + ["R_Recombination"] + argv[argv.index("--catalog"):]
    bad[bad.index("--out") + 1] = str(setup["tmp"] / "cli2")
    assert mp.main(bad) == 2
    assert "--allow-family-mismatch" in capsys.readouterr().err
    constant_run(out, setup["spec"], admitted={"toy_homolysis": (1, 0)})
    convert = [x for x in argv if x != "--write-only"] + ["--convert"]
    assert mp.main(convert) == 2
    assert "inactive_chemistry" in capsys.readouterr().err
    assert mp.main(convert + ["--allow-inactive"]) == 0


def test_species_constraints_keep_the_declared_species_admissible(setup):
    out = setup["tmp"] / "out"
    manifest = produce(setup, out, species_constraints=["maximumRadicalElectrons=0"])
    text = (out / "rmg" / "input.py").read_text()
    assert "generatedSpeciesConstraints(allowed=['input species'], maximumRadicalElectrons=0)" in text
    assert manifest["facts"]["constraints"] == {"allowed": ["input species"], "maximumRadicalElectrons": 0}
    with pytest.raises(mp.ProducerError, match="KEY=VALUE"):
        produce(setup, setup["tmp"] / "bad", species_constraints=["maximumRadicalElectrons"])


# --------------------------------------------------------------------------
# rework 2 (r127)
# --------------------------------------------------------------------------
def test_gas_oligomers_are_booked_on_the_side_kmc_puts_them(tmp_path):
    """Spar's fixture: PE DP6 -> DP3 (kept in the pool) + hexene in MoM's gas
    phase.  kMC counts the undeclared hexene (n = 3 >= 1) as polymer-derived:
    the mass stays 1000 g, the molecules double, the conversion is 0."""
    s = make_setup(tmp_path, {6: 1})
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    u, mu1, mu2 = declared(s["spec"])

    def pool(t):
        f = t / 10.0
        return u, mu1 - 3 * f * u, mu2 * (1 - f) + 9 * f * u

    hexene = Molecule(smiles="C=CCCCC").to_adjacency_list()
    write_run(out, times=[0.0, 10.0], pools=[(LABEL, pool, 0.0)],
              gas=ZERO_GAS + [("C6H12(4)", hexene, lambda t: u * t / 10.0)])
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    assert q["polymer_condensed_mass_g"] == pytest.approx(1000.0, rel=1e-12)
    assert q["polymer_conversion"] == pytest.approx(0.0, abs=1e-12)
    assert q["polymer_derived_molecules_mol"] == pytest.approx(2 * u, rel=1e-12)
    manifest = json.loads((out / "producer.json").read_text())
    assert manifest["read"]["species_classes"]["C6H12(4)"] == "polymer"


def _methyl_capped(mapping):
    mapping["end_groups"]["head"] = {"atoms": [S._atom("m", hydrogens=3)], "bonds": [],
                                     "ports": {"cap": {"atom": "m", "bond_order": 1}},
                                     "electronic_state": {"multiplicity": 1}}


def test_a_species_kmc_cannot_classify_exports_no_polymer_rows(tmp_path):
    """Methyl-capped PE: carbon is on the backbone and in an end group, so
    whether a released propane carries a whole repeat unit's backbone (kMC:
    polymer-derived) depends on ancestry MoM does not keep.  The polymer rows
    are not exported where it is present, even though the booking closes."""
    s = make_setup(tmp_path, BASE, edit=_methyl_capped)
    out = tmp_path / "out"
    produce(s, out)
    mu = declared(s["spec"])
    a, _e = masses()
    m3 = mw("CCC")
    x = 1e-4

    def pool(t):
        return mu[0], mu[1] - m3 * x * t / a, mu[2]

    propane = Molecule(smiles="CCC").to_adjacency_list()
    write_run(out, times=PROFILE_T, pools=[(LABEL, pool, 0.0)], end_groups=("[CH3]", "[H]"),
              gas=ZERO_GAS + [("C3H8(4)", propane, lambda t: x * t)])
    produce(s, out, mode="convert")
    entries = mom_entries(out)
    assert {"polymer_condensed_mass_g", "polymer_derived_molecules_mol", *WHOLE} <= set(entries[0]["quantities"])
    diagnostics = json.loads((out / "mom_diagnostics.json").read_text())
    for entry, diag in list(zip(entries, diagnostics))[1:]:
        for key in ("polymer_condensed_mass_g", "polymer_conversion", "polymer_derived_molecules_mol", *WHOLE):
            assert key not in entry["quantities"]
        assert diag["mass_closure_residual_g"] == pytest.approx(0.0, abs=1e-9)
        assert diag["unclassified_species"] == ["C3H8(4)"]


def test_zero_flux_chemistry_is_inactive(setup):
    """17 refused homolysis rows and one admitted irreversible H abstraction
    that needs H, which nothing admitted makes: inactive.  With one homolysis
    row admitted (it makes H) both rows are live."""
    proxy = proxy_of()
    hab = {"family": "toy_h_abstraction", "admitted": True, "reactants": ["H(1)", proxy],
           "products": ["H2(2)", proxy], "reversible": False}

    def homolysis(admitted):
        return {"family": "toy_homolysis", "admitted": admitted, "reactants": [proxy],
                "products": [proxy, "H(1)"], "reversible": False}

    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], rows=[homolysis(False)] * 17 + [hab])
    with pytest.raises(mp.ProducerError, match="inactive_chemistry"):
        produce(setup, out, mode="convert")
    census = json.loads((out / "producer.json").read_text())["census"]
    assert census["admitted_in_requested_families"] == 1 and census["live_rows"] == []
    live = setup["tmp"] / "live"
    produce(setup, live)
    constant_run(live, setup["spec"], rows=[homolysis(True)] + [homolysis(False)] * 16 + [hab])
    manifest = produce(setup, live, mode="convert")
    assert "inactive_chemistry" not in manifest["flags"]
    assert manifest["census"]["live_rows"] == ["r0", "r17"]


def test_a_cantera_null_row_is_counted_by_its_serialized_family(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    proxy = proxy_of()
    constant_run(out, setup["spec"], rows=[{"family": "toy_homolysis", "admitted": True, "reactants": [proxy],
                                            "products": [proxy, "X0"], "reversible": True,
                                            "cantera_null": True}])
    manifest = produce(setup, out, mode="convert")
    assert manifest["census"]["families"]["toy_homolysis"] == {"rows": 1, "admitted": 1}


def test_an_expanded_multi_arrhenius_row_is_counted_by_name(setup):
    """The Chemkin writer expands a template reaction's MultiArrhenius into
    plain entries without the family: counted by its equation, as a row of
    the requested families (the input admits no other source)."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    proxy = proxy_of()
    constant_run(out, setup["spec"], rows=[{"family": "toy_homolysis", "admitted": True, "reactants": [proxy],
                                            "products": [proxy, "X0"], "reversible": True, "expanded": True}])
    manifest = produce(setup, out, mode="convert")
    census = manifest["census"]
    assert census["unattributed"] == {f"{proxy} <=> {proxy} + X0": {"rows": 1, "admitted": 1}}
    assert census["admitted_in_requested_families"] == 1 and "inactive_chemistry" not in manifest["flags"]


def test_the_bath_does_not_dilute_the_mass_closure(tmp_path):
    """1 ug of PE in 1 mol N2: a declared ethyl appears while the pool keeps
    its mass (the exported polymer mass would be 3.59 % high).  Against the
    bath it closes to 1e-9; against the kMC inventory it does not."""
    def tiny(mapping):
        mapping["initial_mass_kg"] = 1.0e-9

    s = make_setup(tmp_path, BASE, edit=tiny, small={"ethyl": S.ETHYL})
    out = tmp_path / "out"
    produce(s, out)
    mu = declared(s["spec"])
    b = 0.0359e-6 / mw("[CH2]C")
    ethyl = Molecule(smiles="[CH2]C").to_adjacency_list()
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: mu, 0.0)],
              gas=ZERO_GAS + [("C2H5(4)", ethyl, lambda t: b * t / 40.0)])
    produce(s, out, mode="convert")
    entries = mom_entries(out)
    assert "polymer_condensed_mass_g" in entries[0]["quantities"]
    for entry in entries[1:]:
        assert "polymer_condensed_mass_g" not in entry["quantities"]
    assert entries[-1]["quantities"]["species_mol"]["ethyl"] == pytest.approx(b, rel=1e-12)


def test_a_bath_that_does_not_close_exports_no_mass(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: mu, 0.0)], gas=TRACE,
              n2=lambda t: 1.0 - 1e-3 * t)
    produce(setup, out, mode="convert")
    entries = mom_entries(out)
    diagnostics = json.loads((out / "mom_diagnostics.json").read_text())
    for entry, diag in list(zip(entries, diagnostics))[1:]:
        assert "polymer_condensed_mass_g" not in entry["quantities"]
        assert diag["bath_closure_residual_g"] == pytest.approx(-1e-3 * entry["time_s"] * mw("N#N"), rel=1e-9)


def test_the_initial_mobile_inventory_is_checked(tmp_path):
    s = make_setup(tmp_path, BASE, small={"ethylene": S.ETHYLENE}, mobile={"ethylene": 3})
    out = tmp_path / "out"
    manifest = produce(s, out)
    amount = manifest["facts"]["initial_moles"]["ethylene"]
    assert amount > 0
    mu = declared(s["spec"])
    ethylene = Molecule(smiles="C=C").to_adjacency_list()
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: mu, 0.0)],
              gas=ZERO_GAS + [("C2H4(4)", ethylene, lambda t: 0.0)])
    with pytest.raises(mp.ProducerError, match="initial inventory"):
        produce(s, out, mode="convert")
    write_run(out, times=PROFILE_T, pools=[(LABEL, drift(mu), 0.0)],
              gas=ZERO_GAS + [("C2H4(4)", ethylene, lambda t: amount)])
    q = produce(s, out, mode="convert") and mom_entries(out)[-1]["quantities"]
    assert q["species_mol"]["ethylene"] == amount


def test_a_resonance_form_of_a_declared_species_resolves(tmp_path):
    """Declared [CH2]C=CC; RMG's dictionary holds the form C=C[CH]C."""
    s = make_setup(tmp_path, BASE, small={"crotyl": S.CROTYL})
    out = tmp_path / "out"
    produce(s, out)
    mu = declared(s["spec"])
    other = Molecule(smiles="C=C[CH]C")
    assert not other.is_isomorphic(Molecule(smiles="[CH2]C=CC"))
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: mu, 0.0)],
              gas=ZERO_GAS + [("C4H7(4)", other.to_adjacency_list(), lambda t: 1e-6 * t)])
    produce(s, out, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    assert q["species_mol"]["crotyl"] == pytest.approx(4e-5, rel=1e-12)


def test_initial_moments_off_by_a_part_per_million_are_refused(setup):
    """Between INITIAL_RELATIVE (1e-9) and a loose bound such as 0.1."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: [x * (1 + 1e-6) for x in mu], 0.0)], gas=ZERO_GAS)
    with pytest.raises(mp.ProducerError, match="initial inventory"):
        produce(setup, out, mode="convert")


def test_a_doubled_initial_mu2_is_refused(setup):
    """mu0 and mu1 (hence the booked mass) right, mu2 doubled: only the moment
    check sees it."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: [mu[0], mu[1], 2 * mu[2]], 0.0)], gas=ZERO_GAS)
    with pytest.raises(mp.ProducerError, match="initial inventory"):
        produce(setup, out, mode="convert")


def test_an_inert_species_kmc_cannot_hold_blocks_nothing_while_absent(setup):
    """RMG puts Ar in every core; the kMC melt has no argon: unclassifiable,
    harmless at zero amount."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    write_run(out, times=PROFILE_T, pools=[(LABEL, drift(mu), 0.0)],
              gas=ZERO_GAS + [("Ar", "1 Ar u0 p4 c0\n", lambda t: 0.0)])
    manifest = produce(setup, out, mode="convert")
    assert manifest["read"]["species_classes"]["Ar"] == "undetermined"
    q = mom_entries(out)[-1]["quantities"]
    assert q["polymer_condensed_mass_g"] == pytest.approx(1000.0, rel=1e-12)


# --------------------------------------------------------------------------
# rework 3 (r129)
# --------------------------------------------------------------------------
def _refused_census(setup, out):
    with pytest.raises(mp.ProducerError, match="inactive_chemistry"):
        produce(setup, out, mode="convert")
    return json.loads((out / "producer.json").read_text())["census"]


def test_a_zero_rate_row_is_dead(setup):
    """An admitted homolysis with A = 0 on a moving profile: not live."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    proxy = proxy_of()
    constant_run(out, setup["spec"], rows=[{"family": "toy_homolysis", "admitted": True, "reactants": [proxy],
                                            "products": [proxy, "X0"], "reversible": True, "A": 0.0}])
    census = _refused_census(setup, out)
    assert census["live_rows"] == [] and census["zero_rate_rows"] == ["r0"]


def test_a_zero_rate_gas_reaction_unlocks_nothing(setup):
    """The admitted H abstraction needs H; the only maker of H is a gas
    reaction from the bath: with A = 0 it makes nothing, with A = 1 it does."""
    proxy = proxy_of()
    hab = [{"family": "toy_h_abstraction", "admitted": True, "reactants": ["H(1)", proxy],
            "products": ["H2(2)", proxy], "reversible": False}]
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], rows=hab,
                 gas_reactions=[{"reactants": ["N2"], "products": ["H(1)"], "reversible": False, "A": 0.0}])
    assert _refused_census(setup, out)["live_rows"] == []
    live = setup["tmp"] / "live"
    produce(setup, live)
    constant_run(live, setup["spec"], rows=hab,
                 gas_reactions=[{"reactants": ["N2"], "products": ["H(1)"], "reversible": False, "A": 1.0}])
    assert produce(setup, live, mode="convert")["census"]["live_rows"] == ["r0"]


def test_a_flat_profile_is_inactive_whatever_the_census(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], flat=True)
    census = _refused_census(setup, out)
    assert census["live_rows"] == ["r0"] and census["profile_changes"] is False


def test_a_reversible_row_runs_from_its_products(setup):
    """Survivor reversible-not-run: products present at t = 0, reactants absent."""
    proxy = proxy_of()

    def row(reversible):
        return [{"family": "toy_homolysis", "admitted": True, "reactants": ["X5"], "products": [proxy],
                 "reversible": reversible}]

    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], rows=row(True))
    assert produce(setup, out, mode="convert")["census"]["live_rows"] == ["r0"]
    one_way = setup["tmp"] / "one_way"
    produce(setup, one_way)
    constant_run(one_way, setup["spec"], rows=row(False))
    assert _refused_census(setup, one_way)["live_rows"] == []


def test_a_proxy_of_an_empty_pool_is_absent(setup):
    """Survivor pool-mu0-zero: the daughter pool has mu0 = 0 at t = 0; a row
    that needs its proxy is not live."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    proxy = f"{DAUGHTER}(5)"
    moving_run(out, setup["spec"], rows=[{"family": "toy_homolysis", "admitted": True, "reactants": [proxy],
                                          "products": [proxy, "X0"], "reversible": False}])
    assert _refused_census(setup, out)["live_rows"] == []


def test_duplicate_stoichiometries_keep_their_own_families(setup):
    """Spar's fixture: one refused toy_homolysis row and one admitted
    toy_beta_scission row with the same species."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    proxy = proxy_of()
    rows = [{"family": family, "admitted": ok, "reactants": [proxy], "products": [proxy, "X0"],
             "reversible": True} for family, ok in (("toy_homolysis", False), ("toy_beta_scission", True))]
    constant_run(out, setup["spec"], rows=rows)
    census = produce(setup, out, mode="convert")["census"]
    assert census["families"] == {"toy_homolysis": {"rows": 1, "admitted": 0},
                                  "toy_beta_scission": {"rows": 1, "admitted": 1}}


def test_rows_without_a_one_to_one_chemkin_match_are_refused(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    proxy = proxy_of()
    rows = [{"family": "toy_homolysis", "admitted": True, "reactants": [proxy], "products": [proxy, "X0"],
             "reversible": True, "no_chemkin": j == 1} for j in range(2)]
    constant_run(out, setup["spec"], rows=rows)
    with pytest.raises(mp.ProducerError, match="one-to-one"):
        produce(setup, out, mode="convert")


def test_toluene_in_a_ps_deck_is_undetermined_and_withholds_the_polymer_rows(tmp_path):
    """Pins the conservative bracket (r129 item 3): PS carbon is backbone and
    ring, so a released toluene is light or polymer-derived by ancestry MoM
    does not keep.  A later, tighter bracket must change this test."""
    from meltStateTest import ps_mapping
    mapping = ps_mapping({10: 1, 20: 1}, mobile={})
    mapping["small_species"] = {}
    mapping["transport_model"] = {"met": S.ps_met()}
    deck = write_deck(tmp_path / "ps_deck.py", mapping)
    spec = LinearHomopolymerSpec.from_mapping(mapping)
    s = {"deck": deck, "catalog": pe_catalog_file(tmp_path), "db": str(tmp_path / "db"), "tmp": tmp_path,
         "spec": spec}
    (tmp_path / "db" / "input" / "kinetics" / "families").mkdir(parents=True, exist_ok=True)
    out = tmp_path / "out"
    produce(s, out)
    a = spec.repeat_unit.composition.molar_mass_g_mol()
    mu = declared(spec)
    m7, x = mw("Cc1ccccc1"), 1e-4
    toluene = Molecule(smiles="Cc1ccccc1").to_adjacency_list()
    write_run(out, times=PROFILE_T, label="PS", pools=[("PS", lambda t: (mu[0], mu[1] - m7 * x * t / a, mu[2]),
                                                        0.0)],
              mass_override=a, gas=ZERO_GAS + [("C7H8(4)", toluene, lambda t: x * t)])
    manifest = produce(s, out, mode="convert")
    assert manifest["read"]["species_classes"]["C7H8(4)"] == "undetermined"
    entries = mom_entries(out)
    diagnostics = json.loads((out / "mom_diagnostics.json").read_text())
    for entry, diag in list(zip(entries, diagnostics))[1:]:
        for key in ("polymer_condensed_mass_g", "polymer_conversion", "polymer_derived_molecules_mol"):
            assert key not in entry["quantities"]
        assert diag["unclassified_species"] == ["C7H8(4)"]
        assert diag["mass_closure_residual_g"] == pytest.approx(0.0, abs=1e-9)


# --------------------------------------------------------------------------
# I-170: whole-ensemble moments and the header form
# --------------------------------------------------------------------------
def n_mass(smiles, spec):
    """The mass DP (M - e) / a of a molecule, with the deck's a and e."""
    a = spec.repeat_unit.composition.molar_mass_g_mol()
    e = spec.end_composition.molar_mass_g_mol()
    return (mw(smiles) - e) / a


def kmc_quantities(spec, smiles):
    """kMC's reported quantities of a melt of the given molecules (one
    representative each; PE: every carbon has backbone origin)."""
    from rmgpy.kmc import melt_state as ms
    from rmgpy.kmc import rmg_moves as rr
    nodes, origin, k = {}, {}, 0
    for j, smi in enumerate(smiles):
        mol = Molecule(smiles=smi)
        for i, atom in enumerate(mol.atoms):
            atom.id = k + i
        for i, (el, iso, charge, radical, lone, edges) in rr.molecule_rows(mol).items():
            nodes[i] = ms.Node(el, iso, charge, radical, lone, dict(edges))
            origin[i] = ms.Origin("chain", j, "unit 1", str(i), el == "C")
        k += len(mol.atoms)
    return mo.mom_reported_quantities(mo.MeltObserver(spec, xs=4).observe(ms.MeltState(nodes, origin)))


def kmc_whole(spec, smiles):
    """kMC's whole-ensemble keys of a melt of the given molecules."""
    q = kmc_quantities(spec, smiles)
    return {key: q[key] for key in WHOLE}


def test_whole_ensemble_keys_are_the_pools_plus_polymer_class_species(tmp_path):
    """Pool moments plus the hexene (kMC class polymer) at its mass DP."""
    s = make_setup(tmp_path, {6: 1})
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    u, mu1, mu2 = declared(s["spec"])

    def pool(t):
        f = t / 10.0
        return u, mu1 - 3 * f * u, mu2 * (1 - f) + 9 * f * u

    hexene = Molecule(smiles="C=CCCCC").to_adjacency_list()
    write_run(out, times=[0.0, 10.0], pools=[(LABEL, pool, 0.0)],
              gas=ZERO_GAS + [("C6H12(4)", hexene, lambda t: u * t / 10.0)])
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    n6 = n_mass("C=CCCCC", s["spec"])
    whole = [2 * u, 3 * u + n6 * u, 9 * u + n6 * n6 * u]
    assert [q["whole_mu0_mol"], q["whole_mu1_mol"], q["whole_mu2_mol"]] == pytest.approx(whole, rel=1e-12)


def test_propene_ejection_gives_the_same_whole_moments_on_both_sides(tmp_path):
    """Spar r133: PE DP6 + DP4; the DP6 chain ejects propene and keeps 4.5
    monomer-equivalents (nonane).  MoM: pool [2, 8.5, 36.25] u plus propene;
    kMC: nonane, octane and propene by their masses."""
    s = make_setup(tmp_path, {6: 1, 4: 1})
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    mu = declared(s["spec"])
    u = mu[0] / 2

    def pool(t):
        f = t / 10.0
        return mu[0], mu[1] - 1.5 * f * u, mu[2] + f * u * (4.5 ** 2 - 36)

    propene = Molecule(smiles="C=CC").to_adjacency_list()
    write_run(out, times=[0.0, 10.0], pools=[(LABEL, pool, 0.0)],
              gas=ZERO_GAS + [("C3H6(4)", propene, lambda t: u * t / 10.0)])
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    assert "polymer_condensed_mass_g" in q                     # the booking closes
    kmc = kmc_whole(s["spec"], ["CCCCCCCCC", "CCCCCCCC", "C=CC"])
    for key in WHOLE:
        assert q[key] == pytest.approx(kmc[key], rel=1e-12), key


def test_an_h_abstraction_gives_a_fractional_mass_dp_on_both_sides(tmp_path):
    """The DP-10 chain loses one H (to the gas as H): its mass DP is
    10 - m_H / a on both sides."""
    s = make_setup(tmp_path, BASE)
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    mu = declared(s["spec"])
    u = mu[0] / 2
    a, _e = masses()
    d = mw("[H]") / a

    def pool(t):
        f = t / 10.0
        return mu[0], mu[1] - d * f * u, mu[2] + f * u * ((10 - d) ** 2 - 100)

    write_run(out, times=[0.0, 10.0], pools=[(LABEL, pool, 0.0)],
              gas=[("H(1)", ADJ["H"], lambda t: u * t / 10.0), ("H2(2)", ADJ["H2"], lambda t: 0.0)])
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    assert q["whole_mu1_mol"] / u == pytest.approx(30 - d, rel=1e-12) and abs(30 - d - 30) > 1e-3
    kmc = kmc_whole(s["spec"], ["[CH2]" + "C" * 19, "C" * 40])      # DP 10 less one H; DP 20
    for key in WHOLE:
        assert q[key] == pytest.approx(kmc[key], rel=1e-12), key


def test_a_pool_ledger_gives_mu0_mu1_and_mn_and_withholds_mu2_and_mw(tmp_path):
    """Scission DP6 -> two DP3 radicals: each carries one cap (L = -u), mass
    DP 3 - e / 2a.  Sum and count are exact from the ledger; how the missing
    caps spread over chains (mu2, Mw) is not in the profile."""
    s = make_setup(tmp_path, {6: 1})
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    mu = declared(s["spec"])

    def split(t):
        return mu[0] * (1 + t / 10.0), mu[1], mu[2] * (1 - 0.5 * t / 10.0)

    write_run(out, times=[0.0, 10.0], pools=[(LABEL, split, 0.0)], gas=ZERO_GAS)
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    a, e = masses()
    u = mu[0]
    q = mom_entries(out)[-1]["quantities"]
    assert q["whole_mu0_mol"] == pytest.approx(2 * u, rel=1e-12)
    assert q["whole_mu1_mol"] == pytest.approx(6 * u - e * u / a, rel=1e-12)
    assert q["whole_Mn_g_mol"] == pytest.approx((6 * a + e) / 2, rel=1e-12)
    assert "whole_mu2_mol" not in q and "whole_Mw_g_mol" not in q
    assert "whole_mu2_mol" in mom_entries(out)[0]["quantities"]
    reasons = json.loads((out / "mom_diagnostics.json").read_text())[-1]["whole_withheld"]
    assert set(reasons) == {"whole_mu2_mol", "whole_Mw_g_mol"} and "ledger" in reasons["whole_mu2_mol"]


def test_a_pool_with_a_chain_mass_defect_is_shifted_to_the_mass_dp(setup):
    """Feature pool chains lost H2 (defect = M(H2)): mass DP n - M(H2) / a."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    m_h2 = mw("[H][H]")
    a, _e = masses()
    y = 1e-3

    def main(t):
        return mu[0] - y * t * t, mu[1] - 15 * y * t * t, mu[2] - 250 * y * t * t

    def feature(t):
        return y * t * t, 15 * y * t * t, 250 * y * t * t

    write_run(out, times=PROFILE_T, pools=[(LABEL, main, 0.0), (LABEL + "_mod", feature, m_h2)],
              gas=[("H(1)", ADJ["H"], lambda t: 0.0), ("H2(2)", ADJ["H2"], lambda t: y * t * t)])
    produce(setup, out, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    d = m_h2 / a
    f = feature(40.0)
    main40 = main(40.0)
    want = [main40[0] + f[0], main40[1] + f[1] - d * f[0], main40[2] + f[2] - 2 * d * f[1] + d * d * f[0]]
    assert [q["whole_mu0_mol"], q["whole_mu1_mol"], q["whole_mu2_mol"]] == pytest.approx(want, rel=1e-12)


def _boundary_run(tmp_path, scale):
    """A monodisperse DP-1 deck (mu0 = mu1 = mu2 = u); at 10 s the pool row is
    ``scale`` x (u, u, u): the mass and the chain count barely move (a trace of
    H keeps the profile moving), so the booking closes and the ledger is zero,
    and only the moments' realizability is at stake."""
    s = make_setup(tmp_path, {1: 3})
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    mu = declared(s["spec"])
    row = tuple(m * f for m, f in zip(mu, scale))
    write_run(out, times=[0.0, 10.0], pools=[(LABEL, lambda t: mu if t == 0.0 else row, 0.0)], gas=TRACE)
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    return mom_entries(out)[-1]["quantities"], json.loads((out / "mom_diagnostics.json").read_text())[-1]


def test_a_boundary_roundoff_row_keeps_every_quantity(tmp_path):
    """Spar r133 item 2: (1, nextafter(1, 0), 1) x u, a state MoM accepts
    (mu1 a ulp below mu0)."""
    q, _d = _boundary_run(tmp_path, (1.0, float(np.nextafter(1.0, 0.0)), 1.0))
    assert set(WHOLE) <= set(q) and "species_mol" in q and "polymer_derived_molecules_mol" in q


def declared_row_mu2(tmp_path):
    """The unclamped mu2 of the clamp fixture's 10 s row."""
    with open(tmp_path / "out" / "rmg" / "solver" / "simulation_1_8.csv") as handle:
        rows = list(csv.reader(handle))
    return float(rows[-1][rows[0].index(mp.moment_column(f"{LABEL}_mu2"))])


def test_moments_within_the_mom_floors_are_clamped(tmp_path):
    """mu2 a few ulps below mu1^2 / mu0: inside MoM's accepted-state variance
    tolerance (polymer.pyx, floors max(1e-30, 100 atol)): clamped to zero
    variance, Mw = Mn."""
    q, _d = _boundary_run(tmp_path, (1.0, 1.0, 1.0 - 5e-16))
    assert q["whole_mu2_mol"] == q["whole_mu1_mol"] * q["whole_mu1_mol"] / q["whole_mu0_mol"]   # clamped exactly
    assert q["whole_mu2_mol"] != declared_row_mu2(tmp_path)
    assert q["whole_Mw_g_mol"] == pytest.approx(q["whole_Mn_g_mol"], rel=1e-12)


def test_a_pool_variance_beyond_the_mom_census_tolerance_is_corruption(tmp_path):
    """r137 item 4: a pool whose mu0 mu2 - mu1^2 is far below MoM's variance
    census tolerance holds no distribution: corruption, the mass and the
    whole keys are withheld, the species rows kept."""
    q, d = _boundary_run(tmp_path, (1.0, 1.0, 0.5))
    for key in ("polymer_condensed_mass_g", *WHOLE):
        assert key not in q, key
    assert "species_mol" in q and "variance" in d["polymer_withheld"]


def test_the_mom_file_carries_the_flags_and_the_conventions_digest(setup):
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], admitted={"toy_homolysis": (3, 0)})
    manifest = produce(setup, out, mode="convert", allow_inactive=True)
    doc = json.loads((out / "mom.json").read_text())
    header = doc["header"]
    assert sorted(header["flags"]) == ["inactive_chemistry", "kinetics_not_comparable"]
    assert header["flags"] == manifest["flags"]
    assert header["conventions_sha256"] == hashlib.sha256((out / "conventions.json").read_bytes()).hexdigest()
    assert header["producer"]["schema"] == mp.SCHEMA and header["producer"]["deck_sha256"] == manifest["deck"]["sha256"]
    mom, record = kr.load_mom(out / "mom.json", kr.observation_times(end=40.0, snapshots=9),
                              deck_sha256=manifest["deck"]["sha256"], temperature_K=700.0)
    assert record["header"] == header and all(q is not None for q in mom)




@pytest.mark.parametrize("kind", ["chip", "kernel"])
def test_a_mechanism_that_leaves_the_mass_dp_coordinate_withholds_the_whole_keys(setup, kind):
    """Pool-coordinate audit: an admitted discrete_chip/1 row (rounded chip
    units), or a pool kernel with a nonzero rate, takes the pool moments off
    the mass DP after t = 0."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    proxy = proxy_of()
    rows = [{"family": "toy_homolysis", "admitted": True, "reactants": [proxy], "products": [proxy, "X0"],
             "reversible": True, "archetype": "discrete_chip/1" if kind == "chip" else "same_pool/1"}]
    extra = {"channels": {"scission": {"A": 1.0, "n": 0.0, "Ea": 0.0}}} if kind == "kernel" else None
    constant_run(out, setup["spec"], rows=rows, pool_extra=extra)
    produce(setup, out, mode="convert")
    entries = mom_entries(out)
    assert set(WHOLE) <= set(entries[0]["quantities"])
    diagnostics = json.loads((out / "mom_diagnostics.json").read_text())
    for entry, diag in list(zip(entries, diagnostics))[1:]:
        assert not set(WHOLE) & set(entry["quantities"])
        assert ("discrete_chip" if kind == "chip" else "kernel scission") in diag["whole_withheld"]["whole_mu0_mol"]


# --------------------------------------------------------------------------
# I-170 rework 2 (r136)
# --------------------------------------------------------------------------
FLOOR = 1e-14            # MoM's accepted-state floor for atol 1e-16: max(1e-30, 100 atol)
# keeps a profile moving: H (declared) reaches 1e-4 g by 10 s, above the activity floors (1e-8 of the
# deck's 1000 g) and inside the mass closure (1e-6 of the system's mass)
TRACE = [("H(1)", ADJ["H"], lambda t: 1e-5 * t), ("H2(2)", ADJ["H2"], lambda t: 0.0)]


def _ethylene_setup(tmp_path):
    s = make_setup(tmp_path, {1: 1, 6: 1}, monomer=S.ETHYLENE_GRAPH)
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    return s, out, declared(s["spec"])[0] / 2


def test_a_defect_pool_matching_a_declared_species_is_withheld_not_booked(tmp_path):
    """Spar r136/r137: PE DP1 (ethane) + DP6; the DP1 chain loses H2 into the
    feature pool, pinned at DP 1 with ethylene's formula.  A defect pool keeps
    no structure, so its chains are not booked as the declared monomer by
    mass alone (r137 item 1): the polymer rows and the monomer row are
    withheld, never counted as polymer."""
    s, out, u = _ethylene_setup(tmp_path)
    m_h2 = mw("[H][H]")
    ethylene = Molecule(smiles="C=C").to_adjacency_list()

    def main(t):
        f = t / 10.0
        return 2 * u - f * u, 7 * u - f * u, 37 * u - f * u

    def feature(t):
        f = t / 10.0
        return f * u, f * u, f * u

    write_run(out, times=[0.0, 10.0], pools=[(LABEL, main, 0.0), (LABEL + "_mod", feature, m_h2)],
              gas=[("H(1)", ADJ["H"], lambda t: 0.0), ("H2(2)", ADJ["H2"], lambda t: u * t / 10.0),
                   ("C2H4(4)", ethylene, lambda t: 0.0)])
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    for key in ("polymer_condensed_mass_g", "polymer_derived_molecules_mol", *WHOLE):
        assert key not in q, key
    assert "monomer_product" not in q["species_mol"] and "H2" in q["species_mol"]
    assert "structure" in json.loads((out / "mom_diagnostics.json").read_text())[-1]["polymer_withheld"]


def test_a_pool_that_may_hold_a_declared_species_unpinned_is_withheld(tmp_path):
    """The feature pool holds DP 1 and DP 2 chains: some are ethylene, how
    many the moments do not say."""
    s, out, u = _ethylene_setup(tmp_path)
    m_h2 = mw("[H][H]")
    ethylene = Molecule(smiles="C=C").to_adjacency_list()

    def main(t):
        f = t / 10.0
        return 2 * u - f * u, 7 * u - 1.5 * f * u, 37 * u - 2.5 * f * u

    def feature(t):
        f = t / 10.0
        return f * u, 1.5 * f * u, 2.5 * f * u

    write_run(out, times=[0.0, 10.0], pools=[(LABEL, main, 0.0), (LABEL + "_mod", feature, m_h2)],
              gas=[("H(1)", ADJ["H"], lambda t: 0.0), ("H2(2)", ADJ["H2"], lambda t: u * t / 10.0),
                   ("C2H4(4)", ethylene, lambda t: 0.0)])
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    for key in ("polymer_condensed_mass_g", "polymer_conversion", "polymer_derived_molecules_mol", *WHOLE):
        assert key not in q, key
    assert "monomer_product" not in q["species_mol"] and "H2" in q["species_mol"]
    reason = json.loads((out / "mom_diagnostics.json").read_text())[-1]["polymer_withheld"]
    assert "monomer_product" in reason


def test_a_pool_whose_chains_may_lose_their_last_backbone_unit_is_withheld(setup):
    """An admitted volatile ejection of 1.5 units into the pool: chains below
    one unit's backbone (kMC: light) may be in it after t = 0."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    proxy = proxy_of()
    constant_run(out, setup["spec"], rows=[{"family": "toy_homolysis", "admitted": True, "reactants": [proxy],
                                            "products": [proxy, "X0"], "reversible": True,
                                            "archetype": "volatile_ejection/1",
                                            "extra": {"params": {"eject_units": 1.5}, "src_pool": LABEL,
                                                      "dst_pool": LABEL}}])
    produce(setup, out, mode="convert")
    entries = mom_entries(out)
    assert "polymer_condensed_mass_g" in entries[0]["quantities"]
    for entry in entries[1:]:
        assert "polymer_derived_molecules_mol" not in entry["quantities"]
    assert "light" in json.loads((out / "mom_diagnostics.json").read_text())[-1]["polymer_withheld"]


def _octyl_capped(mapping):
    cap = {"atoms": [S._atom(f"c{k}", hydrogens=3 if k == 8 else 2) for k in range(1, 9)],
           "bonds": [{"atoms": [f"c{k}", f"c{k + 1}"], "order": 1} for k in range(1, 8)],
           "ports": {"cap": {"atom": "c1", "bond_order": 1}}, "electronic_state": {"multiplicity": 1}}
    mapping["end_groups"] = {"head": cap, "tail": json.loads(json.dumps(cap))}


def test_a_negative_whole_mu1_is_valid_and_mn_is_exported(tmp_path):
    """Spar r133/r136 item 3: heavy caps (two octyls) and a DP6 cut in two:
    each fragment's mass DP is 3 - e / 2a < 0, so whole mu1 < 0, and Mn =
    (6a + e) / 2 is exported."""
    s = make_setup(tmp_path, {6: 1}, edit=_octyl_capped)
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    mu = declared(s["spec"])

    def split(t):
        return mu[0] * (1 + t / 10.0), mu[1], mu[2] * (1 - 0.5 * t / 10.0)

    octyl = "[CH2]CCCCCCC"
    write_run(out, times=[0.0, 10.0], pools=[(LABEL, split, 0.0)], gas=ZERO_GAS, end_groups=(octyl, octyl))
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    a = s["spec"].repeat_unit.composition.molar_mass_g_mol()
    e = s["spec"].end_composition.molar_mass_g_mol()
    assert q["whole_mu1_mol"] == pytest.approx(mu[0] * (6 - e / a), rel=1e-12) and q["whole_mu1_mol"] < 0
    assert q["whole_Mn_g_mol"] == pytest.approx((6 * a + e) / 2, rel=1e-12)


@pytest.mark.parametrize("keep_mu1", [True, False], ids=["0_6_36", "0_0_36"])
def test_a_pool_with_no_chains_and_units_left_is_corrupt(tmp_path, keep_mu1):
    """Spar r136 item 4 / r137 item 4: the pool ends at (0, 6, 36) u, or at
    (0, 0, 36) u whose variance is zero; chains without a count are
    corruption."""
    s = make_setup(tmp_path, {6: 1})
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    mu = declared(s["spec"])
    end = (0.0, mu[1] if keep_mu1 else 0.0, mu[2])
    write_run(out, times=[0.0, 10.0], pools=[(LABEL, lambda t: mu if t == 0.0 else end, 0.0)], gas=TRACE)
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    for key in ("polymer_condensed_mass_g", "polymer_conversion", *WHOLE):
        assert key not in q, key
    assert "units without chains" in json.loads((out / "mom_diagnostics.json").read_text())[-1]["polymer_withheld"]


def test_a_ledger_just_above_its_tolerance_withholds_mu2_and_mw(setup):
    """Survivor ledger-tol-loose: 10x the ledger tolerance (sum of e_p x 2f)."""
    out = setup["tmp"] / "out"
    produce(setup, out, end=10.0, snapshots=2)
    mu = declared(setup["spec"])
    row = (mu[0] - 20 * FLOOR, mu[1], mu[2])
    write_run(out, times=[0.0, 10.0], pools=[(LABEL, lambda t: mu if t == 0.0 else row, 0.0)], gas=TRACE)
    produce(setup, out, end=10.0, snapshots=2, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    assert "whole_mu2_mol" not in q and "whole_Mw_g_mol" not in q and "whole_Mn_g_mol" in q


@pytest.mark.parametrize("depth, clamped", [(0.5, True), (10.0, False)])
def test_a_negative_pool_moment_is_clamped_only_within_the_floor(setup, depth, clamped):
    """Survivor floor-x1000: a daughter pool at (0, 0, -depth f)."""
    out = setup["tmp"] / "out"
    produce(setup, out, end=10.0, snapshots=2)
    mu = declared(setup["spec"])
    write_run(out, times=[0.0, 10.0],
              pools=[(LABEL, lambda t: mu, 0.0),
                     (DAUGHTER, lambda t: (0.0, 0.0, 0.0) if t == 0.0 else (0.0, 0.0, -depth * FLOOR), 0.0)],
              gas=TRACE)
    produce(setup, out, end=10.0, snapshots=2, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    assert (set(WHOLE) <= set(q)) is clamped and ("polymer_condensed_mass_g" in q) is clamped


def test_an_isomer_of_a_pool_chain_is_not_a_candidate(tmp_path):
    """A declared isobutane has the formula of a capped PE DP-2 chain
    (n-butane) but not its structure: the pool's chains stay polymer."""
    isobutane = {"graph": {"atoms": [S._atom("c", hydrogens=1)] + [S._atom(f"m{k}", hydrogens=3) for k in range(3)],
                           "bonds": [{"atoms": ["c", f"m{k}"], "order": 1} for k in range(3)],
                           "electronic_state": {"multiplicity": 1}}, "disposition": "mobile"}
    s = make_setup(tmp_path, BASE, small={"isobutane": isobutane})
    out = tmp_path / "out"
    produce(s, out)
    write_run(out, times=PROFILE_T, pools=[(LABEL, drift(declared(s["spec"])), 0.0)],
              gas=ZERO_GAS + [("C4H10(4)", Molecule(smiles="CC(C)C").to_adjacency_list(), lambda t: 0.0)])
    produce(s, out, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    assert "polymer_derived_molecules_mol" in q and "isobutane" in q["species_mol"]
    assert json.loads((out / "producer.json").read_text())["read"]["pool_species_candidates"] == {LABEL: []}


# --------------------------------------------------------------------------
# I-170 rework 3 (r137)
# --------------------------------------------------------------------------
def _declared_graph(smiles):
    """A small-species declaration of a closed-shell hydrocarbon from SMILES."""
    mol = Molecule(smiles=smiles)
    heavy = [x for x in mol.atoms if not x.is_hydrogen()]
    key = {x: f"c{k}" for k, x in enumerate(heavy)}
    atoms = [S._atom(key[x], hydrogens=sum(1 for y in x.bonds if y.is_hydrogen())) for x in heavy]
    bonds, seen = [], set()
    for x in heavy:
        for y, bond in x.bonds.items():
            if y.is_hydrogen() or (y, x) in seen:
                continue
            seen.add((x, y))
            bonds.append({"atoms": [key[x], key[y]], "order": int(bond.order)})
    return {"graph": {"atoms": atoms, "bonds": bonds, "electronic_state": {"multiplicity": 1}},
            "disposition": "mobile"}


def _pool_run(tmp_path, small, pools_at_10, *, gas, monomer=None, main_empty=False):
    """Deck: two DP-10 chains.  At 10 s the main pool has given chains to the
    other pools (``pools_at_10``: (label, (mu0, mu1, mu2) at 10 s or a function
    of the deck's chain count giving them, defect)); it keeps its other chains
    at DP 10, or none (``main_empty``: the other pools then take every chain,
    so the end-group ledger stays zero).  The main
    pool keeps H caps, so every n-alkane declared is a candidate of it: a
    species row it may hold stands only when it is provably empty."""
    kw = {"monomer": monomer} if monomer is not None else {}
    s = make_setup(tmp_path, {10: 2}, small=small, **kw)
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    mu = declared(s["spec"])
    pools_at_10 = [(label, m(mu[0]) if callable(m) else m, d) for label, m, d in pools_at_10]
    gone = sum(m[0] for _l, m, _d in pools_at_10)
    assert not main_empty or gone == pytest.approx(mu[0], rel=1e-15)
    left = 0.0 if main_empty else mu[0] - gone
    pools = [(LABEL, lambda t: mu if t == 0.0 else (left, 10 * left, 100 * left), 0.0)]
    for label, m, d in pools_at_10:
        pools.append((label, (lambda m: lambda t: (0.0, 0.0, 0.0) if t == 0.0 else m)(m), d))
    write_run(out, times=[0.0, 10.0], pools=pools, gas=TRACE + list(gas))
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    return s, out, mom_entries(out)[-1]["quantities"], json.loads((out / "mom_diagnostics.json").read_text())[-1]


BUTANE_GAS = [("C4H10(4)", Molecule(smiles="CCCC").to_adjacency_list(), lambda t: 0.0)]
Y = 1e-3


def test_a_pinned_pool_with_a_known_structure_is_booked_as_the_species(tmp_path):
    """A pool without defect pinned at DP 2 (n-butane, declared): booked."""
    s, out, q, d = _pool_run(tmp_path, {"butane": S.BUTANE}, [(DAUGHTER, lambda n: (n, 2 * n, 4 * n), 0.0)],
                                gas=BUTANE_GAS, main_empty=True)
    assert q["species_mol"]["butane"] == declared(s["spec"])[0]
    assert q["polymer_derived_molecules_mol"] == 0.0
    assert d["pools_booked_as_species"] == {DAUGHTER: ["butane", 2]}


@pytest.mark.parametrize("moments", [lambda n: (n, 2 * n, 5 * n), lambda n: (n, 2 * n + 300 * FLOOR, 4 * n + 1200 * FLOOR)],
                         ids=["mean_2_polydisperse", "off_by_100_floors"])
def test_a_pool_not_pinned_at_the_species_dp_is_not_booked(tmp_path, moments):
    """Survivors pinned-mu2-ignored and pinned-loose: a pool whose mean is 2
    but whose mu2 is not 4 mu0, and one a hundred floors off DP 2."""
    _s, out, q, d = _pool_run(tmp_path, {"butane": S.BUTANE}, [(DAUGHTER, moments, 0.0)], gas=BUTANE_GAS,
                                main_empty=True)
    assert "butane" not in q["species_mol"] and "polymer_derived_molecules_mol" not in q
    assert d["pools_booked_as_species"] == {}


def test_a_mass_a_tenth_of_a_dp_off_is_not_a_candidate(setup):
    """Survivor mass-tol-loose: isobutane (C4H10) is 0.07 DP off the H2-loss
    pool's DP-2 chain mass (C4H8); not a candidate, the rows stay."""
    (setup["tmp"] / "iso").mkdir()
    s = make_setup(setup["tmp"] / "iso", BASE, small={"isobutane": _declared_graph("CC(C)C")})
    out = setup["tmp"] / "iso" / "out"
    produce(s, out)
    mu = declared(s["spec"])
    m_h2 = mw("[H][H]")
    y = 1e-3

    def main(t):
        return mu[0] - y * t * t, mu[1] - 15 * y * t * t, mu[2] - 250 * y * t * t

    def feature(t):
        return y * t * t, 15 * y * t * t, 250 * y * t * t

    write_run(out, times=PROFILE_T, pools=[(LABEL, main, 0.0), (LABEL + "_mod", feature, m_h2)],
              gas=[("H(1)", ADJ["H"], lambda t: 0.0), ("H2(2)", ADJ["H2"], lambda t: y * t * t),
                   ("C4H10(4)", Molecule(smiles="CC(C)C").to_adjacency_list(), lambda t: 0.0)])
    manifest = produce(s, out, mode="convert")
    assert manifest["read"]["pool_species_candidates"][LABEL + "_mod"] == []
    q = mom_entries(out)[-1]["quantities"]
    assert "polymer_derived_molecules_mol" in q and "isobutane" in q["species_mol"]


def test_a_pinned_defect_pool_is_never_booked_by_mass(tmp_path):
    """Spar r137 item 1: an H2-loss DP-2 pool (C4H8) with 1-butene declared:
    the pool may hold 2-butene; withheld, not booked as 1-butene."""
    butene = [("C4H8(4)", Molecule(smiles="C=CCC").to_adjacency_list(), lambda t: 0.0)]
    _s, out, q, d = _pool_run(tmp_path, {"butene": _declared_graph("C=CCC")},
                              [(LABEL + "_mod", (Y, 2 * Y, 4 * Y), mw("[H][H]"))], gas=butene)
    assert "butene" not in q["species_mol"] and "polymer_derived_molecules_mol" not in q
    assert d["pools_booked_as_species"] == {} and "structure" in d["polymer_withheld"]


def test_a_species_row_needs_every_contributing_pool_booked_or_empty(tmp_path):
    """Spar r137 item 2: one pool booked as butane next to another that may
    hold butane unpinned: the butane row is withheld, never a partial count."""
    pools = [(DAUGHTER, lambda n: (Y, 2 * Y, 4 * Y), 0.0),
             (DAUGHTER + "_b", lambda n: (n - Y, 2.5 * (n - Y), 6.5 * (n - Y)), 0.0)]
    _s, out, q, d = _pool_run(tmp_path, {"butane": S.BUTANE}, pools, gas=BUTANE_GAS, main_empty=True)
    assert "butane" not in q["species_mol"] and "H" in q["species_mol"]


def test_a_species_row_stands_when_the_blocked_pool_cannot_hold_it(tmp_path):
    """Spar r137 item 2: the butane pool is booked, the main pool is empty;
    an H2-loss pool (formula of ethylene, the
    declared monomer) is blocked but cannot hold butane: the butane row is
    exported in full, the monomer row and the polymer rows are withheld."""
    gas = BUTANE_GAS + [("C2H4(5)", Molecule(smiles="C=C").to_adjacency_list(), lambda t: 0.0)]
    pools = [(DAUGHTER, lambda n: (Y, 2 * Y, 4 * Y), 0.0), (LABEL + "_mod", lambda n: (n - Y,) * 3, mw("[H][H]"))]
    _s, out, q, d = _pool_run(tmp_path, {"butane": S.BUTANE}, pools, gas=gas, monomer=S.ETHYLENE_GRAPH,
                              main_empty=True)
    assert q["species_mol"]["butane"] == pytest.approx(Y, rel=1e-12)
    assert "monomer_product" not in q["species_mol"] and "polymer_derived_molecules_mol" not in q


def test_pool_chains_are_classified_at_t0_too(tmp_path, monkeypatch):
    """Spar r137 item 3: PE DP1 + DP6 with ethane declared.  kMC refuses the
    declaration (a declared species may not be an initial chain), and so does
    the producer; with that refusal bypassed, the t = 0 row withholds the
    rows the DP-1 chains would otherwise put in the polymer class."""
    ethane = _declared_graph("CC")
    s = make_setup(tmp_path, {1: 1, 6: 1}, small={"ethane": ethane})
    out = tmp_path / "out"
    produce(s, out, end=10.0, snapshots=2)
    mu = declared(s["spec"])
    gas = TRACE + [("C2H6(4)", Molecule(smiles="CC").to_adjacency_list(), lambda t: 0.0)]
    write_run(out, times=[0.0, 10.0], pools=[(LABEL, drift(mu), 0.0)], gas=gas)
    with pytest.raises(mp.ProducerError, match="initial chain"):
        produce(s, out, end=10.0, snapshots=2, mode="convert")
    monkeypatch.setattr(mo, "refuse_initial_species_overlap", lambda *args, **kwargs: None)
    produce(s, out, end=10.0, snapshots=2, mode="convert")
    q0 = mom_entries(out)[0]["quantities"]
    for key in ("polymer_condensed_mass_g", "polymer_derived_molecules_mol", *WHOLE):
        assert key not in q0, key
    assert "ethane" not in q0["species_mol"]


def test_a_trace_population_is_not_corruption(setup):
    """Spar r137 item 4: a monodisperse DP-30 trace (5e-15, 1.5e-13, 4.5e-12)
    mol, mu0 under the floor, is a valid state; (0, 6, 36) is not."""
    out = setup["tmp"] / "out"
    produce(setup, out, end=10.0, snapshots=2)
    mu = declared(setup["spec"])
    trace = (5e-15, 1.5e-13, 4.5e-12)
    write_run(out, times=[0.0, 10.0],
              pools=[(LABEL, lambda t: mu, 0.0), (DAUGHTER, lambda t: (0.0, 0.0, 0.0) if t == 0.0 else trace, 0.0)],
              gas=TRACE)
    produce(setup, out, end=10.0, snapshots=2, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    assert {"polymer_condensed_mass_g", "polymer_derived_molecules_mol", *WHOLE} <= set(q)


def test_an_h2_loss_transfer_does_not_make_a_pool_light(setup):
    """r137 item 6: a volatile ejection that moves H2 into the defect pool has
    eject_units = M(H2) / a but moment units a_mass - (d_dst - d_src) / a_src
    = 0 (polymer.pyx:4110-4113): the DP is unchanged, the pool is not light."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    m_h2 = mw("[H][H]")
    a, _e = masses()
    y = 1e-3
    proxy = proxy_of()

    def main(t):
        return mu[0] - y * t * t, mu[1] - 15 * y * t * t, mu[2] - 250 * y * t * t

    def feature(t):
        return y * t * t, 15 * y * t * t, 250 * y * t * t

    rows = [{"family": "toy_h_abstraction", "admitted": True, "reactants": [proxy], "products": [proxy, "X0"],
             "reversible": False, "archetype": "volatile_ejection/1",
             "extra": {"params": {"eject_units": m_h2 / a}, "src_pool": LABEL, "dst_pool": LABEL + "_mod"}}]
    write_run(out, times=PROFILE_T, pools=[(LABEL, main, 0.0), (LABEL + "_mod", feature, m_h2)],
              gas=[("H(1)", ADJ["H"], lambda t: 0.0), ("H2(2)", ADJ["H2"], lambda t: y * t * t)], rows=rows)
    produce(setup, out, mode="convert")
    q = mom_entries(out)[-1]["quantities"]
    assert "polymer_derived_molecules_mol" in q and "polymer_condensed_mass_g" in q


# --------------------------------------------------------------------------
# I-172: the activity gate's rates by RMG, physical motion, colliders
# --------------------------------------------------------------------------
KCAL = 4184.0


def _row(**kw):
    """One admitted homolysis row of the main pool (proxy -> proxy + X0)."""
    proxy = proxy_of()
    return {"family": "toy_homolysis", "admitted": True, "reactants": [proxy], "products": [proxy, "X0"],
            "reversible": False, **kw}


def _census(setup, rows, **kw):
    """(census, refused?) of a moving run with ``rows``."""
    out = setup["tmp"] / f"out{len(list(setup['tmp'].glob('out*')))}"
    produce(setup, out)
    constant_run(out, setup["spec"], rows=rows, **kw)
    try:
        manifest = produce(setup, out, mode="convert")
    except mp.ProducerError as error:
        assert "inactive_chemistry" in str(error), error
        return json.loads((out / "producer.json").read_text())["census"], True
    return manifest["census"], False


def test_a_plog_whose_pressure_rates_are_all_zero_is_dead(setup):
    """Review defect 1: RMG's writer puts a placeholder 1 0 0 on the main
    line of a PLOG; the rate is on the PLOG lines."""
    plog = ["PLOG /  0.1 0.0 0.0 0.0 /", "PLOG /  10.0 0.0 0.0 0.0 /"]
    census, refused = _census(setup, [_row(rate=(1.0, 0.0, 0.0), aux=plog, sidecar_kinetics=None)])
    assert refused and census["zero_rate_rows"] == ["r0"]


def test_a_plog_with_positive_pressure_rates_and_a_zero_main_line_is_live(setup):
    plog = ["PLOG /  0.1 2.0 0.0 0.0 /", "PLOG /  10.0 3.0 0.0 0.0 /"]
    census, refused = _census(setup, [_row(rate=(0.0, 0.0, 0.0), aux=plog, sidecar_kinetics=None)])
    assert not refused and census["live_rows"] == ["r0"]


def test_a_lindemann_falloff_with_a_zero_low_pressure_rate_is_dead(setup):
    proxy = proxy_of()
    row = _row(equation=f"{proxy}(+M)=>{proxy}+X0(+M)", rate=(1.0, 0.0, 0.0), aux=["LOW /  0.0 0.0 0.0 /"],
               sidecar_kinetics=None)
    census, refused = _census(setup, [row])
    assert refused and census["zero_rate_rows"] == ["r0"]


@pytest.mark.parametrize("a", [2.0, 1.0])
def test_cancelling_multi_arrhenius_terms_are_dead(setup, a):
    """Review defect 2: terms (a, 0, 0) and (-a, 0, 0) sum to zero.  (RMG's
    Arrhenius cannot return exactly -1: its except -1 sentinel raises, a
    stated rate error, dead as well.)"""
    census, refused = _census(setup, [_row(expanded=True, terms=[(a, 0.0, 0.0), (-a, 0.0, 0.0)],
                                           sidecar_kinetics=None)])
    assert refused and census["zero_rate_rows"] == ["r0"]


def test_a_multi_arrhenius_with_one_zero_term_is_live(setup):
    census, refused = _census(setup, [_row(expanded=True, terms=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
                                           sidecar_kinetics=None)])
    assert not refused and census["live_rows"] == ["r0"]


OVERFLOW = {"rate": (0.0, 0.0, -1e7 / KCAL), "sidecar_kinetics": {"A": 0.0, "n": 0.0, "Ea": -1e7}}


def test_a_refused_row_whose_rate_overflows_is_never_evaluated(setup):
    """Review defect 3: A = 0, Ea = -1e7 J/mol at 700 K on a refused row:
    the admission filter comes first, so the census converts and records no
    rate error; on an admitted row the error is a stated reason, not a
    crash."""
    proxy = proxy_of()
    refused_row = _row(admitted=False, products=[proxy, "X1"], **OVERFLOW)
    census, refused = _census(setup, [refused_row, _row()])
    assert not refused and census["live_rows"] == ["r1"] and census["rate_errors"] == {}
    census, refused = _census(setup, [_row(**OVERFLOW)])
    assert refused and census["zero_rate_rows"] == ["r0"] and list(census["rate_errors"]) == ["r0"]


@pytest.mark.parametrize("units, dead", [("KCAL/MOLE", True), ("CAL/MOLE", False)])
def test_the_reactions_block_ea_unit_is_honoured(setup, units, dead):
    """Ea = 1e6 in the header's unit: kcal/mol underflows exp(-Ea / RT) to 0 at
    700 K, cal/mol leaves 2e-313."""
    census, refused = _census(setup, [_row(rate=(1.0, 0.0, 1e6), sidecar_kinetics=None)], ea_units=units)
    assert refused is dead and (census["zero_rate_rows"] == ["r0"]) is dead


def test_more_chemkin_reactions_than_sidecar_rows_is_refused(setup):
    proxy = proxy_of()
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], rows=[_row()],
                 gas_reactions=[{"reactants": [proxy], "products": [proxy, "X0"], "reversible": False, "A": 1.0}])
    with pytest.raises(mp.ProducerError, match="no one-to-one match"):
        produce(setup, out, mode="convert")


def test_an_indexed_specific_collider_matches_its_row(setup):
    """Review defect 5: P(3)(+H2(2)) => Q + H(1)(+H2(2))."""
    proxy = proxy_of()
    row = _row(equation=f"{proxy}(+H2(2))=>{proxy}+X0(+H2(2))", rate=(1.0, 0.0, 0.0),
               aux=["LOW /  1.0 0.0 0.0 /"], sidecar_kinetics=None, colliders=["H2(2)"])
    census, refused = _census(setup, [row])
    assert not refused and census["live_rows"] == ["r0"]


def _motion(setup, pool, *, gas=ZERO_GAS, **kw):
    out = setup["tmp"] / "out"
    produce(setup, out)
    write_run(out, times=PROFILE_T, pools=[(LABEL, pool, 0.0)], gas=gas, **kw)
    try:
        return produce(setup, out, mode="convert")["census"]["profile_changes"]
    except mp.ProducerError as error:
        assert "inactive_chemistry" in str(error), error
        return json.loads((out / "producer.json").read_text())["census"]["profile_changes"]


def test_a_real_moment_change_of_5e_7_moves_the_profile(setup):
    """Review defect 4: mu1 falls by 5e-7 of itself by the horizon."""
    mu = declared(setup["spec"])
    assert _motion(setup, lambda t: (mu[0], mu[1] * (1 - 5e-7 * t / 40.0), mu[2])) is True


def test_a_lone_1e_300_appearance_does_not_move_the_profile(setup):
    mu = declared(setup["spec"])
    gas = [("H(1)", ADJ["H"], lambda t: 0.0 if t == 0.0 else 1e-300), ("H2(2)", ADJ["H2"], lambda t: 0.0)]
    assert _motion(setup, lambda t: mu, gas=gas) is False


def test_motion_of_a_bookkeeping_proxy_alone_does_not_move_the_profile(setup):
    mu = declared(setup["spec"])
    assert _motion(setup, lambda t: mu, proxy_amount=lambda t: 7.0 + t) is False


def test_motion_after_the_horizon_does_not_move_the_profile(setup):
    """The solver's last row (45 s) passes the 40 s horizon; only it moves."""
    mu = declared(setup["spec"])
    assert _motion(setup, lambda t: mu if t <= 40.0 else (2 * mu[0], mu[1], mu[2])) is False


def test_a_species_identifier_holding_end_keeps_its_reactions(setup):
    """Spar r140 item 6: RMG's reader ended a REACTIONS block at any line
    holding 'end'; with END recognised as Chemkin's first-token keyword the
    row is read and live."""
    proxy = proxy_of()
    census, refused = _census(setup, [_row(products=[proxy, "Xend"])])
    assert not refused and census["live_rows"] == ["r0"]


def test_a_reaction_rmgs_reader_drops_is_refused(setup, monkeypatch):
    """The census' backstop: a reader that returns fewer reactions than the
    file annotates is refused, never counted short."""
    import rmgpy.chemkin as chemkin
    real = chemkin.load_chemkin_file

    def short(*args, **kwargs):
        species, reactions = real(*args, **kwargs)
        return species, reactions[:-1]

    monkeypatch.setattr(chemkin, "load_chemkin_file", short)
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], rows=[_row()])
    with pytest.raises(mp.ProducerError, match="RMG's Chemkin reader gives 0 reactions"):
        produce(setup, out, mode="convert")


# --------------------------------------------------------------------------
# I-172 rework 1 (spar r140)
# --------------------------------------------------------------------------
def test_a_reverse_rate_mom_disables_does_not_fire(setup):
    """Spar r140 item 1: X5 <=> proxy with X5 absent and the proxy present;
    the proxy's enthalpy makes Keq overflow, so MoM sets kb = 0
    (polymer.pyx 4751-4762): the row carries no flux while the declared gas
    moves.  Refused.  (With equal thermo the reverse runs:
    test_a_reversible_row_runs_from_its_products.)"""
    proxy = proxy_of()
    out = setup["tmp"] / "out"
    produce(setup, out)
    mu = declared(setup["spec"])
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: mu, 0.0)], gas=TRACE,
              rows=[_row(reactants=["X5"], products=[proxy], reversible=True)], enthalpy={proxy: -1.0e6})
    census = _refused_census(setup, out)
    assert census["live_rows"] == [] and census["profile_changes"] is True
    assert "r0" in census["reverse_disabled"]


def test_a_falloff_rate_takes_moms_effective_pressure(setup):
    """Spar r140 item 2: a Lindemann row whose only collider in the gas (the
    N2 bath, H and H2 at zero) has efficiency 0: MoM's effective pressure
    (polymer.pyx 6546-6562) is 0 and so is its rate."""
    proxy = proxy_of()
    row = _row(equation=f"{proxy}(+M)=>{proxy}+X0(+M)", rate=(1.0, 0.0, 0.0),
               aux=["LOW /  1.0 0.0 0.0 /", "N2/0.0/"], sidecar_kinetics=None)
    census, refused = _census(setup, [row])
    assert refused and census["zero_rate_rows"] == ["r0"]
    census, refused = _census(setup, [dict(row, aux=["LOW /  1.0 0.0 0.0 /", "N2/1.0/"])])
    assert not refused and census["live_rows"] == ["r0"]


def test_motion_is_judged_on_the_interpolated_observations(setup):
    """Spar r140 item 3: stored rows at 0 and 45 s, horizon 40 s: conversion
    interpolates mu0 at every observation; so does the gate."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    write_run(out, times=[0.0, 45.0], pools=[(LABEL, drift(declared(setup["spec"])), 0.0)], gas=ZERO_GAS)
    manifest = produce(setup, out, mode="convert")
    assert manifest["census"]["profile_changes"] is True


def test_a_change_of_mu2_alone_moves_the_profile(tmp_path):
    """Spar r140 item 4: DP2 + DP4 -> DP3 + DP3 keeps mu0, mu1 and every
    species but moves Mw."""
    s = make_setup(tmp_path, {2: 1, 4: 1})
    out = tmp_path / "out"
    produce(s, out)
    mu = declared(s["spec"])
    assert mu[2] == pytest.approx(20.0 / 6.0 * mu[1], rel=1e-12)
    write_run(out, times=PROFILE_T, pools=[(LABEL, lambda t: mu if t == 0.0 else (mu[0], mu[1], 9.0 * mu[0]), 0.0)],
              gas=ZERO_GAS)
    assert produce(s, out, mode="convert")["census"]["profile_changes"] is True


def test_a_negative_total_rate_is_a_rate_error(setup):
    """Spar r140 item 5: terms (2, 0, 0) and (-4, 0, 0) sum to -2."""
    census, refused = _census(setup, [_row(expanded=True, terms=[(2.0, 0.0, 0.0), (-4.0, 0.0, 0.0)],
                                           sidecar_kinetics=None)])
    assert refused and census["zero_rate_rows"] == ["r0"] and "negative" in census["rate_errors"]["r0"]


# --------------------------------------------------------------------------
# I-172 rework 2 (spar r141): MoM's gas set and effective pressure
# --------------------------------------------------------------------------
def _falloff_row():
    """A Lindemann row whose N2 collider efficiency is 0."""
    proxy = proxy_of()
    return _row(equation=f"{proxy}(+M)=>{proxy}+X0(+M)", rate=(1.0, 0.0, 0.0),
                aux=["LOW /  1.0 0.0 0.0 /", "N2/0.0/"], sidecar_kinetics=None)


def _falloff_census(setup, **kw):
    out = setup["tmp"] / f"out{len(list(setup['tmp'].glob('out*')))}"
    produce(setup, out)
    kw.setdefault("pools", [(LABEL, drift(declared(setup["spec"])), 0.0)])
    kw.setdefault("gas", ZERO_GAS)
    write_run(out, times=PROFILE_T, rows=[_falloff_row()], **kw)
    try:
        return produce(setup, out, mode="convert")["census"], False
    except mp.ProducerError as error:
        assert "inactive_chemistry" in str(error), error
        return json.loads((out / "producer.json").read_text())["census"], True


RADICAL = Molecule(smiles="[CH2]CCCCCCC").to_adjacency_list()


def test_a_condensed_species_is_no_collider(setup):
    """Spar r141 item 1: 0.1 mol of a daughter MoM condenses
    (conventions.condensed_species) beside the N2 bath of efficiency 0: MoM's
    effective pressure is 0, not 1e5 x 0.1 / 1.1 Pa."""
    census, refused = _falloff_census(setup, gas=ZERO_GAS + [("PD(7)", RADICAL, lambda t: 0.0 if t == 0.0 else 0.1)],
                                      condensed=["PD(7)"])
    assert refused and census["zero_rate_rows"] == ["r0"]


def test_compressed_bookkeeping_columns_are_no_colliders(setup):
    """Spar r141 item 2: a daughter pool with a long label: the profile names
    its mu columns compressed and its proxy by formula, the sidecar by full
    label; mapped by moment_column / species_key, they are not gas."""
    mu = declared(setup["spec"])
    full = [f"{DAUGHTER}_mu{k}" for k in range(3)] + [f"{DAUGHTER}(5)"]
    census, refused = _falloff_census(
        setup, pools=[(LABEL, drift(mu), 0.0), (DAUGHTER, lambda t: (0.0, 0.0, 0.0) if t == 0.0 else (0.1, 0.5, 4.0), 0.0)],
        gas=ZERO_GAS + [("C2H2(5)", RADICAL, lambda t: 0.0 if t == 0.0 else 0.3)], bookkeeping={DAUGHTER: full}, condensed=full)
    assert refused and census["zero_rate_rows"] == ["r0"]


def test_an_empty_gas_keeps_the_previous_rate(setup):
    """Spar r141 item 3: N2 (efficiency 0) at t = 0, then no gas at all: MoM
    skips the kf update on an empty gas (polymer.pyx 6553-6564) and keeps 0."""
    census, refused = _falloff_census(setup, n2=lambda t: 1.0 if t == 0.0 else 0.0)
    assert refused and census["zero_rate_rows"] == ["r0"]


def test_a_gas_empty_from_t0_keeps_the_initial_rate(setup):
    """The first state empty: MoM's kf stays at its initial
    get_rate_coefficient(T, P) at the total pressure (polymer.pyx 4744):
    live."""
    census, refused = _falloff_census(setup, condensed=["N2"])
    assert not refused and census["live_rows"] == ["r0"]


def test_a_profile_column_that_maps_to_nothing_is_refused(setup):
    """Spar r141 item 2: every column is a pool moment or a species of the
    Chemkin file; one that is neither cannot be put in or out of the gas."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], extra_columns=("Q(9)",))
    with pytest.raises(mp.ProducerError, match="neither a pool moment nor a species"):
        produce(setup, out, mode="convert")


# --------------------------------------------------------------------------
# I-172 rework 3 (spar r142): MoM's collider efficiencies and gas mask
# --------------------------------------------------------------------------
HEXANE = Molecule(smiles="CCCCCC").to_adjacency_list()


def test_collider_efficiencies_are_assigned_over_the_full_core(setup):
    """Spar r142 item 1: the efficiency line gives hexane 0; MoM assigns it to
    the FIRST isomorphic core species (model.pyx 415-421 over the core list,
    polymer_input.py 428-434), the condensed hexane-shaped proxy, so the gas
    hexane keeps 1: Peff = 1e5 x 1 / (1 + 1) Pa and the row is live.  (MoM's
    core order is the profile's: the proxy column comes first.)"""
    proxy = proxy_of()
    row = _row(equation=f"{proxy}(+M)=>{proxy}+X0(+M)", rate=(1.0, 0.0, 0.0),
               aux=["LOW /  1.0 0.0 0.0 /", "N2/0.0/ C6H14(4)/0.0/"], sidecar_kinetics=None)
    out = setup["tmp"] / "out"
    produce(setup, out)
    write_run(out, times=PROFILE_T, pools=[(LABEL, drift(declared(setup["spec"])), 0.0)], rows=[row],
              gas=ZERO_GAS + [("C6H14(4)", HEXANE, lambda t: 0.0 if t == 0.0 else 1.0)], proxy_first=True)
    assert produce(setup, out, mode="convert")["census"]["live_rows"] == ["r0"]


def test_a_spawned_pools_phase_is_an_ambiguous_gas_mask(setup):
    """Spar r142 item 2: an unconfigured spawned daughter's phase species are
    in conventions.condensed_species (the serializer's 2.5 closure) though
    the solver keeps them gas unless it configured the pool, which the
    sidecar does not record.  In the gas the row is live, out of it dead:
    refused with the stated ambiguity."""
    mu = declared(setup["spec"])
    full = [f"{DAUGHTER}_mu{k}" for k in range(3)] + [f"{DAUGHTER}(5)"]
    out = setup["tmp"] / "out"
    produce(setup, out)
    write_run(out, times=PROFILE_T, rows=[_falloff_row()],
              pools=[(LABEL, drift(mu), 0.0), (DAUGHTER, lambda t: (0.0, 0.0, 0.0) if t == 0.0 else (0.1, 0.5, 4.0), 0.0)],
              gas=ZERO_GAS + [("C2H2(5)", RADICAL, lambda t: 0.0 if t == 0.0 else 0.3)],
              bookkeeping={DAUGHTER: full}, condensed=full, spawned=[DAUGHTER])
    census = _refused_census(setup, out)
    assert census["zero_rate_rows"] == ["r0"] and "gas mask ambiguous" in census["rate_errors"]["r0"]


# --------------------------------------------------------------------------
# I-176: MoM's core order, reverse-rate ambiguity, the ambiguity cap
# --------------------------------------------------------------------------
def test_collider_efficiencies_follow_moms_core_order_not_the_species_block(setup):
    """I-176 item 1: MoM's core [N2, P(10), C6H14(4)] (the profile's order;
    RMG's Chemkin writer sorts the SPECIES block by index, C6H14(4) first).
    The hexane efficiency 0 goes to the first isomorphic core species, the
    condensed hexane-shaped proxy P(10): the gas hexane keeps 1 and
    Peff = 5e4 Pa, live."""
    proxy = f"{LABEL}(10)"
    row = _row(reactants=[proxy], products=[proxy, "X0"], equation=f"{proxy}(+M)=>{proxy}+X0(+M)",
               rate=(1.0, 0.0, 0.0), aux=["LOW /  1.0 0.0 0.0 /", "N2/0.0/ C6H14(4)/0.0/"], sidecar_kinetics=None)
    out = setup["tmp"] / "out"
    produce(setup, out)
    write_run(out, times=PROFILE_T, pools=[(LABEL, drift(declared(setup["spec"])), 0.0)], rows=[row],
              gas=ZERO_GAS + [("C6H14(4)", HEXANE, lambda t: 0.0 if t == 0.0 else 1.0)], proxy_index=10,
              proxy_first=True)
    assert produce(setup, out, mode="convert")["census"]["live_rows"] == ["r0"]


class _Kinetics:
    """A falloff-like kinetics whose kf is Peff / P (1 at P)."""
    efficiencies = {"collider": 2.0}

    def is_pressure_dependent(self):
        return True

    def get_effective_collider_efficiencies(self, species):
        return [1.0, 2.0]

    def get_rate_coefficient(self, T, P):
        return P / 1.0e5


class _Reaction:
    kinetics = _Kinetics()

    def __init__(self, keq):
        self.keq = keq

    def get_equilibrium_constant(self, T):
        return self.keq


def test_a_reverse_rate_the_gas_masks_disagree_on_is_dead_with_the_reason():
    """I-176 item 2: N2 (efficiency 1) and a spawned pool's species
    (efficiency 2), one mol each: out of the gas kf = 1, in it kf = 1.5; with
    Keq = 7e-309 kb = 1.4e308 is finite under one mask and overflows (MoM: kb
    = 0) under the other.  Forward live; reverse dead, the reason stated."""
    core = [(0, "N2", (None, "gas")), (1, "P2(5)", ("configured:P2", "gas"))]   # P2's first member
    rates = mp.mom_rates([_Reaction(7e-309)], True, 700.0, 1.0e5, core, [[1.0, 1.0]])
    assert rates["forward"] is True and rates["reverse"] is False
    assert "gas mask ambiguous" in rates["error"] and "reverse" in rates["error"]
    rates = mp.mom_rates([_Reaction(1.0)], True, 700.0, 1.0e5, core, [[1.0, 1.0]])
    assert rates["forward"] is True and rates["reverse"] is True and rates["error"] is None


def test_the_ambiguity_cap_counts_spawned_pools_with_a_present_species(setup):
    """I-176 item 3: three spawned pools (12 bookkeeping columns, their mu
    dummies and proxies): every efficiency is 1, so every mask gives
    Peff = P; the ambiguity is three pools, under the cap, and the row is
    live."""
    mu = declared(setup["spec"])
    proxy = proxy_of()
    row = _row(equation=f"{proxy}(+M)=>{proxy}+X0(+M)", rate=(1.0, 0.0, 0.0), aux=["LOW /  1.0 0.0 0.0 /"],
               sidecar_kinetics=None)
    row["aux"] = ["LOW /  1.0 0.0 0.0 /", "X0/1.0/"]             # an efficiency line: MoM's effective pressure
    daughters = [f"{LABEL}_d{k}_rad_primary_end" for k in range(3)]
    pools = [(LABEL, drift(mu), 0.0)] + [(d, lambda t: (0.0, 0.0, 0.0) if t == 0.0 else (0.1, 0.5, 4.0), 0.0)
                                         for d in daughters]
    gas = ZERO_GAS + [(f"C2H2({5 + k})", RADICAL, lambda t: 0.0 if t == 0.0 else 0.3) for k in range(3)]
    full = {d: [f"{d}_mu{j}" for j in range(3)] + [f"{d}({5 + k})"] for k, d in enumerate(daughters)}
    out = setup["tmp"] / "out"
    produce(setup, out)
    write_run(out, times=PROFILE_T, rows=[row], pools=pools, gas=gas, bookkeeping=full,
              condensed=[x for v in full.values() for x in v], spawned=daughters)
    manifest = produce(setup, out, mode="convert")
    assert manifest["census"]["live_rows"] == ["r0"] and manifest["census"]["rate_errors"] == {}


# --------------------------------------------------------------------------
# I-176 rework 1 (spar r147): the solver's own masks for a spawned pool
# --------------------------------------------------------------------------
PROPANE = Molecule(smiles="CCC").to_adjacency_list()
BUTANE_ADJ = Molecule(smiles="CCCC").to_adjacency_list()


def _first_member_run(setup, gas_first):
    """A spawned pool ``dgt`` whose members are its proxy ``dgt`` (propane-
    shaped) and a gas species ``dgt(4)`` (butane), both listed by the sidecar
    (proxy first); ``gas_first``: dgt(4) precedes dgt in MoM's core.  N2 and
    dgt(4) have collider efficiency 0, dgt 1.  A configured pool condenses
    only its FIRST core species with the pool's base label
    (_apply_pool_phase_overrides, polymer.pyx 3204-3212); every other member
    keeps the first stage's phase (gas here)."""
    proxy = proxy_of()
    row = _row(equation=f"{proxy}(+M)=>{proxy}+X0(+M)", rate=(1.0, 0.0, 0.0),
               aux=["LOW /  1.0 0.0 0.0 /", "N2/0.0/ dgt(4)/0.0/"], sidecar_kinetics=None)
    members = [("dgt(4)", BUTANE_ADJ, lambda t: 0.0 if t == 0.0 else 1.0),
               ("dgt", PROPANE, lambda t: 0.0 if t == 0.0 else 1.0)]
    if not gas_first:
        members.reverse()
    sidecar = ["dgt_mu0", "dgt_mu1", "dgt_mu2", "dgt", "dgt(4)"]
    out = setup["tmp"] / f"out{len(list(setup['tmp'].glob('out*')))}"
    produce(setup, out)
    write_run(out, times=PROFILE_T, rows=[row], gas=ZERO_GAS + members,
              pools=[(LABEL, drift(declared(setup["spec"])), 0.0), ("dgt", lambda t: (0.0, 0.0, 0.0), 0.0)],
              bookkeeping={"dgt": sidecar}, condensed=sidecar, spawned=["dgt"])
    try:
        return produce(setup, out, mode="convert")["census"], False
    except mp.ProducerError as error:
        assert "inactive_chemistry" in str(error), error
        return json.loads((out / "producer.json").read_text())["census"], True


def test_a_configured_pool_condenses_only_its_first_member_in_core_order(setup):
    """Spar r147 item 1, the gas species first: configured, dgt(4) is
    condensed and dgt (efficiency 1) is gas; not configured, both are gas:
    live either way, so the row is live (grouping the members would add the
    all-condensed mask, Peff = 0)."""
    census, refused = _first_member_run(setup, gas_first=True)
    assert not refused and census["live_rows"] == ["r0"] and census["rate_errors"] == {}


def test_the_proxy_first_gives_the_solvers_other_verdict(setup):
    """The proxy first: configured, dgt is condensed and only efficiency-0
    colliders hold the gas (Peff = 0, dead); not configured, live: the masks
    disagree and the run is refused with the reason."""
    census, refused = _first_member_run(setup, gas_first=False)
    assert refused and census["zero_rate_rows"] == ["r0"]
    assert "gas mask ambiguous" in census["rate_errors"]["r0"]


def test_a_dead_verdict_states_its_sampled_compositions(setup):
    """Spar r147 item 2: the saved rows are samples of the solver's states;
    a dead verdict says so."""
    census, refused = _falloff_census(setup)
    assert refused and census["zero_rate_rows"] == ["r0"]
    n = len([t for t in PROFILE_T if t < 40.0]) + 1
    assert f"at the {n} sampled compositions" in census["dead_reasons"]["r0"]


# --------------------------------------------------------------------------
# I-176 rework 2 (spar r151): the solver's recorded mask; the reconstruction
# --------------------------------------------------------------------------
HEXYL = Molecule(smiles="[CH2]CCCCC").to_adjacency_list()      # an H-loss daughter of the hexane-shaped proxy


def _phase_run(setup, *, members, aux, bookkeeping, spawned=("dgt",), pools=None, **kw):
    """The falloff row ``proxy(+M) => proxy + X0(+M)`` (A = 1) over MoM's
    core: N2, the main pool's moments, a pool ``dgt`` (no amounts), the main
    proxy (7 mol, condensed) and ``members`` [(column, adjacency, f(t))];
    ``aux``: the efficiency line.  Returns (census, refused)."""
    proxy = proxy_of()
    row = _row(equation=f"{proxy}(+M)=>{proxy}+X0(+M)", rate=(1.0, 0.0, 0.0),
               aux=["LOW /  1.0 0.0 0.0 /", aux], sidecar_kinetics=None)
    out = setup["tmp"] / f"out{len(list(setup['tmp'].glob('out*')))}"
    produce(setup, out)
    if pools is None:
        pools = [(LABEL, drift(declared(setup["spec"])), 0.0), ("dgt", lambda t: (0.0, 0.0, 0.0), 0.0)]
    write_run(out, times=PROFILE_T, rows=[row], gas=ZERO_GAS + list(members), pools=pools,
              bookkeeping=bookkeeping, spawned=list(spawned), **kw)
    try:
        return produce(setup, out, mode="convert")["census"], False
    except mp.ProducerError as error:
        assert "inactive_chemistry" in str(error), error
        return json.loads((out / "producer.json").read_text())["census"], True


DGT = ["dgt_mu0", "dgt_mu1", "dgt_mu2", "dgt"]
LATER = lambda t: 0.0 if t == 0.0 else 1.0      # noqa: E731


def test_section_d_condenses_every_species_of_a_qualifying_base(setup):
    """Spar r151 P1: dgt(4) (hexyl) is an H-loss daughter of the root proxy;
    MoM's first stage then condenses EVERY core species of base ``dgt``
    (get_gas_mask section D), the butane dgt(5) too.  dgt(5) holds the only
    collider weight: in the gas the row is live, out of it dead.  Whether
    dgt(4) qualifies the sidecar cannot settle (gas veto, Polymer-ness), so
    the verdict is the ambiguity, its unit dgt(4)'s."""
    census, refused = _phase_run(
        setup, members=[("dgt", PROPANE, lambda t: 0.0), ("dgt(4)", HEXYL, lambda t: 0.0), ("dgt(5)", BUTANE_ADJ, LATER)],
        aux="N2/0.0/ dgt/0.0/ dgt(4)/0.0/", bookkeeping={"dgt": DGT + ["dgt(4)", "dgt(5)"]}, condensed=DGT)
    assert census["phase_source"] == "reconstructed"
    assert refused and census["zero_rate_rows"] == ["r0"]
    assert "gas mask ambiguous for ['stage1:dgt(4)']" in census["rate_errors"]["r0"]


def test_a_stage1_daughter_is_an_ambiguity_unit(setup):
    """The H-loss daughter dgt(4) itself, a member the configured override
    does not reach (not the pool's first): its first-stage phase is the
    unit stage1:dgt(4)."""
    census, refused = _phase_run(
        setup, members=[("dgt", PROPANE, lambda t: 0.0), ("dgt(4)", HEXYL, LATER)],
        aux="N2/0.0/ dgt/0.0/", bookkeeping={"dgt": DGT + ["dgt(4)"]}, condensed=DGT)
    assert refused and "gas mask ambiguous for ['stage1:dgt(4)']" in census["rate_errors"]["r0"]


@pytest.mark.parametrize("order", [(LABEL, "dgt"), ("dgt", LABEL)])
def test_a_member_of_two_pools_is_condensed_whatever_their_order(setup, order):
    """Spar r151 P1: X(6) is a member of the root pool and of the spawned
    pool dgt.  MoM's overrides only condense, so the root pool's
    condensation stands whatever the order: X(6), the only collider weight,
    is out of the gas and the row dead, with no ambiguity."""
    main = [mp.moment_column(f"{LABEL}_mu{k}") for k in range(3)] + [proxy_of()]
    census, refused = _phase_run(
        setup, members=[("dgt", PROPANE, lambda t: 0.0), ("X(6)", BUTANE_ADJ, LATER)], aux="N2/0.0/ dgt/0.0/",
        bookkeeping={LABEL: main + ["X(6)"], "dgt": DGT + ["X(6)"]}, condensed=DGT + ["X(6)"], pool_order=order)
    assert refused and census["zero_rate_rows"] == ["r0"] and census["rate_errors"] == {}


def test_an_explicit_member_of_a_root_pool_is_condensed(setup):
    """A root pool's member of another base (an explicit-DP species, PD(7))
    is condensed by the pool's override: no collider weight is left in the
    gas, the row is dead, and no phase is in doubt."""
    main = [mp.moment_column(f"{LABEL}_mu{k}") for k in range(3)] + [proxy_of()]
    census, refused = _phase_run(
        setup, members=[("PD(7)", RADICAL, LATER)], aux="N2/0.0/", bookkeeping={LABEL: main + ["PD(7)"]},
        spawned=(), pools=[(LABEL, drift(declared(setup["spec"])), 0.0)])
    assert refused and census["zero_rate_rows"] == ["r0"] and census["rate_errors"] == {}


def test_a_root_pool_is_no_ambiguity(setup):
    """The main (root) pool is configured: its proxy is condensed, the row
    dead for want of gas colliders, and no phase is in doubt."""
    census, refused = _falloff_census(setup)
    assert refused and census["zero_rate_rows"] == ["r0"] and census["rate_errors"] == {}
    assert "ambiguous" not in census["dead_reasons"]["r0"]


class _ManyKinetics(_Kinetics):
    def __init__(self, n):
        self.n = n

    def get_effective_collider_efficiencies(self, species):
        return [1.0] * self.n


class _ManyReaction(_Reaction):
    def __init__(self, n):
        self.kinetics = _ManyKinetics(n)
        self.keq = 1.0


def test_condensed_rules_bring_no_units_to_the_cap():
    """Spar r151 P2: eleven species a configured pool condenses whatever
    their first stage, (True, stage1 unit): simplified, they carry no unit,
    so the cap is not reached and N2 alone keeps the row live."""
    core = [(0, "N2", (None, "gas"))] + [(1 + k, f"S({k})", (True, f"stage1:S({k})")) for k in range(11)]
    rates = mp.mom_rates([_ManyReaction(12)], False, 700.0, 1.0e5, core, [[1.0] * 12])
    assert rates["forward"] is True and rates["error"] is None


def test_the_recorded_solver_mask_is_the_only_phase_source(setup, monkeypatch):
    """With conventions.solver_gas_mask the solver's own verdict decides:
    the r147 fixture the reconstruction calls ambiguous (the spawned pool's
    proxy dgt first) is dead with dgt condensed by the solver and live with
    it gas; neither the reconstruction nor any enumeration runs."""
    def no_reconstruction(*_a, **_k):
        raise AssertionError("the reconstruction ran")
    calls = []
    pressures = mp._effective_pressures

    def counted(*a, **k):
        calls.append(1)
        return pressures(*a, **k)
    monkeypatch.setattr(mp, "_h_loss_radical_daughter", no_reconstruction)
    monkeypatch.setattr(mp, "_effective_pressures", counted)
    verdicts = {}
    for gas in (False, True):
        census, refused = _phase_run(
            setup, members=[("dgt", PROPANE, LATER), ("dgt(4)", BUTANE_ADJ, LATER)], aux="N2/0.0/ dgt(4)/0.0/",
            bookkeeping={"dgt": DGT + ["dgt(4)"]}, condensed=DGT,
            solver_mask={proxy_of(): False, **{c: False for c in DGT[:3]},
                         **{mp.moment_column(f"{LABEL}_mu{k}"): False for k in range(3)}, "dgt": gas})
        assert census["phase_source"] == "solver_mask" and census["rate_errors"] == {}
        verdicts[gas] = (refused, census["zero_rate_rows"])
    assert verdicts == {False: (True, ["r0"]), True: (False, [])}
    assert len(calls) == 2                      # one gas set per run: no enumeration


def test_a_recorded_mask_of_another_core_is_refused(setup):
    """The recorded mask must be the mask of the profile's core, in its
    order: it is static over the simulation that wrote the profile."""
    out = setup["tmp"] / "out"
    produce(setup, out)
    constant_run(out, setup["spec"], solver_mask=[("N2", True), ("H(1)", True)])
    with pytest.raises(mp.ProducerError, match="solver_gas_mask"):
        produce(setup, out, mode="convert")


# --------------------------------------------------------------------------
# I-178 (spar r156): the recorded mask bound to its reactor, profile and
# structures; the reconstruction's base labels
# --------------------------------------------------------------------------
def _bound_run(setup, extra_gas=(), **kw):
    """The r147 fixture with a recorded mask (dgt condensed); ``kw`` binds it;
    ``extra_gas``: more gas columns (no amount)."""
    out = setup["tmp"] / f"out{len(list(setup['tmp'].glob('out*')))}"
    proxy = proxy_of()
    row = _row(equation=f"{proxy}(+M)=>{proxy}+X0(+M)", rate=(1.0, 0.0, 0.0),
               aux=["LOW /  1.0 0.0 0.0 /", "N2/0.0/ dgt(4)/0.0/"], sidecar_kinetics=None)
    produce(setup, out)
    write_run(out, times=PROFILE_T, rows=[row], gas=ZERO_GAS + [("dgt", PROPANE, LATER), ("dgt(4)", BUTANE_ADJ, LATER)]
              + [(name, adj, lambda t: 0.0) for name, adj in extra_gas],
              pools=[(LABEL, drift(declared(setup["spec"])), 0.0), ("dgt", lambda t: (0.0, 0.0, 0.0), 0.0)],
              bookkeeping={"dgt": DGT + ["dgt(4)"]}, condensed=DGT, spawned=["dgt"],
              solver_mask={proxy: False, **{c: False for c in DGT[:3]},
                           **{mp.moment_column(f"{LABEL}_mu{k}"): False for k in range(3)}, "dgt": False}, **kw)
    return out


def test_a_bound_mask_of_reactor_1_and_its_profile_is_used(setup):
    out = _bound_run(setup)
    with pytest.raises(mp.ProducerError, match="inactive_chemistry"):
        produce(setup, out, mode="convert")
    assert json.loads((out / "producer.json").read_text())["census"]["phase_source"] == "solver_mask"


def test_the_mask_of_reactor_2_is_refused_for_reactor_1s_profile(setup):
    """Spar r156 item 1: the polymer reactor is the second reaction system;
    the producer reads reactor 1's profile, which no recorded mask names."""
    out = _bound_run(setup, mask_reactor=2, mask_profile="solver/simulation_2_8.csv")
    with pytest.raises(mp.ProducerError, match="reactor 1"):
        produce(setup, out, mode="convert")


def test_a_mask_naming_another_profile_is_refused(setup):
    out = _bound_run(setup, mask_profile="solver/simulation_1_7.csv")
    with pytest.raises(mp.ProducerError, match="simulation_1_8.csv"):
        produce(setup, out, mode="convert")


def test_an_isomer_under_a_recorded_label_is_refused(setup):
    """Spar r156 item 2: the mask was the solver's for another structure
    under the label dgt(4) (propane, not the butane of the Chemkin file)."""
    out = _bound_run(setup, mask_structures={"dgt(4)": _structure(Molecule(smiles="CCC"))})
    with pytest.raises(mp.ProducerError, match="dgt\(4\) has another structure"):
        produce(setup, out, mode="convert")


def test_an_unbound_mask_record_is_refused(setup):
    """A record without reactor, profile and structures cannot be bound."""
    out = _bound_run(setup, mask_record={"source": "solver", "core": [{"species": "N2", "gas": True}]})
    with pytest.raises(mp.ProducerError, match="solver_gas_mask"):
        produce(setup, out, mode="convert")


def _c6h13_run(setup, comments=None):
    """Spar r156 item 3: a spawned pool named C6H13 (members C6H13, the
    first, and the butane C6H13(9), the only collider weight) beside a
    non-pool hexyl radical, an H-loss daughter of the root proxy, whose
    Chemkin identifier is its formula, C6H13(5) (amount 0)."""
    members = [("C6H13", PROPANE, lambda t: 0.0), ("C6H13(5)", HEXYL, lambda t: 0.0),
               ("C6H13(9)", BUTANE_ADJ, LATER)]
    pool = ["C6H13_mu0", "C6H13_mu1", "C6H13_mu2", "C6H13"]
    return _phase_run(setup, members=members, aux="N2/0.0/ C6H13/0.0/", bookkeeping={"C6H13": pool + ["C6H13(9)"]},
                      spawned=("C6H13",), condensed=pool,
                      pools=[(LABEL, drift(declared(setup["spec"])), 0.0), ("C6H13", lambda t: (0.0, 0.0, 0.0), 0.0)],
                      species_comments=comments)


def test_a_formula_identifier_is_no_base_label(setup):
    """Without its own label the hexyl is its own unit, never the pool's:
    C6H13(9) stays gas and the row live."""
    census, refused = _c6h13_run(setup)
    assert census["phase_source"] == "reconstructed"
    assert not refused and census["live_rows"] == ["r0"] and census["rate_errors"] == {}


def test_a_long_own_label_is_the_base_label(setup):
    """RMG's verbose SPECIES line gives the hexyl's own label: its base is
    that label, not the pool's; the row stays live."""
    census, refused = _c6h13_run(setup, comments={"C6H13(5)": "hexyl_primary_radical(5)"})
    assert not refused and census["live_rows"] == ["r0"] and census["rate_errors"] == {}


def test_an_own_label_of_the_pools_base_propagates(setup):
    """Its own label C6H13(5) has the pool's base: MoM's section D condenses
    every species of base C6H13, so C6H13(9)'s phase is the hexyl's unit."""
    census, refused = _c6h13_run(setup, comments={"C6H13(5)": "C6H13(5)"})
    assert refused and "gas mask ambiguous for ['stage1:C6H13(5)']" in census["rate_errors"]["r0"]


# --------------------------------------------------------------------------
# I-178 rework 1 (spar r160): species identity by RMG's isomorphism
# --------------------------------------------------------------------------
def _structure(molecule):
    return {"adjacency_list": molecule.to_adjacency_list(), "multiplicity": molecule.multiplicity}


def _with_multiplicity(smiles, multiplicity):
    molecule = Molecule(smiles=smiles)
    molecule.multiplicity = multiplicity
    return molecule


@pytest.mark.parametrize("recorded, chemkin", [
    (Molecule(smiles="Oc1ccccn1"), Molecule(smiles="O=c1cccc[nH]1")),                 # tautomers
    (_with_multiplicity("[CH2]C[CH]C", 1), _with_multiplicity("[CH2]C[CH]C", 3)),     # singlet, triplet
], ids=["tautomer", "multiplicity"])
def test_a_distinct_species_under_a_recorded_label_is_refused(setup, recorded, chemkin):
    """Spar r160 item 1: one augmented InChI, two species: refused."""
    out = _bound_run(setup, extra_gas=[("T(10)", chemkin.to_adjacency_list())],
                     mask_structures={"T(10)": _structure(recorded)})
    with pytest.raises(mp.ProducerError, match=r"T\(10\) has another structure"):
        produce(setup, out, mode="convert")


def test_another_resonance_drawing_is_the_recorded_species(setup):
    """Spar r160 item 2: crotyl recorded C=C[CH]C, written [CH2]C=CC in the
    dictionary: one species, the mask is used."""
    out = _bound_run(setup, extra_gas=[("T(10)", Molecule(smiles="[CH2]C=CC").to_adjacency_list())],
                     mask_structures={"T(10)": _structure(Molecule(smiles="C=C[CH]C"))})
    with pytest.raises(mp.ProducerError, match="inactive_chemistry"):
        produce(setup, out, mode="convert")
    assert json.loads((out / "producer.json").read_text())["census"]["phase_source"] == "solver_mask"


@pytest.mark.parametrize("structure", [None, {"adjacency_list": "1 Q u9 p0 c0", "multiplicity": 2}],
                         ids=["absent", "unreadable"])
def test_a_recorded_structure_that_cannot_be_read_is_refused_by_name(setup, structure):
    """Spar r160 item 3: a structure that cannot be read equals nothing."""
    out = _bound_run(setup, mask_structures={"dgt(4)": structure})
    with pytest.raises(mp.ProducerError, match=r"dgt\(4\).*cannot be read"):
        produce(setup, out, mode="convert")

