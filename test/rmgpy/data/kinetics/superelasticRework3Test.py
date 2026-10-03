#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""Regression witnesses for declaration preservation before conversion."""
from copy import deepcopy
import importlib.util
from math import exp
from pathlib import Path
from types import SimpleNamespace

import pytest

from rmgpy.exceptions import ExcitedSpeciesThermoError

from rmgpy.rmg.model import ReactionModel
from rmgpy.kinetics import Arrhenius, TwoTemperaturePlasma, PDepArrhenius
from rmgpy.exceptions import (NonEquilibriumReverseRateError, KineticsError,
                              ResolvedStateTrainingError, SpeciesIdentityError, StateProvenanceError)
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.family import KineticsFamily, TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
from rmgpy.chemkin import (write_kinetics_entry, read_kinetics_entry, get_species_identifier,
                           render_chemkin_file, save_species_dictionary, load_chemkin_file,
                           _process_duplicate_reactions)

# Load helpers by path, independent of pytest's import mode and PYTHONPATH.
ROOT = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location(
    'superelastic_rework2_helpers', ROOT / 'test/rmgpy/data/kinetics/superelasticRework2Test.py')
helpers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)


def exact_state_fixture_equilibrium_constant(temperature):
    """Return Kc for N2v1 -> N2 directly from the fixture NASA coefficients."""
    # N2: (a1, a6, a7) = (4.5, -1000, 3); N2v1: (3.5, 2400, 3).
    return temperature * exp(3400 / temperature - 1)


def test_chemkin_duplicate_roundtrip_retains_independent_reverse(monkeypatch, tmp_path):
    helpers._install_exact_state_thermo(monkeypatch)
    forward, reverse = helpers.explicit_pair(cls=LibraryReaction)
    helpers._load_exact_state_thermo(forward)
    second = deepcopy(forward)
    second.reactants, second.products = forward.reactants, forward.products
    for rxn, a in zip((forward, second, reverse), (1e5, 2e5, 7e5)):
        rxn.library = 'test'
        rxn.kinetics = Arrhenius(A=(a, 's^-1'), n=0, Ea=(0, 'J/mol'), T0=(1, 'K'))
    species = forward.reactants + forward.products
    path, dictionary = tmp_path / 'chem.inp', tmp_path / 'dictionary.txt'
    path.write_text(render_chemkin_file(species, [forward, second, reverse]))
    save_species_dictionary(str(dictionary), species)
    _, loaded = load_chemkin_file(str(path), str(dictionary))
    assert len(loaded) == 2
    by_direction = {}
    for reaction in loaded:
        excited_reactant = reaction.reactants[0].molecule[0].has_resolved_state()
        excited_product = reaction.products[0].molecule[0].has_resolved_state()
        direction = excited_reactant, excited_product
        assert excited_reactant != excited_product
        assert direction not in by_direction
        assert reaction.reversible is False
        by_direction[direction] = reaction.kinetics.get_rate_coefficient(1000)
    assert set(by_direction) == {(True, False), (False, True)}
    assert by_direction[True, False] == pytest.approx(3e5)
    assert by_direction[False, True] == pytest.approx(7e5)


def test_chemkin_duplicate_roundtrip_refuses_unattested_thermo(tmp_path):
    forward, reverse = helpers.explicit_pair(cls=LibraryReaction)
    second = deepcopy(forward)
    second.reactants, second.products = forward.reactants, forward.products
    for rxn, a in zip((forward, second, reverse), (1e5, 2e5, 7e5)):
        rxn.library = 'test'
        rxn.kinetics = Arrhenius(A=(a, 's^-1'), n=0, Ea=(0, 'J/mol'), T0=(1, 'K'))
    species = forward.reactants + forward.products
    path, dictionary = tmp_path / 'chem.inp', tmp_path / 'dictionary.txt'
    with pytest.raises(ExcitedSpeciesThermoError, match='N2v1'):
        render_chemkin_file(species, [forward, second, reverse])
    assert not path.exists() and not dictionary.exists()


def test_unresolved_irreversible_duplicates_remain_irreversible(tmp_path):
    """Correct the base bug that reloaded this unresolved pair as reversible."""
    first, _ = helpers.explicit_pair(resolved=False, cls=LibraryReaction)
    second = deepcopy(first)
    second.reactants, second.products = first.reactants, first.products
    for rxn, rate in zip((first, second), (1e5, 2e5)):
        rxn.library = 'test'
        rxn.kinetics = Arrhenius(A=(rate, 's^-1'), n=0, Ea=(0, 'J/mol'), T0=(1, 'K'))
    species = first.reactants + first.products
    path, dictionary = tmp_path / 'chem.inp', tmp_path / 'dictionary.txt'
    path.write_text(render_chemkin_file(species, [first, second]))
    save_species_dictionary(str(dictionary), species)
    _, loaded = load_chemkin_file(str(path), str(dictionary))
    assert len(loaded) == 1
    assert not loaded[0].reversible
    assert loaded[0].kinetics.get_rate_coefficient(1000) == pytest.approx(3e5)


@pytest.mark.parametrize('reversible', [False, True])
def test_chemkin_plog_duplicate_policy(reversible):
    rxn = helpers.make_reaction(cls=LibraryReaction, reversible=reversible)
    rxn.duplicate = True
    rxn.library = 'test'
    rxn.kinetics = PDepArrhenius(pressures=([1, 100], 'bar'), arrhenius=[rxn.kinetics] * 2)
    second = helpers.make_reaction(cls=LibraryReaction, reversible=reversible)
    second.kinetics = rxn.kinetics
    second.duplicate, second.library = True, rxn.library
    second.reactants, second.products = rxn.reactants, rxn.products
    reactions = [rxn, second]
    if reversible:
        with pytest.raises(NonEquilibriumReverseRateError, match='N2v1.*vibrationallevel 1'):
            _process_duplicate_reactions(reactions)
    else:
        _process_duplicate_reactions(reactions)
        assert len(reactions) == 1
        assert not reactions[0].reversible
        assert reactions[0].kinetics.get_rate_coefficient(1000, 1e5) == pytest.approx(2e8)


def test_chemkin_reader_rejects_original_resolved_tdep():
    rxn = helpers.make_reaction(reversible=False)
    electron = helpers.Species(label='e').from_adjacency_list('1 e u1 p0 c-1')
    species = rxn.reactants + rxn.products + [electron]
    entry = write_kinetics_entry(rxn, species, verbose=False).replace('=>', '<=>')
    entry = '\n'.join(line.split('!')[0] for line in entry.splitlines())
    units = ['', 's^-1', 'cm^3/(mol*s)', 'cm^6/(mol^2*s)']
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1.*vibrationallevel 1'):
        read_kinetics_entry(entry, {get_species_identifier(
            s): s for s in species}, units, units, 'cal/mol')


@pytest.mark.parametrize('exporter', ['chemkin', 'cantera-direct', 'cantera1', 'cantera2', 'rms'])
@pytest.mark.parametrize('declared_subclass', [False, True])
def test_exporter_checks_original_declaration(exporter, declared_subclass):
    rxn = helpers.make_reaction(kinetics=helpers.ElectronArrhenius(
        A=(1e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K')) if declared_subclass else None)
    electron = helpers.Species(label='e').from_adjacency_list('1 e u1 p0 c-1')
    species = rxn.reactants + rxn.products + [electron]
    if exporter == 'chemkin':
        def call(): return write_kinetics_entry(rxn, species, verbose=False)
    elif exporter == 'cantera-direct':
        def call(): return rxn.to_cantera()
    elif exporter == 'cantera1':
        from rmgpy.yaml_cantera1 import reaction_to_dicts
        def call(): return reaction_to_dicts(rxn, species)
    elif exporter == 'cantera2':
        from rmgpy.yaml_cantera2 import reaction_to_dict_list
        def call(): return reaction_to_dict_list(rxn, species)
    else:
        from rmgpy.yaml_rms import obj_to_dict
        def call(): return obj_to_dict(rxn, species, names=[s.label for s in species])
    original = rxn.kinetics
    with pytest.raises(NonEquilibriumReverseRateError, match='N2v1.*vibrationallevel 1'):
        call()
    assert rxn.kinetics is original


def test_rms_preserves_thermal_irreversible_declaration():
    from rmgpy.yaml_rms import obj_to_dict
    rxn = helpers.make_reaction(reversible=False, kinetics=Arrhenius(
        A=(1e5, 's^-1'), n=0, Ea=(0, 'J/mol'), T0=(1, 'K')))
    species = rxn.reactants + rxn.products
    data = obj_to_dict(rxn, species, names=[s.label for s in species])
    assert data['reversible'] is False


@pytest.mark.parametrize('reversible', [False, True])
def test_arkane_refuses_before_generation(reversible):
    from arkane.kinetics import KineticsJob
    from rmgpy.species import TransitionState
    from rmgpy.statmech import Conformer, IdealGasTranslation, LinearRotor, HarmonicOscillator
    rxn = helpers.make_reaction(reversible=reversible)

    def conformer(e0):
        return Conformer(E0=e0, modes=[IdealGasTranslation(mass=(28, 'amu')),
                         LinearRotor(inertia=(10, 'amu*angstrom^2'), symmetry=2),
                         HarmonicOscillator(frequencies=([1000], 'cm^-1'))])
    for species in rxn.reactants + rxn.products:
        species.conformer = conformer(species.thermo.E0)
    rxn.transition_state = TransitionState(conformer=conformer((20000, 'J/mol')))
    original = rxn.kinetics
    job = KineticsJob(rxn)
    with pytest.raises(ExcitedSpeciesThermoError, match='N2v1'):
        job.generate_kinetics()
    assert rxn.kinetics is original
    assert not job.usedTST


@pytest.mark.parametrize('retry', [False, True])
def test_family_multiple_candidates_reach_direction_guard(monkeypatch, retry):
    import rmgpy.data.kinetics.family as module
    forward, reverse = helpers.explicit_pair(cls=TemplateReaction)
    forward.template = reverse.template = ['test']
    # The structural product shortcut aligns candidates on their family-forward
    # basis. Opposite storage orientations therefore need opposite is_forward.
    forward.is_forward, reverse.is_forward = True, False
    family = KineticsFamily(label='test')
    family.own_reverse = True
    family.forbidden = object()
    candidates = iter([[], [reverse, forward]] if retry else [[reverse, forward]])
    monkeypatch.setattr(family, '_generate_reactions', lambda *args, **kwargs: next(candidates))
    monkeypatch.setattr(module, 'find_degenerate_reactions',
                        lambda reactions, *args, **kwargs: reactions)
    with pytest.raises(KineticsError, match='Did not find reverse|Found multiple reverse'):
        family.add_reverse_attribute(forward)
    assert forward.reverse is None


def test_model_comparison_keeps_opposite_irreversible_channels(monkeypatch):
    from rmgpy.tools import diffmodels
    forward, reverse = helpers.explicit_pair()
    monkeypatch.setattr(diffmodels.plt, 'show', lambda: None)
    model1, model2 = ReactionModel(reactions=[forward]), ReactionModel(reactions=[reverse])
    diffmodels.compare_model_kinetics(model1, model2)
    assert model1.reactions == [forward]
    assert model2.reactions == [reverse]
    diffmodels.plt.close('all')


def test_library_marked_duplicates_keep_opposite_channels():
    from rmgpy.data.kinetics.library import KineticsLibrary
    forward, reverse = helpers.explicit_pair()
    forward.duplicate = reverse.duplicate = True
    for reaction in (forward, reverse):
        reaction.kinetics = PDepArrhenius(pressures=([1, 100], 'bar'),
                                          arrhenius=[reaction.kinetics] * 2)
    library = KineticsLibrary(label='marked-pair')
    library.entries = {i: Entry(index=i, item=rxn, data=rxn.kinetics)
                       for i, rxn in enumerate((forward, reverse), 1)}
    library.convert_duplicates_to_multi()
    assert len(library.entries) == 2
    assert {id(entry.item) for entry in library.entries.values()} == {id(forward), id(reverse)}


@pytest.mark.parametrize('own_reverse', [False, True])
def test_family_training_extraction_reaches_reverse_helper(monkeypatch, own_reverse):
    from rmgpy.molecule import Group
    forward, reverse = helpers.explicit_pair()
    forward.kinetics = helpers.ElectronArrhenius(A=(1e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K'))
    for s in forward.reactants + forward.products:
        s.molecule[0].atoms[0].label = '*1'
    # State-aware molecular templates deliberately refuse ordinary ground-state
    # groups. Stub only that structural matching decision so this fixture enters
    # the actual training extraction/reverse branches without adding group support.

    class MatchingMolecule:
        def __init__(self, molecule, matches):
            self.molecule, self.matches = molecule, matches

        def __getattr__(self, name):
            return getattr(self.molecule, name)

        def __deepcopy__(self, memo):
            return type(self)(deepcopy(self.molecule, memo), self.matches)

        def copy(self, deep=False):
            return type(self)(self.molecule.copy(deep=deep), self.matches)

        def is_subgraph_isomorphic(self, *args, **kwargs):
            return self.matches

    for species, matches in ((forward.reactants[0], own_reverse), (forward.products[0], True)):
        species.molecule = [MatchingMolecule(species.molecule[0], matches)]
    root = Group().from_adjacency_list('1 *1 N u0')
    family = KineticsFamily(label='test')
    family.own_reverse = own_reverse
    family.reverse_map = {'*1': '*1'}
    monkeypatch.setattr(family, 'get_root_template', lambda: [Entry(item=root)])
    monkeypatch.setattr(family, 'get_training_depository', lambda: SimpleNamespace(
        entries={1: Entry(item=forward, data=forward.kinetics)}))
    # Group training is intentionally unavailable for resolved species on the
    # rebased plasma, so its more specific refusal wins before reversal logic.
    with pytest.raises(ResolvedStateTrainingError, match='test.*resolved species'):
        family.get_training_set(estimate_thermo=False, remove_degeneracy=False,
                                fix_labels=True, get_reverse=own_reverse)


@pytest.mark.parametrize('exporter', ['library-entry', 'legacy-library', 'arkane-chemkin'])
def test_other_serializers_check_original_before_erasure(exporter, tmp_path, monkeypatch):
    import io
    from rmgpy.data.kinetics.common import save_entry
    from rmgpy.data.kinetics.library import KineticsLibrary
    from arkane.kinetics import KineticsJob

    class FormattableReaction(helpers.Reaction):
        def __format__(self, spec):
            return format(str(self), spec)
    rxn = helpers.make_reaction(cls=FormattableReaction, kinetics=helpers.ElectronArrhenius(
        A=(1e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K')))
    entry = Entry(index=1, item=rxn, data=rxn.kinetics)
    if exporter == 'library-entry':
        def call(): return save_entry(io.StringIO(), entry)
    elif exporter == 'legacy-library':
        library = KineticsLibrary(label='declared')
        library.entries = {1: entry}
        monkeypatch.setattr(library, 'get_species',
                            lambda: {sp.label: sp for sp in rxn.reactants + rxn.products})

        def call(): return library.save_old(str(tmp_path))
    else:
        prepare_transition_state(rxn)
        def call(): return KineticsJob(rxn).write_chemkin(str(tmp_path))
    error = SpeciesIdentityError if exporter == 'legacy-library' else NonEquilibriumReverseRateError
    match = 'Old-style kinetics libraries' if exporter == 'legacy-library' else 'N2v1.*vibrationallevel 1'
    with pytest.raises(error, match=match):
        call()


def test_arkane_chemkin_checks_and_writes_declared_irreversible_direction(tmp_path, monkeypatch):
    from arkane.kinetics import KineticsJob

    class FormattableReaction(helpers.Reaction):
        def __format__(self, spec):
            return format(str(self), spec)

    rxn = helpers.make_reaction(
        reversible=False,
        cls=FormattableReaction,
        kinetics=helpers.ElectronArrhenius(
            A=(1e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K')),
    )
    prepare_transition_state(rxn)
    # I-308 independently refuses resolved structures in ordinary Chemkin.
    # Isolate this writer's direction contract so the I-309 chokepoint and
    # emitted arrow are exercised together.
    monkeypatch.setattr('rmgpy.export.refuse_resolved_species', lambda *args: None)

    KineticsJob(rxn).write_chemkin(str(tmp_path))

    assert ' => ' in (tmp_path / 'chem.inp').read_text()


def test_policy_thermodynamic_property_remains_available(monkeypatch):
    # Thermodynamic properties require an exact-state library source; the
    # independent Te reversal policy still refuses the original declaration.
    helpers._install_exact_state_thermo(monkeypatch)
    rxn = helpers.make_reaction()
    helpers._load_exact_state_thermo(rxn)
    scalar = rxn.get_equilibrium_constant(300)
    temperatures = helpers.np.array([300., 1000.])
    vector = rxn.get_equilibrium_constants(temperatures)
    expected = [exact_state_fixture_equilibrium_constant(temperature)
                for temperature in temperatures]
    assert scalar == pytest.approx(expected[0])
    assert vector == pytest.approx(expected)
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        rxn.generate_reverse_rate_coefficient()


def test_policy_thermodynamic_property_refuses_unattested_thermo():
    rxn = helpers.make_reaction()
    with pytest.raises(ExcitedSpeciesThermoError, match='N2v1'):
        rxn.get_equilibrium_constant(300)
    with pytest.raises(ExcitedSpeciesThermoError, match='N2v1'):
        rxn.get_equilibrium_constants(helpers.np.array([300., 1000.]))
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        rxn.generate_reverse_rate_coefficient()


def test_serialized_irreversible_declaration_survives_admission():
    import pickle
    forward, reverse = helpers.explicit_pair()
    from rmgpy.rmg.model import CoreEdgeReactionModel
    for original in (forward, reverse):
        restored = pickle.loads(pickle.dumps(original))
        assert restored.reversible is False
        assert restored.has_resolved_species()
        restored.check_resolved_species_reversibility()
        model = CoreEdgeReactionModel()
        model.add_reaction_to_core(restored)
        assert model.core.reactions == [restored]


def test_family_rules_from_training_reaches_reverse_helper(monkeypatch):
    from rmgpy.exceptions import UndeterminableKineticsError
    import rmgpy.rmg.main as main
    rxn = helpers.make_reaction(reversible=False, kinetics=helpers.ElectronArrhenius(
        A=(1e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K')))
    family = KineticsFamily(label='test')
    family.auto_generated = False
    family.rules = SimpleNamespace(get_entries=lambda: [], entries={})
    monkeypatch.setattr(main, 'determine_procnum_from_ram', lambda: 1)
    monkeypatch.setattr(family, 'get_training_depository', lambda: SimpleNamespace(
        entries={1: Entry(index=1, item=rxn, data=rxn.kinetics)}))

    from rmgpy.molecule import Group
    calls = []

    def reverse_orientation(reaction):
        calls.append(reaction)
        if len(calls) == 1:
            raise UndeterminableKineticsError(reaction)
        return [Entry(label='test', item=Group().from_adjacency_list('1 N u0'))]
    monkeypatch.setattr(family, 'get_reaction_template', reverse_orientation)
    monkeypatch.setattr(family, 'calculate_degeneracy', lambda reaction: 1)
    thermo = SimpleNamespace(get_thermo_data=lambda *args, **kwargs: rxn.reactants[0].thermo)
    # Upstream refuses resolved training before any automatic orientation.
    with pytest.raises(ResolvedStateTrainingError, match='test.*resolved species'):
        family.add_rules_from_training(thermo_database=thermo)


def test_reconstruct_source_reaches_reverse_helper():
    from rmgpy.data.kinetics.database import KineticsDatabase
    rxn = helpers.make_reaction(reversible=False)
    training = Entry(index=1, item=rxn, data=rxn.kinetics)
    # Source reconstruction has no state-matched provenance, and therefore
    # refuses before it could attempt an equilibrium reverse.
    with pytest.raises(ExcitedSpeciesThermoError, match='N2v1'):
        KineticsDatabase().reconstruct_kinetics_from_source(
            rxn, {'Training': ['test', training, True]})


def test_family_source_extraction_keeps_direction(monkeypatch):
    forward, reverse = helpers.explicit_pair()
    forward.kinetics.comment = reverse.kinetics.comment = 'Matched reaction 1 training-pair in training'
    training = Entry(index=1, label='training-pair', item=forward, data=forward.kinetics)
    family = KineticsFamily(label='test')
    monkeypatch.setattr(family, 'get_training_depository',
                        lambda: SimpleNamespace(entries={1: training}))
    assert family.extract_source_from_comments(forward)[1][2] is False
    assert family.extract_source_from_comments(reverse)[1][2] is True


@pytest.mark.parametrize('stage', ['initial', 'perturbed'])
def test_arkane_sensitivity_reaches_reverse_helper(tmp_path, monkeypatch, stage):
    from arkane.sensitivity import KineticsSensitivity
    from rmgpy.quantity import Temperature
    rxn = helpers.make_reaction(kinetics=helpers.ElectronArrhenius(
        A=(1e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K')))
    job = SimpleNamespace(reaction=rxn, sensitivity_conditions=[
                          Temperature(300, 'K')], execute=lambda **kwargs: None)
    if stage == 'initial':
        with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
            KineticsSensitivity(job, str(tmp_path))
    else:
        sensitivity = KineticsSensitivity.__new__(KineticsSensitivity)
        sensitivity.job = job
        sensitivity.conditions = job.sensitivity_conditions
        sensitivity.f_sa_rates = {}
        sensitivity.r_sa_rates = {}
        monkeypatch.setattr(sensitivity, 'perturb', lambda species: None)
        with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
            sensitivity.execute()


def test_collision_limit_reaches_reverse_helper():
    from rmgpy.transport import TransportData
    rxn = helpers.make_reaction()
    rxn.products.append(helpers.Species(label='Ar').from_smiles('[Ar]'))
    for species in rxn.products:
        species.transport_data = TransportData(sigma=(3.5, 'angstrom'), epsilon=(100, 'K'))
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        rxn.check_collision_limit_violation(300, 1000, 1e5, 1e5)


def test_html_export_reaches_reverse_helper(tmp_path):
    from rmgpy.rmg.output import save_diff_html
    forward, reverse = helpers.explicit_pair(cls=LibraryReaction)
    forward.library = reverse.library = 'test'
    with pytest.raises(SpeciesIdentityError, match='N2v1'):
        save_diff_html(str(tmp_path / 'diff.html'), [], [], [], [[forward, reverse]], [], [])


def test_direct_rms_export_checks_original_before_julia(monkeypatch):
    from rmgpy.rmg import reactionmechanismsimulator_reactors as module
    rxn = helpers.make_reaction(kinetics=helpers.ElectronArrhenius(
        A=(1e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K')))
    touched = []

    class JuliaBarrier:
        def __getattr__(self, name):
            touched.append(name)
            raise AssertionError('Julia reached before original declaration check')
    monkeypatch.setattr(module, 'Main', JuliaBarrier())
    with pytest.raises(SpeciesIdentityError, match='N2v1'):
        module.to_rms(rxn, species_names=['N2v1', 'N2'], rms_species_list=[None, None])
    assert touched == []


def test_explicit_cantera_pair_reload_preserves_directions(monkeypatch):
    import cantera as ct
    import yaml
    import rmgpy.data.rmg as data_rmg
    from rmgpy.yaml_cantera2 import generate_cantera_data
    thermo = ThermoDatabase()
    library = ThermoLibrary(label='ExactStateFixture')
    library.load(str(ROOT / 'test/rmgpy/test_data/excited_states/thermo.py'), thermo.local_context, {})
    thermo.libraries = {library.label: library}
    thermo.library_order = [library.label]
    monkeypatch.setattr(data_rmg, 'database', SimpleNamespace(thermo=thermo, solvation=None))
    forward, reverse = helpers.explicit_pair()
    for species in forward.reactants + forward.products:
        species.thermo = None
        species.get_thermo_data()
    electron = helpers.Species(label='e').from_adjacency_list('1 e u1 p0 c-1')
    electron.thermo = forward.products[0].thermo
    for reaction, a in ((forward, 1e5), (reverse, 2e5)):
        reaction.reactants.insert(0, electron)
        reaction.products.insert(0, electron)
        reaction.kinetics = TwoTemperaturePlasma(A=(a, 'm^3/(mol*s)'), n=1)
    species = [electron, forward.reactants[1], forward.products[1]]
    gas = ct.Solution(yaml=yaml.safe_dump(generate_cantera_data(species, [forward, reverse], is_plasma=True)),
                      name='gas', transport_model=None)
    gas.TP = 300, 1e5
    gas.Te = 12000
    assert len(gas.reactions()) == 2
    assert all(not reaction.reversible for reaction in gas.reactions())
    assert list(gas.reverse_rate_constants) == [0., 0.]
    for reaction, rate in zip(gas.reactions(), gas.forward_rate_constants):
        excited_reactant = any(label.startswith('N2v1') for label in reaction.reactants)
        excited_product = any(label.startswith('N2v1') for label in reaction.products)
        assert excited_reactant != excited_product
        assert rate == pytest.approx(1.2e12 if excited_reactant else 2.4e12)


def test_explicit_cantera_pair_reload_refuses_unattested_thermo():
    from rmgpy.yaml_cantera2 import generate_cantera_data
    forward, reverse = helpers.explicit_pair()
    electron = helpers.Species(label='e').from_adjacency_list('1 e u1 p0 c-1')
    electron.thermo = forward.products[0].thermo
    for reaction, a in ((forward, 1e5), (reverse, 2e5)):
        reaction.reactants.insert(0, electron)
        reaction.products.insert(0, electron)
        reaction.kinetics = TwoTemperaturePlasma(A=(a, 'm^3/(mol*s)'), n=1)
    species = [electron, forward.reactants[1], forward.products[1]]
    with pytest.raises(ExcitedSpeciesThermoError, match='N2v1'):
        generate_cantera_data(species, [forward, reverse], is_plasma=True)


@pytest.mark.parametrize('read_comments', [False, True])
def test_chemkin_comment_reader_checks_returned_declaration(read_comments):
    from rmgpy.chemkin import read_reaction_comments
    rxn = helpers.make_reaction()
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        read_reaction_comments(rxn, '', read=read_comments)


def test_arkane_direct_output_checks_original_before_writing(tmp_path):
    from arkane.kinetics import KineticsJob
    rxn = helpers.make_reaction(kinetics=helpers.ElectronArrhenius(
        A=(1e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K')))
    prepare_transition_state(rxn)
    job = KineticsJob(rxn)
    assert not job.usedTST
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        job.write_output(str(tmp_path))
    assert not (tmp_path / 'output.py').exists()


def test_arkane_kinetics_library_checks_exact_rate(tmp_path):
    from arkane.output import save_kinetics_lib

    reaction = helpers.make_reaction(kinetics=helpers.ElectronArrhenius(
        A=(1e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K')))
    output = tmp_path / 'kinetics-library'
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        save_kinetics_lib([reaction], str(output), 'test', '')
    assert not output.exists()


def test_arkane_get_libraries_checks_exact_rate_before_entry(monkeypatch):
    import arkane.main

    helpers._install_exact_state_thermo(monkeypatch)
    irreversible, _ = helpers.explicit_pair()
    helpers._load_exact_state_thermo(irreversible)
    expected_rate = irreversible.kinetics.get_rate_coefficient(300, 12000)
    assert irreversible.reversible is False
    app = arkane.main.Arkane()
    app.species_dict = {
        species.label: species
        for species in irreversible.reactants + irreversible.products
    }
    app.reaction_dict = {'channel': irreversible}

    _, kinetics_library, _ = app.get_libraries()
    entry = kinetics_library.entries[1]
    assert entry.item is irreversible
    assert entry.data is irreversible.kinetics
    assert entry.item.reactants[0].molecule[0].has_resolved_state()
    assert not entry.item.products[0].molecule[0].has_resolved_state()
    assert entry.item.reversible is False
    assert entry.data.get_rate_coefficient(300, 12000) == pytest.approx(expected_rate)


def test_arkane_get_libraries_refuses_unattested_thermo(monkeypatch):
    import arkane.main

    irreversible = helpers.make_reaction(reversible=False)
    app = arkane.main.Arkane()
    app.species_dict = {
        species.label: species
        for species in irreversible.reactants + irreversible.products
    }
    app.reaction_dict = {'channel': irreversible}

    with pytest.raises(ExcitedSpeciesThermoError, match='N2v1'):
        app.get_libraries()
    assert app.reaction_dict['channel'] is irreversible
    assert irreversible.reversible is False

    reversible = helpers.make_reaction(reversible=True)
    app.species_dict = {
        species.label: species
        for species in reversible.reactants + reversible.products
    }
    app.reaction_dict = {'channel': reversible}

    def fail_entry_construction(*args, **kwargs):
        raise AssertionError('kinetics Entry built before serialization check')

    monkeypatch.setattr(arkane.main, 'Entry', fail_entry_construction)
    with pytest.raises(
            NonEquilibriumReverseRateError,
            match='N2v1.*vibrationallevel 1'):
        app.get_libraries()


@pytest.mark.parametrize('transport', ['repr', 'pickle', 'copy'])
@pytest.mark.parametrize(
    'cls', [helpers.Reaction, LibraryReaction, TemplateReaction])
def test_object_serializers_check_original_declaration(transport, cls):
    import pickle
    rxn = helpers.make_reaction(cls=cls, kinetics=helpers.ElectronArrhenius(
        A=(1e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K')))
    call = {'repr': lambda: repr(rxn), 'pickle': lambda: pickle.dumps(rxn),
            'copy': lambda: rxn.copy()}[transport]
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        call()


def test_generic_dictionary_reader_checks_original_declaration():
    from rmgpy.rmgobject import recursive_make_object
    rxn = helpers.make_reaction()
    data = {'class': 'Reaction', 'reactants': rxn.reactants, 'products': rxn.products,
            'kinetics': rxn.kinetics, 'reversible': True}
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        recursive_make_object(data, {'Reaction': helpers.Reaction})
    data['reversible'] = False
    restored = recursive_make_object(data, {'Reaction': helpers.Reaction})
    assert restored.reversible is False
    assert isinstance(restored.kinetics, TwoTemperaturePlasma)


def test_generic_dictionary_reader_checks_evaluated_declaration():
    from rmgpy.rmgobject import recursive_make_object
    rxn = helpers.make_reaction()
    # The generic reader also supports repr-style dictionary keys. Supply an
    # explicit constructor scope so this branch receives the original rate.
    text = "Reaction(reactants=reactants, products=products, kinetics=kinetics)"
    scope = {'Reaction': helpers.Reaction, 'reactants': rxn.reactants,
             'products': rxn.products, 'kinetics': rxn.kinetics}
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        recursive_make_object(text, scope)


def prepare_transition_state(reaction):
    from rmgpy.species import TransitionState
    from rmgpy.statmech import Conformer
    reaction.transition_state = TransitionState(conformer=Conformer(E0=(20000, 'J/mol')))


@pytest.mark.parametrize('declared_subclass', [False, True])
def test_arkane_saved_pdep_input_refuses_resolved_electron_rates(tmp_path, declared_subclass):
    from arkane.pdep import PressureDependenceJob
    from rmgpy.pdep.network import Network
    reaction = helpers.make_reaction(reversible=False)
    if declared_subclass:
        reaction.kinetics = helpers.ElectronArrhenius(
            A=(1e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K'))
        reaction.network_kinetics = Arrhenius(
            A=(1e5, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K'))
    prepare_transition_state(reaction)
    reaction.label = 'explicit-channel'
    reaction.transition_state.label = 'TS'
    network = Network(label='resolved-electron', path_reactions=[reaction])
    network.energy_correction = 0.0
    job = PressureDependenceJob(
        network, Tlist=([300, 1000], 'K'), Plist=([1e4, 1e5], 'Pa'))
    with pytest.raises(
            NonEquilibriumReverseRateError,
            match='reaction N2v1.*N2.*Resolved species: N2v1.*vibrationallevel 1'):
        job.save_input_file(str(tmp_path / 'input.py'))


@pytest.mark.parametrize('has_tst_modes', [False, True], ids=['ilt', 'rrkm'])
def test_arkane_saved_pdep_input_checks_serialized_surrogate(tmp_path, has_tst_modes):
    from arkane.pdep import PressureDependenceJob
    from rmgpy.pdep.network import Network
    from rmgpy.statmech import Conformer, HarmonicOscillator

    reactions = list(helpers.explicit_pair())
    for reaction, label, rate in zip(reactions, ('forward', 'reverse'), (1e5, 2e5)):
        reaction.label = label
        reaction.kinetics = Arrhenius(A=(rate, 's^-1'), n=0, Ea=(0, 'J/mol'), T0=(1, 'K'))
        reaction.network_kinetics = helpers.ElectronArrhenius(
            A=(rate, 's^-1'), n=1, Ea=(0, 'J/mol'), T0=(1, 'K'))
        modes = [HarmonicOscillator(frequencies=([1000], 'cm^-1'))] if has_tst_modes else []
        prepare_transition_state(reaction)
        reaction.transition_state.label = 'TS_' + label
        reaction.transition_state.conformer = Conformer(E0=(20000, 'J/mol'), modes=modes)

    network = Network(label='resolved-surrogate', path_reactions=reactions)
    network.energy_correction = 0.0
    job = PressureDependenceJob(
        network, Tlist=([300, 1000], 'K'), Plist=([1e4, 1e5], 'Pa'))
    output = tmp_path / 'input.py'

    with pytest.raises(
            NonEquilibriumReverseRateError,
            match='reaction N2v1.*N2.*Resolved species: N2v1.*vibrationallevel 1'):
        job.save_input_file(str(output))

    assert not output.exists()


@pytest.mark.parametrize('conversion', [False, True])
def test_barrier_height_checks_original_before_conversion(conversion):
    from rmgpy.kinetics import ArrheniusEP

    class ElectronEP(ArrheniusEP):
        uses_electron_temperature = True
    rate = (ElectronEP(A=(1e5, 's^-1'), n=0, alpha=0, E0=(0, 'J/mol'))
            if conversion else helpers.ElectronArrhenius(
                A=(1e5, 's^-1'), n=0, Ea=(0, 'J/mol'), T0=(1, 'K')))
    reaction = helpers.make_reaction(reversible=not conversion, kinetics=rate)
    with pytest.raises(NonEquilibriumReverseRateError, match='vibrationallevel 1'):
        reaction.fix_barrier_height()
    assert reaction.kinetics is rate


def test_direct_cantera_resolved_thermal_irreversible_refuses_state_blind_route():
    reaction = helpers.make_reaction(reversible=False, kinetics=Arrhenius(
        A=(1e5, 's^-1'), n=0, Ea=(0, 'J/mol'), T0=(1, 'K')))
    with pytest.raises(SpeciesIdentityError, match='N2v1'):
        reaction.to_cantera()
