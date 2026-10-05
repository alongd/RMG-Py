#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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

"""Closed kinetics schemas and current-state pdep computation boundaries."""

import importlib
import importlib.util
import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest

import rmgpy.kinetics as kinetics
import rmgpy.rmg.pdep as pdep
from rmgpy.exceptions import NetworkError
from rmgpy.kinetics.uncertainties import RateUncertainty
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.species import Species

_spec = importlib.util.spec_from_file_location('electron_reference_fixtures',
    Path(pdep.__file__).parents[2] / 'test/rmgpy/rmg/electronReferenceTest.py')
helpers = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(helpers)

FORMS = ('method_default', 'method_kwdefault', 'static_default', 'class_default',
         'dunder_class_data')


def carrier_rate(form, reference):
    if form == 'method_default':
        class Rate(kinetics.Arrhenius):
            def material(self, value=reference):
                return value
    elif form == 'method_kwdefault':
        class Rate(kinetics.Arrhenius):
            def material(self, *, value=reference):
                return value
    elif form == 'static_default':
        class Rate(kinetics.Arrhenius):
            @staticmethod
            def material(value=reference):
                return value
    elif form == 'class_default':
        class Rate(kinetics.Arrhenius):
            @classmethod
            def material(cls, value=reference):
                return value
    else:
        class Rate(kinetics.Arrhenius):
            __material__ = reference
    return Rate(A=(1, 's^-1'), n=0, Ea=(0, 'J/mol'))


@pytest.mark.parametrize('form', FORMS)
@pytest.mark.parametrize('electron', (True, False))
def test_review_defaults_and_class_data_refuse(form, electron):
    reaction = helpers.cell_reaction(helpers.CENSUS_CELLS[0], False)
    reference = (Species().from_adjacency_list('1 e u0 p0 c-1') if electron
                 else Species().from_smiles('N#N'))
    reaction.kinetics = carrier_rate(form, reference)
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    before = helpers.admission_state(model)
    with pytest.raises(NetworkError, match='exact type') as error:
        model.register_reaction(reaction)
    assert reaction.library in str(error.value)
    assert helpers.admission_state(model) == before


def shipped_kinetics_classes():
    classes = set()
    for source in Path(kinetics.__file__).parent.glob('*.pyx'):
        module = importlib.import_module('rmgpy.kinetics.' + source.stem)
        classes.update(cls for cls in vars(module).values()
                       if inspect.isclass(cls) and cls.__module__ == module.__name__)
    return classes


def assert_complete_allowlist():
    missing = shipped_kinetics_classes() - set(getattr(pdep, '_KINETICS_REFERENCE_FIELDS', {}))
    assert not missing, 'Missing shipped kinetics classes: ' + ', '.join(
        sorted(cls.__name__ for cls in missing))


def test_shipped_kinetics_allowlist_is_complete():
    assert_complete_allowlist()


def test_static_check_detects_removed_class(monkeypatch):
    monkeypatch.delitem(pdep._KINETICS_REFERENCE_FIELDS, kinetics.Arrhenius)
    with pytest.raises(AssertionError, match='Arrhenius'):
        assert_complete_allowlist()


@pytest.mark.parametrize('cls', helpers.KINETICS_CLASSES, ids=lambda cls: cls.__name__)
def test_shipped_rate_ordinary_network_control(cls):
    reaction = helpers.cell_reaction(helpers.CENSUS_CELLS[0], False)
    reaction.kinetics = helpers.make_rate(cls)
    network = pdep.PDepNetwork(source=reaction.reactants[:])
    if cls is kinetics.EEDFChannel:
        with pytest.raises(NetworkError, match='Electron reactions cannot enter pressure-dependent networks'):
            network.add_path_reaction(reaction)
        assert network.path_reactions == []
        assert pdep._has_electron_participant(reaction)
        return
    network.add_path_reaction(reaction)
    assert network.path_reactions == [reaction]
    assert not pdep._has_electron_participant(reaction)


@pytest.mark.parametrize('cls', (kinetics.TunnelingModel, kinetics.Wigner,
                                kinetics.Eckart, RateUncertainty),
                         ids=lambda cls: cls.__name__)
def test_shipped_auxiliary_kinetics_record_control(cls):
    reaction = helpers.cell_reaction(helpers.CENSUS_CELLS[0], False)
    if cls is RateUncertainty:
        record = cls(mu=0, var=1, Tref=1000, N=1)
        reaction.kinetics.uncertainty = record
    else:
        from rmgpy.species import TransitionState
        if cls is kinetics.Wigner:
            record = cls(frequency=(-1000, 'cm^-1'))
        elif cls is kinetics.Eckart:
            record = cls(frequency=(-1000, 'cm^-1'), E0_reac=(0, 'kJ/mol'),
                         E0_TS=(50, 'kJ/mol'), E0_prod=(0, 'kJ/mol'))
        else:
            record = cls()
        reaction.transition_state = TransitionState(tunneling=record)
    assert list(pdep._reaction_species_references(record)) == []
    network = pdep.PDepNetwork(source=reaction.reactants[:])
    network.add_path_reaction(reaction)
    assert network.path_reactions == [reaction]
    assert not pdep._has_electron_participant(reaction)


ENTRY_POINTS = ('initialize', 'calculate_rate_coefficients', 'set_conditions',
                'calculate_microcanonical_rates', 'solve_full_me', 'solve_reduced_me',
                'get_leak_coefficient', 'get_leak_branching_ratios', 'solve_ss_network',
                'update', 'update_configurations', 'add_path_reaction', 'add_net_reaction')


@pytest.mark.parametrize('entry', ENTRY_POINTS)
@pytest.mark.parametrize('slot', ('path_reactions', 'net_reactions'))
def test_current_reactions_checked_before_each_entry(entry, slot):
    network = pdep.PDepNetwork()
    reaction = helpers.cell_reaction(helpers.CENSUS_CELLS[0], True)
    admission = entry in ('add_path_reaction', 'add_net_reaction')
    if not admission:
        getattr(network, slot).append(reaction)
    before = dict(network.__dict__)
    with pytest.raises(NetworkError, match=reaction.library):
        ordinary = helpers.cell_reaction(helpers.CENSUS_CELLS[0], False)
        args = {
            'initialize': (300, 1000, 1e4, 1e6),
            'calculate_rate_coefficients': ([300], [1e5], 'modified strong collision'),
            'set_conditions': (300, 1e5), 'calculate_microcanonical_rates': (),
            'solve_full_me': ([], []), 'solve_reduced_me': ([], []),
            'get_leak_coefficient': (300, 1e5), 'get_leak_branching_ratios': (300, 1e5),
            'solve_ss_network': (300, 1e5), 'update': (None, None),
            'update_configurations': (None,), 'add_path_reaction': (reaction,),
            'add_net_reaction': (reaction,),
        }
        getattr(network, entry)(*args[entry])
    assert network.__dict__ == before


def test_no_per_attribute_guard():
    assert '__getattribute__' not in vars(pdep.PDepNetwork)


def test_unknown_reaction_refuses_without_reading_properties():
    class Proxy:
        @property
        def electrons(self):
            raise AssertionError('unknown reaction state must not be read')
    with pytest.raises(NetworkError, match='ProxyLibrary'):
        pdep._check_electron_channel_routing(Proxy(), library='ProxyLibrary')


def test_numeric_quantity_subclass_refuses():
    from rmgpy.quantity import ScalarQuantity
    class Quantity(ScalarQuantity):
        pass
    reaction = helpers.cell_reaction(helpers.CENSUS_CELLS[0], False)
    reaction.kinetics._A = Quantity(1, 's^-1')
    with pytest.raises(NetworkError, match=reaction.library):
        pdep.PDepNetwork().add_path_reaction(reaction)


def test_unknown_network_collection_refuses_without_iteration():
    class Reactions(list):
        def __iter__(self):
            raise AssertionError('unknown collection must not be inspected')
    network = pdep.PDepNetwork()
    network.path_reactions = Reactions()
    with pytest.raises(NetworkError, match='unsupported reaction-list type'):
        network.initialize(300, 1000, 1e4, 1e6)


def test_unknown_network_cannot_override_the_entry_check():
    class Network(pdep.PDepNetwork):
        def _check_reactions(self):
            raise AssertionError('the trusted entry must reject this type')
    with pytest.raises(NetworkError, match='unsupported network type'):
        Network().get_leak_coefficient(300, 1e5)


def test_unknown_participant_diagnostic_does_not_execute_str():
    class Participant:
        def __str__(self):
            raise AssertionError('unknown participant must not be formatted')
    reaction = helpers.cell_reaction(helpers.CENSUS_CELLS[0], False)
    reaction.reactants = [Participant()]
    with pytest.raises(NetworkError, match=reaction.library):
        pdep.PDepNetwork().add_path_reaction(reaction)


def test_numeric_value_is_not_a_material_participant():
    reaction = helpers.cell_reaction(helpers.CENSUS_CELLS[0], False)
    reaction.reactants = [0]
    with pytest.raises(NetworkError, match=reaction.library):
        pdep.PDepNetwork().add_path_reaction(reaction)


def test_unknown_type_diagnostic_does_not_format_class_metadata():
    class Metadata:
        def __str__(self):
            raise AssertionError('unknown class metadata must not be formatted')
    class Rate(kinetics.Arrhenius):
        pass
    Rate.__module__ = Metadata()
    reaction = helpers.cell_reaction(helpers.CENSUS_CELLS[0], False)
    reaction.kinetics = Rate(A=(1, 's^-1'))
    with pytest.raises(NetworkError, match=reaction.library):
        pdep.PDepNetwork().add_path_reaction(reaction)


def test_unknown_type_diagnostic_does_not_format_name_subclass():
    class Name(str):
        def __format__(self, spec):
            raise AssertionError('unsupported type name must not be formatted')
    class Rate(kinetics.Arrhenius):
        pass
    Rate.__name__ = Name('Rate')
    reaction = helpers.cell_reaction(helpers.CENSUS_CELLS[0], False)
    reaction.kinetics = Rate(A=(1, 's^-1'))
    with pytest.raises(NetworkError, match=reaction.library):
        pdep.PDepNetwork().add_path_reaction(reaction)


@pytest.mark.parametrize('cls', pdep._reaction_types(), ids=lambda cls: cls.__name__)
def test_shipped_reaction_ordinary_network_control(cls):
    reaction = cls(reactants=[Species().from_smiles('C')],
                   products=[Species().from_smiles('[CH3]')],
                   kinetics=kinetics.Arrhenius(A=(1, 's^-1')))
    network = pdep.PDepNetwork(source=reaction.reactants[:])
    network.add_path_reaction(reaction)
    assert network.path_reactions == [reaction]
    assert not pdep._has_electron_participant(reaction)


@pytest.mark.parametrize('form', ('empty', 'single', 'multiple'))
def test_template_atom_labels_are_aliases_of_checked_material(form):
    from rmgpy.data.kinetics.family import TemplateReaction
    reaction = TemplateReaction(reactants=[Species().from_smiles('CC')],
                                products=[Species().from_smiles('[CH2]C')],
                                kinetics=kinetics.Arrhenius(A=(1, 's^-1')))
    if form != 'empty':
        atoms = reaction.reactants[0].molecule[0].atoms
        reaction.labeled_atoms['reactants']['*1'] = atoms[0] if form == 'single' else atoms[:2]
        reaction.labeled_atoms['products']['*2'] = reaction.products[0].molecule[0].atoms[0]
    network = pdep.PDepNetwork(source=reaction.reactants[:])
    network.add_path_reaction(reaction)
    assert network.path_reactions == [reaction]


@pytest.mark.parametrize('form', ('foreign', 'detached_electron', 'wrong_side'))
def test_template_atom_labels_refuse_unsupported_or_detached_material(form):
    from rmgpy.data.kinetics.family import TemplateReaction
    reaction = TemplateReaction(reactants=[Species().from_smiles('C')],
                                products=[Species().from_smiles('[CH3]')],
                                kinetics=kinetics.Arrhenius(A=(1, 's^-1')),
                                family='LabelControl')
    if form == 'foreign':
        class Foreign:
            def __str__(self):
                raise AssertionError('unknown label reference must not be formatted')
        atom = Foreign()
    elif form == 'detached_electron':
        atom = Species().from_adjacency_list('1 e u0 p0 c-1').molecule[0].atoms[0]
    else:
        atom = reaction.products[0].molecule[0].atoms[0]
    reaction.labeled_atoms['reactants']['*1'] = atom
    network = pdep.PDepNetwork(source=reaction.reactants[:])
    with pytest.raises(NetworkError, match='LabelControl'):
        network.add_path_reaction(reaction)
    assert network.path_reactions == []


@pytest.mark.parametrize('form,native', [
    (form, native) for form in ('estimator', 'rank', 'transition_state', 'conformer', 'mode', 'tunneling')
    for native in (False, True) if form != 'estimator' or not native
])
def test_unknown_nested_reaction_record_refuses(form, native):
    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.species import TransitionState
    from rmgpy.statmech import Conformer, HarmonicOscillator
    class Foreign:
        electron = Species().from_adjacency_list('1 e u0 p0 c-1')
    from rmgpy.reaction import Reaction
    cls = Reaction if native else TemplateReaction
    options = {} if native else {'family': 'RecordControl'}
    reaction = cls(reactants=[Species().from_smiles('C')],
                   products=[Species().from_smiles('[CH3]')],
                   kinetics=kinetics.Arrhenius(A=(1, 's^-1')), **options)
    if form == 'estimator':
        reaction.estimator = Foreign()
    elif form == 'rank':
        reaction.rank = Foreign()
    else:
        state = TransitionState(conformer=Conformer(), tunneling=kinetics.Wigner(frequency=(-1000, 'cm^-1')))
        if form == 'transition_state':
            class State(TransitionState, Foreign):
                pass
            state = State()
        elif form == 'conformer':
            class Geometry(Conformer, Foreign):
                pass
            state.conformer = Geometry()
        elif form == 'mode':
            class Mode(HarmonicOscillator, Foreign):
                pass
            state.conformer.modes = [Mode()]
        else:
            class Tunneling(kinetics.Wigner, Foreign):
                pass
            state.tunneling = Tunneling(frequency=(-1000, 'cm^-1'))
        reaction.transition_state = state
    network = pdep.PDepNetwork(source=reaction.reactants[:])
    with pytest.raises(NetworkError, match='pressure-dependent networks'):
        network.add_path_reaction(reaction)
    assert network.path_reactions == []


def test_native_transition_state_statmech_network_control():
    from rmgpy.species import TransitionState
    from rmgpy.statmech import Conformer, HarmonicOscillator, IdealGasTranslation, LinearRotor
    from rmgpy.reaction import Reaction
    state = TransitionState(
        conformer=Conformer(E0=(10, 'kJ/mol'), modes=[
            IdealGasTranslation(mass=(16, 'amu')),
            LinearRotor(inertia=(1, 'amu*angstrom^2')),
            HarmonicOscillator(frequencies=([1000], 'cm^-1'))]),
        tunneling=kinetics.Wigner(frequency=(-1000, 'cm^-1')))
    reaction = Reaction(reactants=[Species().from_smiles('C')],
                        products=[Species().from_smiles('[CH3]')],
                        kinetics=kinetics.Arrhenius(A=(1, 's^-1')),
                        transition_state=state)
    network = pdep.PDepNetwork(source=reaction.reactants[:])
    network.add_path_reaction(reaction)
    assert network.path_reactions == [reaction]
    assert not pdep._has_electron_participant(reaction)



@pytest.mark.parametrize('cls', (pdep._reaction_types()[0], pdep._reaction_types()[2]),
                         ids=lambda cls: cls.__name__)
def test_unknown_record_metaclass_does_not_execute_hash_or_equality(cls):
    class Metadata(type):
        def __hash__(self):
            raise AssertionError('unsupported class must not be hashed')
        def __eq__(self, other):
            raise AssertionError('unsupported class must not be compared')
    class Foreign(metaclass=Metadata):
        pass
    reaction = cls(reactants=[Species().from_smiles('C')],
                   products=[Species().from_smiles('[CH3]')],
                   kinetics=kinetics.Arrhenius(A=(1, 's^-1')))
    reaction.rank = Foreign()
    with pytest.raises(NetworkError, match='pressure-dependent networks'):
        pdep.PDepNetwork().add_path_reaction(reaction)


@pytest.mark.parametrize('location', ('arrhenius_solute', 'multi_arrhenius_child'))
def test_eedf_prescan_refuses_unknown_kinetics_child_without_metaclass_callbacks(location):
    from rmgpy.reaction import Reaction

    calls = []

    class Metadata(type):
        def __hash__(self):
            calls.append('hash')
            return type.__hash__(self)

        def __eq__(self, other):
            calls.append('eq')
            return type.__eq__(self, other)

    class Foreign(metaclass=Metadata):
        pass

    rate = kinetics.Arrhenius(A=(1, 's^-1'))
    if location == 'arrhenius_solute':
        rate.solute = Foreign()
    else:
        rate = kinetics.MultiArrhenius(arrhenius=[rate])
        rate.arrhenius.append(Foreign())
    reaction = Reaction(
        reactants=[Species().from_smiles('C')],
        products=[Species().from_smiles('[CH3]')], kinetics=rate)

    with pytest.raises(NetworkError, match='pressure-dependent networks'):
        pdep.PDepNetwork().add_path_reaction(reaction)
    assert calls == []



def test_fast_physical_schema_covers_every_declared_object_slot():
    _assert_native_slot_inventory()


@pytest.mark.parametrize('field', ('value_si', 'uncertainty_si', 'array_subclass'))
@pytest.mark.parametrize('electron', (False, True))
def test_native_quantity_refuses_object_arrays_and_array_subclasses(field, electron):
    import numpy as np
    from rmgpy.reaction import Reaction
    reference = (Species().from_adjacency_list('1 e u0 p0 c-1') if electron
                 else Species().from_smiles('N#N'))
    rate = kinetics.KineticsData(Tdata=([300], 'K'), kdata=([1], 's^-1'))
    if field == 'array_subclass':
        class Array(np.ndarray):
            carrier = reference
            @property
            def dtype(self):
                raise AssertionError('unknown array must not be inspected')
        rate._Tdata.value_si = np.array([300.]).view(Array)
    else:
        setattr(rate._Tdata, field, np.array([reference], dtype=object))
    reaction = Reaction(reactants=[Species().from_smiles('C')],
                        products=[Species().from_smiles('[CH3]')], kinetics=rate)
    network = pdep.PDepNetwork(source=reaction.reactants[:])
    with pytest.raises(NetworkError, match='pressure-dependent networks'):
        network.add_path_reaction(reaction)
    assert network.path_reactions == []


@pytest.mark.parametrize('reaction_class', (pdep._reaction_types()[0], pdep._reaction_types()[1]))
@pytest.mark.parametrize('electron', (False, True))
@pytest.mark.parametrize('location', ('conformer', 'props', 'molecule_props', 'atom_coords',
                                     'thermo', 'transport', 'mass_transfer'))
def test_native_species_physical_graph_refuses_object_arrays(reaction_class, electron, location):
    import numpy as np
    from rmgpy.statmech import Conformer, HinderedRotor
    from rmgpy.thermo import ThermoData
    from rmgpy.transport import TransportData
    from rmgpy.data.vaporLiquidMassTransfer import HenryLawConstantData
    a = Species(label='ethane').from_smiles('CC')
    b = Species(label='methyl').from_smiles('[CH3]')
    reference = Species().from_adjacency_list('1 e u0 p0 c-1') if electron else 1.0
    array = np.array([reference], dtype=object)
    if location == 'conformer':
        rotor = HinderedRotor(inertia=(1, 'amu*angstrom^2'), barrier=(1, 'kJ/mol'))
        rotor.energies = array
        a.conformer = Conformer(E0=(0, 'kJ/mol'), modes=[rotor])
    elif location == 'props':
        a.props['carrier'] = array
    elif location == 'molecule_props':
        a.molecule[0].props['carrier'] = array
    elif location == 'atom_coords':
        a.molecule[0].atoms[0].coords = array
    elif location == 'thermo':
        a.thermo = ThermoData(Tdata=([300], 'K'), Cpdata=([30], 'J/(mol*K)'))
        a.thermo._Tdata.value_si = array
    elif location == 'transport':
        a.transport_data = TransportData()
        a.transport_data.epsilon = array
    else:
        a.henry_law_constant_data = HenryLawConstantData(Ts=[300], kHs=[1])
        a.henry_law_constant_data.kHs = array
    reaction = reaction_class(reactants=[a], products=[b, b],
        kinetics=kinetics.Arrhenius(A=(1, 's^-1')), elementary_high_p=True)
    if reaction_class is pdep._reaction_types()[1]:
        reaction.library = 'NativeSpeciesCarrier'
    network = pdep.PDepNetwork(source=[a])
    with pytest.raises(NetworkError, match='pressure-dependent networks'):
        network.add_path_reaction(reaction)
    assert network.path_reactions == []


def test_select_energy_grains_checks_current_initialized_path():
    spec = importlib.util.spec_from_file_location('native_grain_fixture',
        Path(pdep.__file__).parents[2] / 'test/rmgpy/rmg/pdepTest.py')
    fixture_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixture_module)
    fixture = fixture_module.TestPdep()
    fixture.setup_class()
    network = fixture.pdepnetwork
    network.energy_correction = None
    network.products_cache = []
    network.initialize(500., 1500., 1e4, 1e6, maximum_grain_size=2000., minimum_grain_count=100)
    network.path_reactions[0].specific_collider = Species().from_adjacency_list('1 e u0 p0 c-1')
    with pytest.raises(NetworkError, match='pressure-dependent networks'):
        network.select_energy_grains(1000., grain_count=100)


def _declared_object_slots(cls):
    """Inventory all Cython object slots, including private and inherited slots."""
    import ast
    import re
    root = Path(pdep.__file__).parents[2]
    primitive = {'bint', 'char', 'short', 'int', 'long', 'float', 'double', 'size_t'}
    slots = set()
    for base in cls.__mro__:
        if base is object:
            continue
        source = root.joinpath(*base.__module__.split('.'))
        pxd = source.with_suffix('.pxd')
        if pxd.exists():
            text = pxd.read_text()
            block = re.search(r'^cdef class ' + base.__name__ + r'(?:\([^\n]*\))?:\n((?:[ \t].*\n|\n)*)', text, re.MULTILINE)
            assert block is not None, (base.__name__, pxd)
            for declaration in re.findall(r'^    cdef\s+(?:public\s+|readonly\s+)?([^\n(]+)$', block.group(1), re.MULTILINE):
                parts = declaration.strip().split(None, 1)
                kind, names = parts if len(parts) == 2 else ('object', parts[0])
                if kind not in primitive:
                    slots.update(name.strip() for name in names.split(','))
        else:
            tree = ast.parse(source.with_suffix('.py').read_text())
            definition = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == base.__name__)
            for node in ast.walk(definition):
                if isinstance(node, (ast.Assign, ast.AnnAssign, ast.AugAssign)):
                    targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                    for target in targets:
                        if isinstance(target, ast.Attribute) and isinstance(target.value, ast.Name) and target.value.id == 'self':
                            slots.add(target.attr)
    return slots


_NATIVE_SLOT_EXCLUSIONS = {
    ('rmgpy.molecule.molecule', 'Molecule', '__weakref__'):
        'Weak-reference support is interpreter bookkeeping, not a material record slot.',
}


def _assert_native_slot_inventory():
    from rmgpy.quantity import ScalarQuantity, ArrayQuantity
    from rmgpy.molecule import Molecule
    schemas = dict(pdep._NATIVE_RECORD_FIELDS)
    for cls in pdep._KINETICS_TYPES:
        schemas[cls] = tuple(pdep._KINETICS_REFERENCE_FIELDS[cls]) + tuple(pdep._KINETICS_DATA_FIELDS[cls])
    for cls, fields in getattr(pdep, '_PYTHON_RECORD_FIELDS', {}).items():
        schemas[cls] = fields
    for cls, fields in getattr(pdep, '_NUMERIC_RECORD_FIELDS', {}).items():
        schemas[cls] = fields
    # Material roots and numeric leaves must participate in the same inventory.
    for cls in (Species, Molecule, ScalarQuantity, ArrayQuantity):
        schemas.setdefault(cls, ())
    aliases = getattr(pdep, '_NUMERIC_SLOT_ALIASES', {})
    for cls, fields in schemas.items():
        declared = {aliases.get(name, name) for name in _declared_object_slots(cls)}
        excluded = {slot for (module, name, slot), reason in _NATIVE_SLOT_EXCLUSIONS.items()
                    if (module, name) == (cls.__module__, cls.__name__)}
        assert excluded <= declared
        assert all(len(reason) > 30 for key, reason in _NATIVE_SLOT_EXCLUSIONS.items()
                   if key[:2] == (cls.__module__, cls.__name__))
        missing = declared - set(fields) - excluded
        assert not missing, f'Missing native record slots for {cls.__name__}: {sorted(missing)}'
    for cls, fields in schemas.items():
        if (cls not in pdep._NUMERIC_RECORD_FIELDS and cls not in pdep._PYTHON_RECORD_FIELDS
                and any(not hasattr(cls, name) for name in fields)):
            assert cls in pdep._NATIVE_RECORD_READERS, f'Missing native slot reader: {cls.__name__}'
    import ast
    source = Path(pdep.__file__).parents[2]
    for relative, name in (('rmgpy/rmg/pdep.py', '_is_numeric_data'),
                           ('rmgpy/reaction.py', '_is_native_channel_number')):
        tree = ast.parse((source / relative).read_text())
        definition = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
        reads = {aliases.get(n.attr, n.attr) for n in ast.walk(definition) if isinstance(n, ast.Attribute)}
        for cls, fields in pdep._NUMERIC_RECORD_FIELDS.items():
            assert set(fields) <= reads, (name, cls.__name__, set(fields) - reads)
    # Private slots use named compiled readers. Audit their actual tuple reads,
    # so adding a declaration to the schema alone cannot conceal an omitted read.
    import ast
    tree = ast.parse((Path(pdep.__file__).parents[2] / 'rmgpy/reaction.py').read_text())
    for cls, name in getattr(pdep, '_FAST_TYPED_RECORD_CHECKERS', {}).items():
        definition = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)
        reads = {n.attr for n in ast.walk(definition) if isinstance(n, ast.Attribute)
                 and isinstance(n.value, ast.Name) and n.value.id == 'atomtype'}
        assert set(schemas[cls]) == reads, (name, set(schemas[cls]) - reads)
    for cls, reader in getattr(pdep, '_NATIVE_RECORD_READERS', {}).items():
        definition = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == reader.__name__)
        expression = next(n.value for n in ast.walk(definition) if isinstance(n, ast.Return))
        assert isinstance(expression, ast.Tuple), reader.__name__
        actual = [n.attr for n in expression.elts if isinstance(n, ast.Attribute)]
        assert tuple(actual) == tuple(schemas[cls]), (cls.__name__, actual, schemas[cls])


def test_native_record_slot_inventory_is_complete():
    _assert_native_slot_inventory()


def test_native_record_slot_inventory_detects_removed_slot(monkeypatch):
    from rmgpy.statmech import Conformer
    _assert_native_slot_inventory()
    monkeypatch.setitem(pdep._NATIVE_RECORD_FIELDS, Conformer,
                        tuple(f for f in pdep._NATIVE_RECORD_FIELDS[Conformer] if f != 'modes'))
    with pytest.raises(AssertionError, match='Conformer.*modes'):
        _assert_native_slot_inventory()


def _assert_network_entry_inventory():
    exclusions = getattr(pdep, '_NETWORK_ACCESSOR_EXCLUSIONS', {})
    methods = {name: method for name, method in inspect.getmembers(pdep.PDepNetwork, inspect.isroutine)
               if not name.startswith('_')}
    assert set(exclusions) <= set(methods)
    for name, method in methods.items():
        if name in exclusions:
            assert len(exclusions[name]) > 30, name
        else:
            assert getattr(method, '_electron_channel_entry', False), f'Unguarded network entry: {name}'


def test_public_network_entry_inventory_is_complete():
    _assert_network_entry_inventory()


def test_public_network_entry_inventory_detects_removed_guard(monkeypatch):
    _assert_network_entry_inventory()
    monkeypatch.setattr(pdep.PDepNetwork, 'select_energy_grains',
                        pdep.PDepNetwork.select_energy_grains.__wrapped__)
    with pytest.raises(AssertionError, match='select_energy_grains'):
        _assert_network_entry_inventory()


_INVENTORIED_NETWORK_ENTRIES = tuple(name for name, method in
    inspect.getmembers(pdep.PDepNetwork, inspect.isroutine)
    if not name.startswith('_') and name not in pdep._NETWORK_ACCESSOR_EXCLUSIONS)


@pytest.mark.parametrize('entry', _INVENTORIED_NETWORK_ENTRIES)
@pytest.mark.parametrize('slot', ('path_reactions', 'net_reactions'))
def test_every_discovered_network_entry_checks_before_computation(entry, slot):
    network = pdep.PDepNetwork()
    reaction = helpers.cell_reaction(helpers.CENSUS_CELLS[0], True)
    admission = entry in ('add_path_reaction', 'add_net_reaction')
    if not admission:
        getattr(network, slot).append(reaction)
    before = dict(network.__dict__)
    function = getattr(network, entry)
    signature = inspect.signature(function)
    args, kwargs = [], {}
    for parameter in signature.parameters.values():
        if parameter.default is not inspect.Parameter.empty:
            continue
        if parameter.kind in (parameter.POSITIONAL_ONLY, parameter.POSITIONAL_OR_KEYWORD):
            args.append(reaction if admission else None)
        elif parameter.kind is parameter.KEYWORD_ONLY:
            kwargs[parameter.name] = None
    with pytest.raises(NetworkError, match=reaction.library):
        function(*args, **kwargs)
    assert network.__dict__ == before


@pytest.mark.parametrize('field', ('units', 'uncertainty_type'))
def test_numeric_leaf_metadata_must_be_exact_native_strings(field):
    from rmgpy.quantity import ScalarQuantity
    class Carrier(str):
        material = Species().from_adjacency_list('1 e u0 p0 c-1')
    quantity = ScalarQuantity(1, 's^-1')
    with pytest.raises(TypeError, match='Expected str'):
        setattr(quantity, field, Carrier('s^-1' if field == 'units' else '+|-'))
    reaction = helpers.cell_reaction(helpers.CENSUS_CELLS[0], False)
    reaction.kinetics._A = quantity
    network = pdep.PDepNetwork()
    network.add_path_reaction(reaction)
    assert network.path_reactions == [reaction]


@pytest.mark.parametrize('field', ('generic', 'single'))
@pytest.mark.parametrize('electron', (False, True))
def test_typed_atomtype_walk_refuses_arrays_in_each_list_kind(field, electron):
    import numpy as np
    from rmgpy.molecule.atomtype import AtomType
    reference = Species().from_adjacency_list('1 e u0 p0 c-1') if electron else 1.0
    atomtype = AtomType(label='NativeDescriptor')
    setattr(atomtype, field, [np.array([reference], dtype=object)])
    reaction = helpers.cell_reaction(helpers.CENSUS_CELLS[0], False)
    reaction.reactants[0].molecule[0].atoms[0].atomtype = atomtype
    with pytest.raises(NetworkError, match=reaction.library):
        pdep.PDepNetwork().add_path_reaction(reaction)


def test_private_record_slot_inventory_detects_missing_reader(monkeypatch):
    _assert_native_slot_inventory()
    monkeypatch.delitem(pdep._NATIVE_RECORD_READERS, Species)
    with pytest.raises(AssertionError, match='Missing native slot reader: Species'):
        _assert_native_slot_inventory()


def _module_network_entries():
    """Discover module-level network parameters across both production packages."""
    import re
    import runpy
    import io
    import tokenize
    root = Path(pdep.__file__).parents[2]
    census = runpy.run_path(str(root / 'scripts/state_data_census.py'))
    entries = []
    for key, site in census['collect_sites'](root).items():
        relative, name = key.split(':', 1)
        if '.' in name or '<' in name:
            continue
        lines = (root / relative).read_text().splitlines()
        remaining = lines[site['line'] - 1:]
        depth = 0
        started = False
        for token in tokenize.generate_tokens(io.StringIO('\n'.join(remaining)).readline):
            if token.type == tokenize.OP and token.string == '(':
                depth += 1
                started = True
            elif token.type == tokenize.OP and token.string == ')':
                depth -= 1
                if started and depth == 0:
                    header = '\n'.join(remaining[:token.end[0] - 1] + [remaining[token.end[0] - 1][:token.end[1]]])
                    break
        else:
            raise AssertionError(f'Unparsed function signature: {key}')
        if re.search(r'(?:\(|,)\s*(?:[\w.]+\s+)?network\b', header):
            module = relative.rsplit('.', 1)[0].replace('/', '.')
            entries.append((module, name))
    return tuple(sorted(entries))


_MODULE_NETWORK_EXCLUSIONS = {
    ('rmgpy.rmg.pdep', '_check_network_reactions'):
        'The validation boundary itself checks the original current channels; it performs no numerical computation or registration.',
    ('rmgpy.thermo.state', 'require_network_thermo_allowed'):
        'This thermo validation boundary performs no electron-routing numerical computation or registration.',
}


def _assert_module_network_inventory():
    entries = _module_network_entries()
    assert set(_MODULE_NETWORK_EXCLUSIONS) <= set(entries)
    for module_name, name in entries:
        if (module_name, name) in _MODULE_NETWORK_EXCLUSIONS:
            assert len(_MODULE_NETWORK_EXCLUSIONS[module_name, name]) > 30
            continue
        function = getattr(importlib.import_module(module_name), name)
        for network_class in (pdep.PDepNetwork, pdep.rmgpy.pdep.network.Network):
            for slot in ('path_reactions', 'net_reactions'):
                network = network_class()
                reaction = helpers.cell_reaction(helpers.CENSUS_CELLS[0], True)
                getattr(network, slot).append(reaction)
                parameters = inspect.signature(function).parameters.values()
                args = [network if parameter.name == 'network' else None
                        for parameter in parameters
                        if parameter.default is inspect.Parameter.empty]
                try:
                    function(*args)
                except NetworkError as error:
                    assert reaction.library in str(error), (module_name, name, str(error))
                except Exception as error:
                    raise AssertionError(f'Unguarded module entry: {module_name}.{name}') from error
                else:
                    raise AssertionError(f'Unguarded module entry: {module_name}.{name}')


def test_module_network_inventory_is_complete():
    _assert_module_network_inventory()


def test_module_network_inventory_detects_removed_guard(monkeypatch):
    import rmgpy.pdep.me as me
    _assert_module_network_inventory()
    monkeypatch.setattr(me, 'generate_full_me_matrix', lambda network: None)
    with pytest.raises(AssertionError, match='generate_full_me_matrix'):
        _assert_module_network_inventory()


@pytest.mark.parametrize('size', (4, 8, 16))
@pytest.mark.parametrize('reaction_class', pdep._reaction_types())
def test_admission_record_reads_are_linear(size, reaction_class, monkeypatch):
    original = pdep._NATIVE_RECORD_READERS[Species]
    reads, channels = [], []
    check = pdep._check_electron_channel_routing
    def channel_check(reaction, *args, **kwargs):
        channels.append(id(reaction))
        return check(reaction, *args, **kwargs)
    monkeypatch.setattr(pdep, '_check_electron_channel_routing', channel_check)
    def reader(value):
        reads.append(id(value))
        return original(value)
    monkeypatch.setitem(pdep._NATIVE_RECORD_READERS, Species, reader)
    a, b = Species(label='methane').from_smiles('C'), Species(label='methyl').from_smiles('[CH3]')
    network = pdep.PDepNetwork(source=[a])
    for index in range(size):
        reaction = reaction_class(reactants=[a], products=[b], kinetics=kinetics.Arrhenius(A=(1, 's^-1')))
        network.add_net_reaction(reaction)
    assert len(reads) == 2 * size
    assert len(channels) == size
    reads.clear()
    channels.clear()
    network._check_reactions()
    assert sorted(reads) == sorted((id(a), id(b)))
    assert len(channels) == size
    # A new call must inspect a mutation despite shared identity in the prior call.
    import numpy as np
    a.props['carrier'] = np.array([1.0], dtype=object)
    with pytest.raises(NetworkError, match='pressure-dependent networks'):
        network._check_reactions()


@pytest.mark.parametrize('size', (4, 8, 16))
@pytest.mark.parametrize('explicit_network', (False, True))
def test_model_admission_record_reads_are_linear(size, explicit_network, monkeypatch):
    from rmgpy.reaction import Reaction
    import rmgpy.rmg.model as model_module

    original = pdep._NATIVE_RECORD_READERS[Species]
    reads, channels = [], []
    check = pdep._check_electron_channel_routing

    def channel_check(reaction, *args, **kwargs):
        channels.append(id(reaction))
        return check(reaction, *args, **kwargs)

    def reader(value):
        reads.append(id(value))
        return original(value)

    monkeypatch.setattr(pdep, '_check_electron_channel_routing', channel_check)
    monkeypatch.setattr(model_module, '_check_electron_channel_routing', channel_check)
    monkeypatch.setitem(pdep._NATIVE_RECORD_READERS, Species, reader)

    source = Species(label='n-C34H70').from_smiles('C' * 34)

    def alkyl(length):
        smiles = '[CH3]' if length == 1 else '[CH2]' + 'C' * (length - 1)
        return Species(label=f'C{length}H{2 * length + 1}').from_smiles(smiles)

    reactions = [
        Reaction(
            reactants=[source],
            products=[alkyl(index), alkyl(34 - index)],
            kinetics=kinetics.Arrhenius(A=(1, 's^-1')),
        )
        for index in range(1, size + 1)
    ]
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    monkeypatch.setattr(model, 'make_new_reaction',
                        lambda reaction, **kwargs: (reaction, True))
    network = pdep.PDepNetwork(source=[source])
    for reaction in reactions:
        model.process_new_reactions(
            [reaction],
            source,
            pdep_network=network if explicit_network else None,
            generate_thermo=False,
            generate_kinetics=False,
        )
    admitted = network if explicit_network else model.network_list[0]

    assert len(admitted.path_reactions) == size
    assert len(reads) == 9 * size
    assert len(channels) == 2 * size


def test_fast_fallback_does_not_commit_unchecked_records():
    from rmgpy.reaction import Reaction
    import numpy as np
    a, b = Species().from_smiles('C'), Species().from_smiles('[CH3]')
    # Shipped kinetics in physical metadata uses the general walker after a
    # native fast attempt. Its other unprocessed roots must still be checked.
    a.props['nested'] = kinetics.Arrhenius(A=(1, 's^-1'))
    b.props['carrier'] = np.array([1.0], dtype=object)
    network = pdep.PDepNetwork()
    network.net_reactions.append(Reaction(reactants=[a], products=[b], kinetics=kinetics.Arrhenius(A=(1, 's^-1'))))
    with pytest.raises(NetworkError, match='pressure-dependent networks'):
        network._check_reactions()


@pytest.mark.parametrize('reaction_class', pdep._reaction_types())
def test_compiled_material_walk_preserves_reference_order_and_shared_molecule(reaction_class):
    from rmgpy.molecule import Molecule
    molecule = Molecule().from_smiles('C')
    a, b = Species(label='left', molecule=[molecule]), Species(label='right', molecule=[molecule])
    reaction = reaction_class(reactants=[a], products=[b], kinetics=kinetics.Arrhenius(A=(1, 's^-1')))
    assert list(pdep._reaction_species_references(reaction)) == [b, molecule, a]
