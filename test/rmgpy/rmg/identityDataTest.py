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

"""Plain identity data cannot introduce material or provenance object graphs."""

import io
import pickle
from types import SimpleNamespace

import pytest

from rmgpy.reaction import Reaction
from rmgpy.rmg import pdep
from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.common import save_entry
from rmgpy.data.kinetics.library import LibraryReaction, KineticsLibrary
from rmgpy.kinetics import Arrhenius, ThirdBody, EEDFChannel
from rmgpy.species import Species
from rmgpy.exceptions import DatabaseError, NetworkError

callbacks = []


class StringCarrier(str):
    def __new__(cls, reference):
        obj = super().__new__(cls, 'Control')
        obj.reference = reference
        return obj
    def __reduce__(self):
        callbacks.append('string_reduce')
        return (StringCarrier, (self.reference,))


class ReducerCarrier:
    def __init__(self, reference):
        self.reference = reference
    def __reduce__(self):
        callbacks.append('reduce')
        return (tuple, ([self.reference],))


def species(electron=False):
    return (Species(label='hidden_e').from_adjacency_list('1 e u0 p0 c-1')
            if electron else Species(label='hidden_N2').from_smiles('N#N'))


def ordinary():
    source = Species(label='H2').from_smiles('[H][H]')
    atom = Species(label='H').from_smiles('[H]')
    r = LibraryReaction(reactants=[source], products=[atom, atom],
        kinetics=Arrhenius(A=(1., 's^-1')), reversible=False,
        library='Control', source_identity=('Control', 1), elementary_high_p=True)
    r.library_source = '/external/Control'
    return r


def payload(kind):
    if kind == 'Species': return species()
    if kind == 'Molecule': return species().molecule[0]
    if kind == 'electron_Species': return species(True)
    if kind == 'electron_Molecule': return species(True).molecule[0]
    if kind == 'Reaction':
        return Reaction(reactants=[species()], products=[],
            kinetics=Arrhenius(A=(1., 's^-1')), reversible=False)
    if kind == 'electron_Reaction':
        return Reaction(reactants=[species(True)], products=[],
            kinetics=Arrhenius(A=(1., 's^-1')), reversible=False)
    if kind == 'EEDFChannel': return EEDFChannel('Ar -> Ar*', 'argon-lxcat-v1', 'ine')
    if kind == 'EEDF_Reaction':
        return Reaction(reactants=[Species(label='Ar').from_smiles('[Ar]')],
            products=[Species(label='Ar*').from_smiles('[Ar]')],
            kinetics=payload('EEDFChannel'), reversible=False)
    if kind in ('ThirdBody_Species', 'ThirdBody_electron'):
        return ThirdBody(arrheniusLow=Arrhenius(A=(1., 'm^3/(mol*s)')),
            efficiencies={species(kind.endswith('electron')).molecule[0]: 2.})
    if kind == 'LibraryReaction_entry_electron':
        r = ordinary()
        r.entry = Entry(item=species(True))
        return r
    if kind == 'Entry_electron': return Entry(item=species(True))
    raise AssertionError(kind)


def wrap(kind, reference):
    if kind == 'direct': return reference
    if kind == 'tuple': return ('owner', ('nested', (reference,)))
    if kind == 'dict': return {'nested': ('owner', {'hidden': reference})}
    if kind == 'str_subclass': return StringCarrier(reference)
    if kind == 'reducer': return ReducerCarrier(reference)
    raise AssertionError(kind)


FIELDS = ('source_identity', 'library_source')
KINDS = ('Species', 'Molecule', 'electron_Species', 'electron_Molecule',
         'Reaction', 'electron_Reaction', 'EEDFChannel', 'EEDF_Reaction',
         'ThirdBody_Species', 'ThirdBody_electron',
         'LibraryReaction_entry_electron', 'Entry_electron')
CARRIERS = ('direct', 'tuple', 'dict', 'str_subclass', 'reducer')


@pytest.fixture(autouse=True)
def no_callbacks():
    callbacks.clear()
    yield
    assert not callbacks


@pytest.mark.parametrize('field', FIELDS)
@pytest.mark.parametrize('kind', KINDS)
@pytest.mark.parametrize('carrier', CARRIERS)
@pytest.mark.parametrize('operation', ('guard', 'register', 'network', 'revalidate'))
def test_identity_probe_refuses_object_graphs(field, kind, carrier, operation):
    """The reviewer's 480 cases, including bypassed setters and later mutation."""
    reaction = ordinary()
    value = wrap(carrier, payload(kind))
    original = getattr(reaction, field)
    with pytest.raises(ValueError, match=field + '.*plain data'):
        setattr(reaction, field, value)
    assert getattr(reaction, field) == original

    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    before = pickle.dumps(model.reaction_dict)
    network = pdep.PDepNetwork(source=reaction.reactants[:])
    if operation == 'revalidate':
        network.add_path_reaction(reaction)
    # Old/restored state or direct dictionary mutation must not bypass routing.
    reaction.__dict__[field] = value
    with pytest.raises(NetworkError, match='reaction[.]' + field + '.*plain data'):
        if operation == 'guard':
            pdep._check_electron_channel_routing(reaction)
        elif operation == 'register':
            model.register_reaction(reaction)
        elif operation == 'network':
            network.add_path_reaction(reaction)
        else:
            pdep._check_network_reactions(network)
    assert pickle.dumps(model.reaction_dict) == before
    assert network.path_reactions == ([reaction] if operation == 'revalidate' else [])


@pytest.mark.parametrize('field', FIELDS)
@pytest.mark.parametrize('kind', KINDS)
@pytest.mark.parametrize('carrier', CARRIERS)
@pytest.mark.parametrize('transport', ('read', 'copy', 'pickle', 'restore'))
def test_identity_transport_refuses_corrupted_state(field, kind, carrier, transport):
    reaction = ordinary()
    value = wrap(carrier, payload(kind))
    if transport == 'restore':
        before = dict(reaction.__dict__)
        with pytest.raises(ValueError, match=field + '.*plain data'):
            reaction.__setstate__({'library': 'changed', field: value})
        assert reaction.__dict__ == before
    else:
        reaction.__dict__[field] = value
        with pytest.raises(ValueError, match=field + '.*plain data'):
            if transport == 'read':
                getattr(reaction, field)
            elif transport == 'copy':
                reaction.copy()
            else:
                pickle.dumps(reaction)


@pytest.mark.parametrize('kind', KINDS)
@pytest.mark.parametrize('carrier', CARRIERS)
def test_identity_constructor_refuses_objects(kind, carrier):
    with pytest.raises(ValueError, match='source_identity.*plain data'):
        LibraryReaction(source_identity=wrap(carrier, payload(kind)))


@pytest.mark.parametrize('field', FIELDS)
@pytest.mark.parametrize('kind', ('Reaction', 'EEDFChannel', 'LibraryReaction_entry_electron'))
def test_legacy_pickle_load_refuses_object_identity(field, kind):
    class LegacyState:
        def __reduce__(self):
            return LibraryReaction, (), {field: payload(kind)}

    # Produce the pre-validation state format without using the new reducer.
    data = pickle.dumps(LegacyState())
    with pytest.raises(ValueError, match=field + '.*plain data'):
        pickle.loads(data)


@pytest.mark.parametrize('field', FIELDS)
@pytest.mark.parametrize('value', (None, 'Control', 1, 1.5, True, (),
                                 ('Control', 1), ('plain', (None, 2.5, False))))
@pytest.mark.parametrize('transport', ('direct', 'copy', 'pickle', 'restore'))
def test_plain_identity_data_admitted(field, value, transport):
    reaction = ordinary()
    setattr(reaction, field, value)
    if transport == 'copy':
        reaction = reaction.copy()
    elif transport == 'pickle':
        reaction = pickle.loads(pickle.dumps(reaction))
    elif transport == 'restore':
        reaction.__setstate__({field: value})
    assert getattr(reaction, field) == value
    pdep._check_electron_channel_routing(reaction)
    model = CoreEdgeReactionModel()
    model.pressure_dependence = SimpleNamespace(maximum_atoms=None)
    model.register_reaction(reaction)
    network = pdep.PDepNetwork(source=reaction.reactants[:])
    network.add_path_reaction(reaction)
    pdep._check_network_reactions(network)
    assert network.path_reactions == [reaction]


@pytest.mark.parametrize('transport', ('direct', 'copy', 'pickle', 'restore'))
def test_missing_library_source_reads_as_none_after_transport(transport):
    reaction = ordinary()
    del reaction.__dict__['library_source']
    if transport == 'copy':
        reaction = reaction.copy()
    elif transport == 'pickle':
        reaction = pickle.loads(pickle.dumps(reaction))
    elif transport == 'restore':
        reaction.__setstate__({'source_identity': ('Control', 1)})
    # Public-descriptor serializers discover this property even without a key
    # in the old state dictionary. Reading it must not invent stored provenance.
    assert 'library_source' in dir(reaction)
    assert reaction.library_source is None
    assert 'library_source' not in reaction.__dict__
    pdep._check_electron_channel_routing(reaction)


@pytest.mark.parametrize('kind', ('electron_Reaction', 'EEDF_Reaction'))
def test_direct_electron_channels_still_refused(kind):
    reaction = payload(kind)
    network = pdep.PDepNetwork(source=reaction.reactants[:])
    with pytest.raises(NetworkError):
        network.add_path_reaction(reaction)
    assert not network.path_reactions


@pytest.mark.parametrize('field', FIELDS)
def test_identity_depth_bound_and_shared_tuples(field):
    reaction = ordinary()
    value = 'leaf'
    for _ in range(32):
        value = (value,)
    setattr(reaction, field, value)
    pdep._check_electron_channel_routing(reaction)
    # The shared value is valid alone but too deep inside either extra wrapper.
    for invalid in ((value,), (value, (value,))):
        with pytest.raises(ValueError, match='maximum tuple depth'):
            setattr(reaction, field, invalid)
        reaction.__dict__[field] = invalid
        with pytest.raises(NetworkError, match='maximum tuple depth'):
            pdep._check_electron_channel_routing(reaction)


@pytest.mark.parametrize('base', (str, int, float, tuple))
def test_identity_rejects_scalar_subclasses_without_callbacks(base):
    class Carrier(base):
        def __iter__(self):
            callbacks.append('iter')
            raise AssertionError

        def __reduce__(self):
            callbacks.append('reduce')
            raise AssertionError

    value = Carrier()
    value.reference = species(True)
    for field in FIELDS:
        with pytest.raises(ValueError, match='plain data'):
            setattr(ordinary(), field, value)


def test_identity_rejects_custom_metaclass_without_comparison():
    class Meta(type):
        def __eq__(self, other):
            callbacks.append('eq')
            raise AssertionError

        def __hash__(self):
            callbacks.append('hash')
            raise AssertionError

    class Carrier(metaclass=Meta):
        pass

    with pytest.raises(ValueError, match='plain data'):
        ordinary().source_identity = Carrier()


@pytest.mark.parametrize('keyword,value', (
    ('sourceIdentity', (StringCarrier(species()), 1)),
    ('sourceProvenance', (StringCarrier(species()), None)),
    ('sourceProvenance', ('Control', StringCarrier(species()))),
))
def test_seed_loader_refuses_string_subclasses(keyword, value):
    library = KineticsLibrary(label='seed', auto_generated=True)
    with pytest.raises(DatabaseError, match='Invalid source'):
        library.load_entry(index=1, label='H2 => H + H', kinetics=None, **{keyword: value})
    assert not library.entries


@pytest.mark.parametrize('field', ('source_identity', 'source_provenance'))
@pytest.mark.parametrize('operation', ('load', 'save'))
def test_seed_entry_readback_refuses_mutated_provenance(field, operation):
    library = KineticsLibrary(label='seed', auto_generated=True)
    library.load_entry(index=1, label='H2 => H + H', kinetics=Arrhenius(A=(1., 's^-1')),
                       sourceIdentity=('Control', 1), sourceProvenance=('Control', None))
    entry = library.entries[1]
    setattr(entry, field, ('Control', ReducerCarrier(species(True))))
    with pytest.raises(ValueError, match=field + '.*plain data'):
        if operation == 'load':
            library.get_library_reactions()
        else:
            save_entry(io.StringIO(), entry)
