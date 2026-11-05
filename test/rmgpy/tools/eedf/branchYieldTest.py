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

"""Sourced cascade vectors: fail closed before a collision can be branched."""

import copy
import hashlib
import math
import os
import runpy
from pathlib import Path
import subprocess
import sys

import pytest
import yaml

from rmgpy.tools.eedf.schema import FingerprintMismatch, SpecError


ID6 = {
    'Ar(1s5)': 0.7420494699646644,
    'Ar(1s4)': 0.21201413427561838,
    'Ar(1s3)': 0.038476639183352965,
    'Ar(1s2)': 0.00745975657636435,
}
PROCESS = 'e + Ar(1S0) -> e + Ar(4p[1/2]1), Excitation'


def database_cascade_channels(root):
    """Independent source regeneration supports both pre/post cascade databases."""
    generator = runpy.run_path(str(root / 'scripts/gen_argon_cascades.py'))
    import inspect
    args = [(root / 'input/eedf/argon_radiative_data.yaml').read_bytes()]
    if 'channel_map_raw' in inspect.signature(generator['generate']).parameters:
        args.append((root / 'input/eedf/argon_channel_map.yaml').read_bytes())
    channels = generator['generate'](*args)['channels']
    assert channels and any(row['status'] == 'complete_for_E1' for row in channels)
    return channels


def branch_reactions(snapshot, channel_id=6, collision_set='test-Ar'):
    """Actual state-resolved reactions; the source supplies the fractions only."""
    from rmgpy.data.base import Entry
    from rmgpy.data.kinetics.library import LibraryReaction
    from rmgpy.kinetics import EEDFChannel
    from rmgpy.species import Species
    row = next(row for row in snapshot['channels'] if row['id'] == channel_id)
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    argon = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    reactions = []
    # Independent literal state correspondence, including the database alias.
    for offset, (product, label, state, spin, radicals, pairs) in enumerate([
        ('Ar(1s5)', 'Ars', '3P2', 3, 2, 3),
        ('Ar(1s4)', 'Ar(1s4)', '3P1', 3, 2, 3),
        ('Ar(1s3)', 'Ar(1s3)', '3P0', 3, 2, 3),
        ('Ar(1s2)', 'Ar(1s2)', '1P1', 1, 0, 4),
    ]):
        excited = Species(label=label).from_adjacency_list(
            f'multiplicity {spin}\nelectronicstate {state}\n1 Ar u{radicals} p{pairs} c0')
        marker = EEDFChannel(process=row['process'], collision_set=collision_set,
                             branch_product=product, branch_yield=row['product_yields'][product],
                             yield_table_sha256=snapshot['yield_table_sha256'])
        entry = Entry(index=1000 + channel_id * 4 + offset, data=marker)
        reaction = LibraryReaction(reactants=[electron, argon], products=[electron, excited],
                                   reversible=False, kinetics=marker, library='PlasmaArgon', entry=entry,
                                   source_identity=('PlasmaArgon', entry.index))
        entry.item = reaction
        reactions.append(reaction)
    return reactions


def family_from_reactions(snapshot, reactions, channel_id=6):
    from rmgpy.tools.eedf.ownership import make_branch_family
    row = next(row for row in snapshot['channels'] if row['id'] == channel_id)
    return make_branch_family(snapshot, channel_id, row['process'], row['threshold_eV'], 'test-Ar',
                              [(r.library, r.entry.index, r, r.kinetics) for r in reactions])


@pytest.fixture
def branch_fixture(source_files):
    from rmgpy.tools.eedf.ownership import load_cascade_yields
    snapshot = load_cascade_yields(*source_files[:2])
    reactions = branch_reactions(snapshot)
    return snapshot, reactions, family_from_reactions(snapshot, reactions)


def test_family_freezes_sourced_vector_and_molecular_alias(branch_fixture):
    from rmgpy.tools.eedf.ownership import iter_reaction_owners, reaction_product_state
    snapshot, reactions, family = branch_fixture
    assert family['schema'] == 'eedf-branch-family/1'
    assert family['product_yields'] == ID6
    assert family['yield_table_sha256'] == snapshot['yield_table_sha256']
    assert reactions[0].products[1].label == 'Ars'
    assert reaction_product_state(reactions[0]) == 'Ar(1s5)'
    assert list(iter_reaction_owners(family)) == [
        {'library': r.library, 'index': r.entry.index, 'repr': repr(r.kinetics)} for r in reactions]
    snapshot['channels'][0]['product_yields']['Ar(1s5)'] = 0.
    assert family['product_yields'] == ID6


@pytest.mark.parametrize('problem', [
    'missing', 'extra', 'duplicate_product', 'duplicate_entry', 'product_state',
    'reactant_state', 'stoichiometry', 'reversible', 'collision_set', 'source', 'yield',
])
def test_family_refuses_invalid_reactions(branch_fixture, problem):
    from rmgpy.kinetics import EEDFChannel
    snapshot, reactions, _ = branch_fixture
    expected = {
        'missing': 'all four', 'extra': 'all four', 'duplicate_product': 'unique',
        'duplicate_entry': 'duplicate branch source', 'product_state': 'molecular product',
        'reactant_state': 'ground-state', 'stoichiometry': 'one electron',
        'reversible': 'irreversible', 'collision_set': 'sourced family',
        'source': 'sourced family', 'yield': 'sourced family',
    }[problem]
    if problem == 'missing':
        reactions.pop()
    elif problem == 'extra':
        reactions.append(copy.deepcopy(reactions[0]))
    elif problem == 'duplicate_product':
        from rmgpy.data.base import Entry
        duplicate = copy.deepcopy(reactions[0])
        duplicate.entry = Entry(index=reactions[-1].entry.index)
        reactions[-1] = duplicate
    elif problem == 'duplicate_entry':
        reactions[-1].entry.index = reactions[0].entry.index
    elif problem == 'product_state':
        reactions[0].products[1].molecule[0].electronic_state = '3P1'
    elif problem == 'reactant_state':
        reactions[0].reactants[1] = reactions[0].products[1]
    elif problem == 'stoichiometry':
        reactions[0].products.append(reactions[0].products[0])
    elif problem == 'reversible':
        reactions[0].reversible = True
    else:
        marker = reactions[0].kinetics
        fields = dict(process=marker.process, collision_set=marker.collision_set,
                      branch_product=marker.branch_product, branch_yield=marker.branch_yield,
                      yield_table_sha256=marker.yield_table_sha256)
        fields[{'collision_set': 'collision_set', 'source': 'yield_table_sha256',
                'yield': 'branch_yield'}[problem]] = {
                    'collision_set': 'other', 'source': '0' * 64, 'yield': .5}[problem]
        reactions[0].kinetics = EEDFChannel(**fields)
    with pytest.raises((SpecError, FingerprintMismatch), match=expected):
        family_from_reactions(snapshot, reactions)


def test_missing_zero_branch_refuses(source_files):
    from rmgpy.tools.eedf.ownership import load_cascade_yields
    snapshot = load_cascade_yields(*source_files[:2])
    snapshot['channels'][0]['product_yields'] = dict(zip(ID6, [1., 0., 0., 0.]))
    reactions = branch_reactions(snapshot)
    family = family_from_reactions(snapshot, reactions)
    assert family['product_yields']['Ar(1s2)'] == 0.
    with pytest.raises(SpecError, match='all four'):
        family_from_reactions(snapshot, reactions[:-1])


@pytest.mark.parametrize('problem', ['schema', 'policy', 'extra_field', 'marker_code', 'marker_sup',
                                     'sum', 'digest', 'threshold', 'multiplicity', 'charged_product'])
def test_family_envelope_and_species_refusals(branch_fixture, problem):
    from rmgpy.tools.eedf.ownership import validate_family
    snapshot, reactions, family = branch_fixture
    if problem in ('multiplicity', 'charged_product'):
        molecule = reactions[0].products[1].molecule[0]
        if problem == 'multiplicity':
            molecule.multiplicity = 1
        else:
            molecule.atoms[0].charge = 1
        with pytest.raises(SpecError, match='argon'):
            family_from_reactions(snapshot, reactions)
        return
    if problem == 'schema':
        family['schema'] = 'eedf-branch-family/2'
    elif problem == 'policy':
        family['sum_atol'] = .01
    elif problem == 'extra_field':
        family['unrecognized'] = True
    elif problem == 'marker_code':
        family['branches'][0]['repr'] = "__import__('os').getcwd()"
    elif problem == 'marker_sup':
        family['branches'][0]['repr'] = family['branches'][0]['repr'].replace("side='ine'", "side='sup'")
    elif problem == 'sum':
        family['product_yields']['Ar(1s5)'] += .001
    elif problem == 'digest':
        family['source_sha256'] = None
    else:
        family['threshold_eV'] = 12.
    with pytest.raises(SpecError):
        validate_family(family, process=PROCESS, threshold_eV=12.907)


def test_sum_preserving_source_mutation_refuses_stale_marker(branch_fixture):
    snapshot, reactions, _ = branch_fixture
    values = snapshot['channels'][0]['product_yields']
    values['Ar(1s5)'] += .001
    values['Ar(1s4)'] -= .001
    with pytest.raises(FingerprintMismatch, match='sourced family'):
        family_from_reactions(snapshot, reactions)


@pytest.fixture
def source_files(tmp_path):
    radiative = tmp_path / 'radiative.yaml'
    radiative.write_text('sourced: fixture\n')
    document = {
        'schema_version': 1, 'radiative_model': 'E1',
        'source_sha256': hashlib.sha256(radiative.read_bytes()).hexdigest(),
        'channels': [
            {'id': 6, 'process': PROCESS, 'threshold_eV': 12.907,
             'upper_nist_id': '018001.000006', 'paschen': '2p10',
             'status': 'complete_for_E1', 'product_yields': dict(ID6),
             'missing_transitions': [], 'transition_ids': ['upper->lower']},
            {'id': 11, 'process': 'unresolved collision', 'threshold_eV': 13.273,
             'status': 'unresolved', 'product_yields': None},
        ],
    }
    yields = tmp_path / 'yields.yaml'
    yields.write_text(yaml.safe_dump(document))
    return yields, radiative, document


def test_id6_literal_vector_and_source_identity(source_files):
    from rmgpy.tools.eedf.ownership import YIELD_SUM_ATOL, load_cascade_yields, require_cascade_yields
    yields, radiative, _ = source_files
    snapshot = load_cascade_yields(yields, radiative)
    row = require_cascade_yields(snapshot, 6, process=PROCESS, threshold_eV=12.907)
    assert row['product_yields'] == ID6
    assert snapshot['yield_table_sha256'] == hashlib.sha256(yields.read_bytes()).hexdigest()
    assert snapshot['source_sha256'] == hashlib.sha256(radiative.read_bytes()).hexdigest()
    assert abs(math.fsum(row['product_yields'].values()) - 1.) <= YIELD_SUM_ATOL


@pytest.mark.parametrize('bad', [None, -0.1, float('nan'), float('inf'), 10 ** 1000, 1.1, True, '0.5'])
def test_invalid_yield_refuses(bad):
    from rmgpy.tools.eedf.ownership import validate_yield_vector
    values = dict(ID6, **{'Ar(1s5)': bad})
    with pytest.raises(SpecError, match='yield'):
        validate_yield_vector(values)


def test_yield_sum_plus_1e_minus3_refuses():
    from rmgpy.tools.eedf.ownership import validate_yield_vector
    values = dict(ID6)
    values['Ar(1s5)'] += 1e-3
    with pytest.raises(SpecError, match='sum'):
        validate_yield_vector(values)


def test_yield_sum_minus_1e_minus3_refuses():
    from rmgpy.tools.eedf.ownership import validate_yield_vector
    values = dict(ID6)
    values['Ar(1s5)'] -= 1e-3
    with pytest.raises(SpecError, match='sum'):
        validate_yield_vector(values)


@pytest.mark.parametrize('offset, accepted', [(5e-13, True), (-5e-13, True), (2e-12, False), (-2e-12, False)])
def test_fixed_normalization_tolerance(offset, accepted):
    from rmgpy.tools.eedf.ownership import validate_yield_vector
    values = dict(ID6)
    values['Ar(1s5)'] += offset
    if accepted:
        assert validate_yield_vector(values) == values  # no renormalization
    else:
        with pytest.raises(SpecError, match='sum'):
            validate_yield_vector(values)


def test_explicit_zero_is_valid_but_missing_zero_is_not():
    from rmgpy.tools.eedf.ownership import validate_yield_vector
    values = dict(zip(ID6, [1., 0., 0., 0.]))
    assert validate_yield_vector(values) == values
    del values['Ar(1s2)']
    with pytest.raises(SpecError, match='products'):
        validate_yield_vector(values)


def test_extra_product_refuses():
    from rmgpy.tools.eedf.ownership import validate_yield_vector
    with pytest.raises(SpecError, match='products'):
        validate_yield_vector(dict(ID6, extra=0.))


def test_unresolved_family_refuses(source_files):
    from rmgpy.tools.eedf.ownership import load_cascade_yields, require_cascade_yields
    snapshot = load_cascade_yields(*source_files[:2])
    with pytest.raises(SpecError, match='unresolved'):
        require_cascade_yields(snapshot, 11)


def test_missing_source_family_refuses(source_files):
    from rmgpy.tools.eedf.ownership import load_cascade_yields, require_cascade_yields
    snapshot = load_cascade_yields(*source_files[:2])
    with pytest.raises(SpecError, match='missing'):
        require_cascade_yields(snapshot, 7)


@pytest.mark.parametrize('kwargs', [{'process': 'different'}, {'threshold_eV': 12.}])
def test_source_collision_identity_refuses(source_files, kwargs):
    from rmgpy.tools.eedf.ownership import load_cascade_yields, require_cascade_yields
    with pytest.raises(SpecError, match='collision identity'):
        require_cascade_yields(load_cascade_yields(*source_files[:2]), 6, **kwargs)


@pytest.mark.parametrize('target', ['yield', 'radiative'])
def test_source_pin_mismatch_refuses(source_files, target):
    from rmgpy.tools.eedf.ownership import load_cascade_yields
    kwargs = {'yield_table_sha256' if target == 'yield' else 'source_sha256': '0' * 64}
    with pytest.raises(FingerprintMismatch, match='source hash'):
        load_cascade_yields(*source_files[:2], **kwargs)


def test_radiative_digest_mismatch_refuses(source_files):
    from rmgpy.tools.eedf.ownership import load_cascade_yields
    yields, radiative, _ = source_files
    radiative.write_text('changed source')
    with pytest.raises(FingerprintMismatch, match='radiative source'):
        load_cascade_yields(yields, radiative)


@pytest.mark.parametrize('mutation', ['duplicate_id', 'duplicate_process', 'unknown_status', 'unresolved_with_numbers', 'complete_with_null', 'missing_transition', 'negative_threshold'])
def test_malformed_source_refuses(source_files, mutation):
    from rmgpy.tools.eedf.ownership import load_cascade_yields
    yields, radiative, document = source_files
    row = document['channels'][0]
    if mutation == 'duplicate_id':
        document['channels'].append(copy.deepcopy(row))
    elif mutation == 'duplicate_process':
        document['channels'][1]['process'] = row['process']
    elif mutation == 'unknown_status':
        row['status'] = 'probably_complete'
    elif mutation == 'unresolved_with_numbers':
        row['status'] = 'unresolved'
    elif mutation == 'complete_with_null':
        row['product_yields'] = None
    elif mutation == 'missing_transition':
        row['missing_transitions'] = ['missing']
    else:
        row['threshold_eV'] = -1.
    yields.write_text(yaml.safe_dump(document))
    with pytest.raises(SpecError):
        load_cascade_yields(yields, radiative)


def test_duplicate_yaml_key_refuses(source_files):
    from rmgpy.tools.eedf.ownership import load_cascade_yields
    yields, radiative, _ = source_files
    yields.write_text(yields.read_text() + 'schema_version: 1\n')
    with pytest.raises(SpecError, match='duplicate'):
        load_cascade_yields(yields, radiative)


def test_source_replacement_uses_verified_snapshot(source_files, monkeypatch):
    from rmgpy.tools.eedf.ownership import load_cascade_yields, require_cascade_yields
    yields, radiative, _ = source_files
    expected = hashlib.sha256(yields.read_bytes()).hexdigest()
    original = Path.read_bytes
    def replace_after_read(path):
        data = original(path)
        if path == yields:
            path.write_text('replaced after read')
        return data
    monkeypatch.setattr(Path, 'read_bytes', replace_after_read)
    snapshot = load_cascade_yields(yields, radiative, yield_table_sha256=expected)
    assert require_cascade_yields(snapshot, 6)['product_yields'] == ID6


@pytest.mark.parametrize('mutation', ['none', 'sum_error', 'sum_preserving'])
def test_database_generator_and_engine_check_same_file(tmp_path, mutation):
    from rmgpy.tools.eedf.ownership import YIELD_SUM_ATOL, load_cascade_yields, require_cascade_yields
    database = os.environ.get('EEDF_DATABASE_ROOT')
    if not database:
        pytest.skip('EEDF_DATABASE_ROOT supplies the sourced generator and yield document')
    root = Path(database)
    source = root / 'input/eedf/argon_radiative_data.yaml'
    original = root / 'input/eedf/argon_cascade_yields.yaml'
    content = original.read_bytes()
    pin = hashlib.sha256(content).hexdigest()
    candidate = tmp_path / 'candidate.yaml'
    candidate.write_bytes(content)
    if mutation != 'none':
        data = yaml.safe_load(content)
        data['channels'][0]['product_yields']['Ar(1s5)'] += 1e-3
        if mutation == 'sum_preserving':
            data['channels'][0]['product_yields']['Ar(1s4)'] -= 1e-3
        candidate.write_text(yaml.safe_dump(data))
    run = subprocess.run([sys.executable, str(root / 'scripts/gen_argon_cascades.py'),
                          '--check', '--source', str(source), '--output', str(candidate)],
                         capture_output=True, text=True)
    print('generator-check', mutation, 'exit', run.returncode, run.stdout, run.stderr)
    if mutation == 'none':
        assert run.returncode == 0
        snapshot = load_cascade_yields(candidate, source, yield_table_sha256=pin)
        assert require_cascade_yields(snapshot, 6)['product_yields'] == ID6
        assert snapshot['channels'] == database_cascade_channels(root)
        for row in snapshot['channels']:
            if row['status'] == 'complete_for_E1':
                assert abs(math.fsum(row['product_yields'].values()) - 1.) <= YIELD_SUM_ATOL
            else:
                with pytest.raises(SpecError, match='unresolved'):
                    require_cascade_yields(snapshot, row['id'])
    else:
        assert run.returncode != 0
        assert 'generated yields differ' in run.stderr
        with pytest.raises(FingerprintMismatch, match='source hash'):
            load_cascade_yields(candidate, source, yield_table_sha256=pin)
        if mutation == 'sum_error':
            with pytest.raises(SpecError, match='sum'):
                load_cascade_yields(candidate, source)


@pytest.mark.parametrize('state,product', [
    ('3P2', 'Ar(1s5)'), ('3P1', 'Ar(1s4)'),
    ('3P0', 'Ar(1s3)'), ('1P1', 'Ar(1s2)'),
])
def test_four_s_cascade_matches_molecular_source_at_binding_and_runtime(branch_fixture, state, product):
    from rmgpy.tools.eedf.ownership import validate_bound_branch
    snapshot, _, _ = branch_fixture
    snapshot['channels'][0]['process'] = PROCESS.replace('Ar(1S0)', 'Ar(' + state + ')')
    reactions = branch_reactions(snapshot)
    source = next(r.products[1] for r in reactions if r.kinetics.branch_product == product)
    source.label = 'source-alias'
    for reaction in reactions:
        reaction.reactants[1] = source
    family = family_from_reactions(snapshot, reactions)
    assert family['process'] == snapshot['channels'][0]['process']
    for reaction in reactions:
        assert validate_bound_branch(family, reaction) == reaction.kinetics.branch_product
    wrong_source = next(r.products[1] for r in reactions if r.kinetics.branch_product != product)
    reactions[0].reactants[1] = wrong_source
    with pytest.raises(SpecError, match='source'):
        family_from_reactions(snapshot, reactions)
    with pytest.raises(SpecError, match='source'):
        validate_bound_branch(family, reactions[0])


@pytest.mark.parametrize('mutation', [None, 'missing_declaration', 'missing_member', 'wrong_member', 'thermal'])
def test_self_return_admission_requires_declared_complete_family(branch_fixture, monkeypatch, mutation):
    from rmgpy.data.kinetics.library import KineticsLibrary
    from rmgpy.rmg.model import CoreEdgeReactionModel
    from rmgpy.kinetics import EEDFChannel, Arrhenius
    import rmgpy.rmg.model as model_module
    _, reactions, _ = branch_fixture
    source = reactions[0].products[1]
    indices = tuple(r.entry.index for r in reactions)
    for r in reactions:
        k = r.kinetics
        r.reactants = [r.reactants[0], source]
        r.kinetics = r.entry.data = EEDFChannel(
            process=k.process.replace('Ar(1S0)', 'Ar(3P2)'), collision_set=k.collision_set,
            branch_product=k.branch_product, branch_yield=k.branch_yield,
            yield_table_sha256=k.yield_table_sha256)
        r.entry.long_desc = f'EEDF branch family = {indices!r}'
        r.duplicate = True
    library = KineticsLibrary(label='PlasmaArgon')
    library.entries = {r.entry.index: r.entry for r in reactions}
    model = CoreEdgeReactionModel()
    monkeypatch.setattr(model_module, 'get_family_library_object', lambda name: library)
    target = reactions[0]
    if mutation == 'missing_declaration':
        target.entry.long_desc = ''
    elif mutation == 'missing_member':
        del library.entries[indices[-1]]
    elif mutation == 'wrong_member':
        target.entry.long_desc = 'EEDF branch family = (1, 2, 3, 4)'
    elif mutation == 'thermal':
        target.kinetics = Arrhenius(A=(1., 'm^3/(mol*s)'))
    found, reaction = model.check_for_existing_reaction(target)
    assert found is (mutation is not None)
    assert reaction is None
