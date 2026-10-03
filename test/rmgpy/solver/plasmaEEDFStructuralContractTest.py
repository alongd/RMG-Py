#!/usr/bin/env python3
###############################################################################
# RMG - Reaction Mechanism Generator                                          #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),              #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)       #
# Permission is hereby granted, free of charge, to any person obtaining a copy #
# of this software and associated documentation files (the "Software"), to     #
# deal in the Software without restriction, including without limitation the   #
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or   #
# sell copies of the Software, and to permit persons to whom the Software is   #
# furnished to do so, subject to the following conditions:                     #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER          #
# DEALINGS IN THE SOFTWARE.                                                  #
###############################################################################
"""Old/tip witnesses for the structural Report-10 database contract."""

import json
import random
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest
import rmgpy
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.species import Species

import plasmaEEDFRedTest as red
from plasmaEEDFStateIsolationTest import _database_view, _configured_view


def _private_library(view, name):
    library = view / 'kinetics/libraries' / name
    if library.is_symlink():
        source = library.resolve()
        library.unlink()
        shutil.copytree(source, library)
    return library


def _failure(call, pattern):
    with pytest.raises(pytest.fail.Exception, match=pattern) as caught:
        call()
    print(str(caught.value))
    return str(caught.value)


@pytest.mark.database
def test_report10_contract_pins_recombination_entry_zero(tmp_path, monkeypatch):
    """An inactive pinned entry cannot silently become an active argon channel."""
    view = _database_view(tmp_path)
    library = _private_library(view, 'PlasmaRadiativeRecombination')
    file = library / 'reactions.py'
    original = file.read_text()
    start, end = original.index('entry(\n    index = 0,'), original.index('entry(\n    index = 1,')
    entry = original[start:end]
    changed = entry.replace('[Lip] => [Li]', '[Ars] + [Ar] => [Ar] + [Ar]', 1)
    changed = changed.replace('BadnellRRArrhenius(Z=3, N=2)',
                              "Arrhenius(A=(3e-15, 'cm^3/(molecule*s)'), n=0, Ea=(0, 'J/mol'), T0=(1, 'K'))", 1)
    assert changed != entry
    file.write_text(original[:start] + changed + original[end:])
    dictionary = library / 'dictionary.txt'
    dictionary.write_text(dictionary.read_text() + '\n[Ars]\nmultiplicity 3\n1 Ar u2 p3 c0\n')
    with _configured_view(view, monkeypatch):
        message = _failure(red._report10_database,
                           r'PlasmaRadiativeRecombination:0:.*serialisation differs:.*reaction.reactants')
        assert 'kinetics.class' in message


@pytest.mark.database
def test_report10_contract_names_engine_dispatch_flags(tmp_path, monkeypatch):
    """Two disabled Te dispatch flags are named together, with engine kf values."""
    view = _database_view(tmp_path)
    file = view / 'kinetics/libraries/PlasmaArgon/reactions.py'
    original = file.read_text()
    changed = original.replace('1.187793909e-13', '1.3065732999e-13', 1)
    assert changed != original
    file.write_text(changed)
    load = KineticsLibrary.load

    def changed_load(library, *args, **kwargs):
        result = load(library, *args, **kwargs)
        if library.label == 'PlasmaArgon':
            for index in (86, 88):
                library.entries[index].data.uses_electron_temperature = False
        return result

    monkeypatch.setattr(KineticsLibrary, 'load', changed_load)
    with _configured_view(view, monkeypatch):
        message = _failure(red._report10_database,
                           r'PlasmaArgon:86: rate differs.*PlasmaArgon:88: rate differs')
        assert message.count('kinetics.uses_electron_temperature') == 2
        assert 'expected=' in message and 'actual=0' in message and 'relative difference=1' in message
        assert 'expected=2.0830700188811362e-21, actual=0' in message
        assert 'actual=1.561498061157321e-88' in message


@pytest.mark.database
@pytest.mark.parametrize('consumer', ['preflight', 'view'])
@pytest.mark.parametrize('problem', ['library_file', 'required_directory', 'unreadable_file',
                                    'unreadable_library', 'missing_file', 'root_file', 'recombination_file',
                                    'absent_library', 'absent_root'])
def test_report10_contract_classifies_input_paths(tmp_path, monkeypatch, consumer, problem):
    """Only genuine absence skips; present corrupt input fails with its path."""
    source = tmp_path / 'configured'
    library = source / 'kinetics/libraries/PlasmaArgon'
    library.mkdir(parents=True)
    reactions = library / 'reactions.py'
    reactions.write_text('')
    (library / 'dictionary.txt').write_text('')
    damaged = library
    if problem in ('library_file', 'absent_library'):
        shutil.rmtree(library)
        if problem == 'library_file':
            library.write_text('not a directory')
    elif problem == 'recombination_file':
        damaged = source / 'kinetics/libraries/PlasmaRadiativeRecombination'
        damaged.write_text('not a directory')
    elif problem == 'required_directory':
        reactions.unlink()
        reactions.mkdir()
        damaged = reactions
    elif problem == 'unreadable_file':
        reactions.chmod(0)
        damaged = reactions
    elif problem == 'unreadable_library':
        library.chmod(0)
    elif problem == 'missing_file':
        reactions.unlink()
        damaged = reactions
    elif problem in ('root_file', 'absent_root'):
        shutil.rmtree(source)
        if problem == 'root_file':
            source.write_text('not a directory')
        damaged = source
    monkeypatch.setattr(rmgpy.settings, 'require_database_directory', lambda: str(source))
    call = red._report10_database if consumer == 'preflight' else lambda: _database_view(tmp_path / 'scratch')
    try:
        if problem.startswith('absent_'):
            with pytest.raises(pytest.skip.Exception) as caught:
                call()
            assert str(damaged) in str(caught.value) or str(source) in str(caught.value)
            print(str(caught.value))
        else:
            message = _failure(call, r'{0}.*(expected|unreadable|lacks required content)'.format(damaged))
            assert str(damaged) in message
    finally:
        if problem == 'unreadable_library':
            library.chmod(0o755)
        elif problem == 'unreadable_file':
            reactions.chmod(0o644)


def _random_draws():
    # Choose an entry, then discover scalar/quantity fields from its entire
    # pinned snapshot. No kinetics-class parameter list drives these draws.
    rng = random.Random(320)
    libraries = json.loads((red.DATA / 'report10_database_contract.json').read_text())['libraries']
    entries = [(name, int(index)) for name, library in sorted(libraries.items())
               for index in sorted(library, key=int)]
    draws = []
    for draw in range(24):
        name, index = rng.choice(entries)
        entry = libraries[name][str(index)]['entry']
        fields = []
        for owner, schema in (('entry', entry), ('reaction', entry['reaction']),
                              ('kinetics', entry['kinetics'])):
            for attribute, value in sorted(schema.items()):
                if attribute == 'class':
                    continue
                if owner == 'kinetics' and attribute == 'T0' and entry['kinetics'].get('n', {}).get('value_si') == 0:
                    continue
                if isinstance(value, (bool, int, float)) or (
                        isinstance(value, dict) and 'value_si' in value):
                    fields.append(owner + '.' + attribute)
                elif owner == 'kinetics' and attribute in ('Tmin', 'Tmax', 'Pmin', 'Pmax'):
                    fields.append(owner + '.' + attribute)
                elif owner == 'reaction' and attribute == 'specific_collider':
                    fields.append(owner + '.' + attribute)
        draws.append((draw, name, index, rng.choice(fields)))
    return draws


@pytest.mark.database
@pytest.mark.parametrize('draw, name, index, field', _random_draws())
def test_report10_contract_random_non_neutral_field(tmp_path, monkeypatch, draw, name, index, field):
    """Each draw must name the changed entry and the randomly selected field."""
    load = KineticsLibrary.load
    mutation_applied = []

    def changed_load(library, *args, **kwargs):
        result = load(library, *args, **kwargs)
        if library.label == name:
            entry = library.entries[index]
            owner, attribute = field.split('.')
            target = entry.data if owner == 'kinetics' else entry.item if owner == 'reaction' else entry
            original = getattr(target, attribute)
            if attribute in ('Tmin', 'Tmax', 'Pmin', 'Pmax'):
                units = 'K' if attribute.startswith('T') else 'Pa'
                value = ((original.value_si * 1.1 if original.value_si else 1.)
                         if original is not None else 321.)
                setattr(target, attribute, (value, units))
            elif attribute == 'specific_collider':
                target.specific_collider = Species(label='witness').from_adjacency_list('1 Ar u0 p4 c0')
            elif hasattr(original, 'value'):
                values = np.asarray(original.value).copy()
                if values.ndim:
                    position = np.flatnonzero(values)[0]
                    values.flat[position] *= 1.1
                else:
                    values = float(values) * 1.1 if values else 0.1
                setattr(target, attribute, (values, original.units))
            elif isinstance(original, bool):
                setattr(target, attribute, not original)
            else:
                setattr(target, attribute, original + 1)
            mutation_applied.append((name, index, field))
        return result

    monkeypatch.setattr(KineticsLibrary, 'load', changed_load)
    message = _failure(red._report10_database, r'{0}:{1}:'.format(name, index))
    assert mutation_applied == [(name, index, field)]
    serial_field = field.removeprefix('entry.')
    assert serial_field in message, 'draw {0}: {1}'.format(draw, message)


@pytest.mark.database
def test_report10_contract_snapshot_is_deterministic_and_complete():
    """Reordering participants and neutral labels/prose preserves serialisation."""
    libraries, problems = red._report10_libraries(Path(rmgpy.settings.require_database_directory()))
    assert not problems
    for name, library in libraries.items():
        for entry in library.entries.values():
            before = red._canonical_entry(entry, name)
            entry.item.reactants.reverse()
            entry.item.products.reverse()
            entry.label = 'neutral entry label'
            entry.item.label = 'neutral reaction label'
            entry.short_desc = 'neutral description'
            entry.data.comment = 'neutral kinetics comment'
            for spc in entry.item.reactants + entry.item.products:
                spc.label = 'neutral species label'
            assert red._canonical_entry(entry, name) == before
            # No kinetics data descriptor may escape the dynamic field list.
            for attribute in dir(entry.data):
                if not attribute.startswith('_') and attribute != 'comment':
                    if not callable(getattr(entry.data, attribute)):
                        assert attribute in before['kinetics']



def _metadata_witness(view):
    file = view / 'kinetics/libraries/PlasmaArgon/reactions.py'
    original = file.read_text()
    changed = original.replace('name = "PlasmaArgon"',
                               'name = "PlasmaRadiativeRecombination"\nautoGenerated = True', 1)
    start = changed.index('entry(\n    index = 91,')
    end = changed.index('entry(\n    index = 92,', start)
    entry = changed[start:end]
    changed_entry = entry.replace('longDesc = u"""', 'longDesc = u"""rate rule [Attacher]\n', 1)
    assert changed != original and changed_entry != entry
    file.write_text(changed[:start] + changed_entry + changed[end:])


@pytest.mark.database
@pytest.mark.parametrize('boundary', ['preflight', 'initialized'])
def test_report10_contract_names_metadata_removed_entry(tmp_path, monkeypatch, boundary):
    """Library metadata and parsed prose must name the entry they remove."""
    view = _database_view(tmp_path)
    _metadata_witness(view)
    with _configured_view(view, monkeypatch):
        if boundary == 'preflight':
            message = _failure(red._report10_database, r'PlasmaArgon:91:')
            assert 'library.name' in message and 'library.auto_generated' in message
            assert 'loader_provenance' in message
        else:
            # Exercise the second boundary independently of the preflight.
            monkeypatch.setattr(red, '_report10_database', lambda: None)
            output = tmp_path / 'output'
            output.mkdir()
            _failure(lambda: red._report10_deck(output, monkeypatch),
                     r'PlasmaArgon:91: missing initialized reaction')


@pytest.mark.database
@pytest.mark.parametrize('consumer', ['preflight', 'view'])
@pytest.mark.parametrize('component', ['', 'kinetics', 'kinetics/libraries',
                                     'kinetics/libraries/PlasmaArgon'])
@pytest.mark.parametrize('damage', ['dangling', 'file', 'unreadable', 'absent'])
def test_report10_contract_classifies_every_component(tmp_path, monkeypatch, consumer, component, damage):
    """Ancestor corruption names the component; genuine absence alone skips."""
    source = tmp_path / 'configured'
    damaged = source / component if component else source
    damaged.parent.mkdir(parents=True, exist_ok=True)
    if damage == 'dangling':
        damaged.symlink_to(tmp_path / 'nonexistent')
    elif damage == 'file':
        damaged.write_text('not a directory')
    elif damage == 'unreadable':
        damaged.mkdir()
        damaged.chmod(0)
    monkeypatch.setattr(rmgpy.settings, 'require_database_directory', lambda: str(source))
    call = red._report10_database if consumer == 'preflight' else lambda: _database_view(tmp_path / 'scratch')
    try:
        if damage == 'absent':
            with pytest.raises(pytest.skip.Exception) as caught:
                call()
            message = str(caught.value)
            print(message)
            assert message == 'Report-10 input is absent: {0}'.format(damaged)
        else:
            message = _failure(call, r'{0}:.*(dangling|expected|unreadable)'.format(damaged))
            assert 'input {0}:'.format(damaged) in message
    finally:
        if damage == 'unreadable':
            damaged.chmod(0o755)



@pytest.mark.database
def test_report10_contract_names_all_missing_and_extra_initialized_reactions(tmp_path, monkeypatch):
    """Missing and extra identities are reported together, including new libraries."""
    import copy

    reactor, _, reactions = red._report10_deck(tmp_path, monkeypatch)
    for reaction in reactions:
        if reaction.library == 'PlasmaArgon' and reaction.entry.index in (86, 91):
            old_index = reaction.entry.index
            reaction.entry = copy.copy(reaction.entry)
            reaction.entry.index = 94 if old_index == 86 else 5
            if old_index == 91:
                reaction.library = 'UnexpectedLibrary'
    message = _failure(lambda: red._report10_initialized_rates(reactor, reactions),
                       r'PlasmaArgon:86: missing initialized reaction')
    for identity, problem in [('PlasmaArgon:86', 'missing'), ('PlasmaArgon:91', 'missing'),
                              ('PlasmaArgon:94', 'extra'), ('UnexpectedLibrary:5', 'extra')]:
        assert identity + ': ' + problem + ' initialized reaction' in message
    assert 'rate evaluation failed' not in message


@pytest.mark.database
@pytest.mark.parametrize('name', ['PlasmaArgon', 'PlasmaRadiativeRecombination'])
@pytest.mark.parametrize('field', ['name', 'auto_generated', 'label', 'solvent', 'metal',
                                 'site', 'facet', 'top', 'thermo_convention', 'new_loader_metadata'])
def test_report10_contract_discovers_library_metadata(tmp_path, monkeypatch, name, field):
    """Existing and newly introduced library fields cannot escape discovery."""
    load = KineticsLibrary.load

    def changed_load(library, *args, **kwargs):
        result = load(library, *args, **kwargs)
        if library.label == name:
            setattr(library, field, True if field == 'auto_generated' else
                    ['changed'] if field == 'top' else 'changed')
        return result

    monkeypatch.setattr(KineticsLibrary, 'load', changed_load)
    message = _failure(red._report10_database, r'{0}:.*library\.{1}'.format(name, field))
    pinned = json.loads((red.DATA / 'report10_database_contract.json').read_text())['libraries'][name]
    for index in pinned:
        assert '{0}:{1}: library serialisation differs: library.{2}'.format(name, index, field) in message



@pytest.mark.database
def test_report10_contract_names_initialized_structure_callback(tmp_path, monkeypatch):
    """A loaded callback changes products while leaving raw entries/rates alone."""
    view = _database_view(tmp_path)
    file = view / 'kinetics/libraries/PlasmaArgon/reactions.py'
    file.write_text(file.read_text() + """
def initialized_reactions(original=entry.__self__.get_library_reactions):
    reactions = original()
    for reaction in reactions:
        if reaction.entry.index == 92:
            reaction.products[:] = [reaction.reactants[0], reaction.reactants[0]]
    return reactions
entry.__self__.get_library_reactions = initialized_reactions
""")
    with _configured_view(view, monkeypatch):
        red._report10_database()  # The raw entry and rate contract is unchanged.
        output = tmp_path / 'output'
        output.mkdir()
        message = _failure(lambda: red._report10_deck(output, monkeypatch),
                           r'PlasmaArgon:92: initialized structure differs: products')
        assert 'rate differs' not in message


@pytest.mark.database
@pytest.mark.parametrize('consumer', ['preflight', 'view'])
@pytest.mark.parametrize('ancestor', ['ancestor', 'ancestor/nested'])
@pytest.mark.parametrize('damage', ['dangling', 'file', 'unreadable', 'absent'])
def test_report10_contract_classifies_above_root(tmp_path, monkeypatch, consumer, ancestor, damage):
    """The first absent or corrupt component above the configured root is named."""
    damaged = tmp_path / 'configured' / ancestor
    damaged.parent.mkdir(parents=True)
    if damage == 'dangling':
        damaged.symlink_to(tmp_path / 'nonexistent')
    elif damage == 'file':
        damaged.write_text('not a directory')
    elif damage == 'unreadable':
        damaged.mkdir()
        damaged.chmod(0)
    configured = damaged / 'input'
    call = red._report10_database if consumer == 'preflight' else lambda: _database_view(tmp_path / 'scratch')
    try:
        with _configured_view(configured, monkeypatch):
            if damage == 'absent':
                with pytest.raises(pytest.skip.Exception) as caught:
                    call()
                print(str(caught.value))
                assert str(caught.value) == 'Report-10 input is absent: {0}'.format(damaged)
            else:
                _failure(call, r'input {0}:.*(dangling|expected|unreadable)'.format(damaged))
    finally:
        if damage == 'unreadable':
            damaged.chmod(0o755)



@pytest.mark.database
def test_report10_contract_names_every_initialized_structure_change(tmp_path, monkeypatch):
    """Two changed solver product rows are both named without changing rates."""
    reactor, core, reactions = red._report10_deck(tmp_path, monkeypatch)
    metastable = next(index for index, spc in enumerate(core) if spc.label == 'Ars')
    for reaction in reactions:
        if reaction.library == 'PlasmaArgon' and reaction.entry.index in (92, 93):
            row = reactor.product_indices[reactor.reaction_index[reaction]]
            row[row >= 0] = metastable
    message = _failure(lambda: red._report10_initialized_rates(reactor, reactions),
                       r'PlasmaArgon:92: initialized structure differs: products')
    assert 'PlasmaArgon:93: initialized structure differs: products' in message
    assert 'rate differs' not in message


@pytest.mark.database
@pytest.mark.parametrize('name', ['PlasmaArgon', 'PlasmaRadiativeRecombination'])
def test_report10_contract_names_thermo_convention_edit(tmp_path, monkeypatch, name):
    """A file-level convention edit names every affected entry without rate errors."""
    red._report10_database()
    view = _database_view(tmp_path)
    library = _private_library(view, name)
    file = library / 'reactions.py'
    original = file.read_text()
    # KineticsLibrary.load does not read the thermoConvention DSL keyword.
    # Edit its actual loaded field through the bound entry owner instead.
    mutation = "\nentry.__self__.thermo_convention = 'ion'\n"
    assert mutation not in original
    file.write_text(original + mutation)
    with _configured_view(view, monkeypatch):
        message = _failure(red._report10_database,
                           r'{0}:.*library\.thermo_convention'.format(name))
        pinned = json.loads((red.DATA / 'report10_database_contract.json').read_text())['libraries'][name]
        for index in pinned:
            assert '{0}:{1}: library serialisation differs: library.thermo_convention'.format(name, index) in message
        assert 'rate evaluation failed' not in message
        assert 'rate differs' not in message


@pytest.mark.database
def test_report10_contract_uses_library_verified_thermo():
    """Every rate probe is value-matched; none relies on a caller assertion."""
    path, _ = red._report10_database()
    libraries, problems = red._report10_libraries(path)
    assert not problems
    with red._report10_thermo(path) as thermo:
        assert thermo.library_order == list(red.REPORT10_THERMO_LIBRARIES)
        for name, library in libraries.items():
            for entry in library.entries.values():
                reactor, _ = red._entry_reactor(entry, name, thermo)
                assert reactor.thermo_source_assertions == {}
                assert reactor.thermo_provenance_diagnostics
                for diagnostic in reactor.thermo_provenance_diagnostics.values():
                    assert diagnostic.startswith('value-matched to library PlasmaThermo/')
                    assert 'caller-asserted' not in diagnostic


@pytest.mark.database
@pytest.mark.parametrize('change_kinetics, thermo_failure', [
    pytest.param(False, 'syntax', id='False'),
    pytest.param(True, 'syntax', id='True'),
    pytest.param(True, 'cause_cycle', id='cause-cycle'),
    pytest.param(True, 'context_cycle', id='context-cycle'),
    pytest.param(True, 'depth_limit', id='depth-limit'),
])
def test_report10_contract_keeps_entry_differences_after_thermo_load_failure(
        tmp_path, monkeypatch, change_kinetics, thermo_failure):
    """Broken thermo must fail alongside every independently readable raw edit."""
    view = _database_view(tmp_path)
    thermo = view / 'thermo'
    source = thermo.resolve()
    thermo.unlink()
    thermo.mkdir()
    for child in source.iterdir():
        if child.name != 'libraries':
            (thermo / child.name).symlink_to(child, target_is_directory=child.is_dir())
    libraries = thermo / 'libraries'
    libraries.mkdir()
    for child in (source / 'libraries').iterdir():
        target = libraries / child.name
        if child.name == 'PlasmaThermo.py':
            shutil.copy2(child, target)
        else:
            target.symlink_to(child, target_is_directory=child.is_dir())
    broken = libraries / 'PlasmaThermo.py'
    cycle = """
Failure = entry.__self__.load_entry.__globals__['DatabaseError']
try:
    raise Failure('original thermo failure')
except Failure as original:
    try:
        raise Failure('thermo recovery failure') from original
    except Failure as recovery:
        raise original from recovery
"""
    context_cycle = """
Failure = entry.__self__.load_entry.__globals__['DatabaseError']
original = Failure('original thermo failure')
recovery = Failure('thermo recovery failure')
original.__context__ = recovery
recovery.__context__ = original
raise original
"""
    depth = """
Failure = entry.__self__.load_entry.__globals__['DatabaseError']
failure = Failure('original thermo failure')
number = 0
while number < 80:
    outer = Failure('nested thermo failure')
    outer.__cause__ = failure
    failure = outer
    number += 1
raise failure
"""
    failures = {'syntax': '\n! invalid syntax\n', 'cause_cycle': cycle,
                'context_cycle': context_cycle, 'depth_limit': depth}
    broken.write_text(broken.read_text() + failures[thermo_failure])
    if change_kinetics:
        file = view / 'kinetics/libraries/PlasmaArgon/reactions.py'
        content = file.read_text()
        for old, new in [('1.187793909e-13', '1.3065732999e-13'), ('2.0e-13', '2.2e-13')]:
            assert old in content
            content = content.replace(old, new, 1)
        file.write_text(content)

    def unexpected_rate_probe(*args):
        raise AssertionError('thermo-dependent rates must not be evaluated')

    monkeypatch.setattr(red, '_entry_reactor', unexpected_rate_probe)
    # Bound the witness itself so a lost termination guard fails promptly.
    class UnterminatedChain(BaseException):
        pass

    walker = red._report10_thermo.__wrapped__.__code__
    steps = 0

    def limit_walk(frame, event, arg):
        nonlocal steps
        if frame.f_code is not walker:
            return None
        if event == 'line':
            steps += 1
            if steps > 1024:
                raise UnterminatedChain('thermo diagnostic did not terminate')
        return limit_walk

    saved_database = red.rmg_data_module.database
    saved_trace = sys.gettrace()
    if thermo_failure != 'syntax':
        sys.settrace(limit_walk)
    try:
        with _configured_view(view, monkeypatch):
            message = _failure(red._report10_database,
                               r'thermo library PlasmaThermo: loading failed:')
    finally:
        sys.settrace(saved_trace)
    assert red.rmg_data_module.database is saved_database
    if thermo_failure == 'syntax':
        assert "AttributeError: 'NoneType' object has no attribute 'tb_lineno'" in message
        assert 'SyntaxError: invalid syntax' in message
    elif thermo_failure.endswith('cycle'):
        assert message.count('DatabaseError: original thermo failure') == 1
        assert message.count('DatabaseError: thermo recovery failure') == 1
        assert message.count('exception chain cycle detected') == 1
    else:
        assert 'exception chain depth limit (64) reached' in message
        assert message.count('DatabaseError: nested thermo failure') == 64
    assert 'rate evaluation failed' not in message
    assert 'rate differs' not in message
    assert message.count('rate not evaluated because thermo failed to load') == 10
    if change_kinetics:
        for index in (88, 91):
            assert 'PlasmaArgon:{0}: serialisation differs: kinetics.A.value_si'.format(index) in message
    else:
        assert 'serialisation differs' not in message
