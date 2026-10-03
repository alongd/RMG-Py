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
"""Regression for the real-deck helper's process-global isolation."""

from contextlib import contextmanager
import logging
import shutil
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
import rmgpy
import rmgpy.data.rmg as rmg_data_module
import rmgpy.rmg.input as rmg_input
from rmgpy import electron_placement  # Load lazy dependencies before the globals baseline.
from rmgpy.molecule import symmetry
from rmgpy.data.kinetics import quarantine
from rmgpy.exceptions import SettingsError
from rmgpy.rmg.main import RMG

from plasmaEEDFRedTest import (
    _report10_database,
    _report10_deck,
    test_report10_omits_resonance_power_and_states_but_closes as _characterise_report10,
)


def _freeze_global(value, seen=None):
    """Detect in-place container changes as well as changed module bindings."""
    seen = set() if seen is None else seen
    if id(value) in seen:
        return ('cycle', id(value))
    seen = seen | {id(value)}
    if isinstance(value, dict):
        return tuple(sorted((repr(key), _freeze_global(item, seen)) for key, item in value.items()))
    if isinstance(value, (set, frozenset)):
        return tuple(sorted(repr(item) for item in value))
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_global(item, seen) for item in value)
    if isinstance(value, np.ndarray):
        return value.shape, value.dtype.str, value.tobytes()
    return repr(value)


def _rmg_globals():
    """Snapshot every loaded rmgpy module, without a cache-name allowlist."""
    return {name + '.' + key: (id(value), _freeze_global(value))
            for name, module in sys.modules.copy().items()
            if name == 'rmgpy' or name.startswith('rmgpy.')
            for key, value in vars(module).copy().items()}


def _globals_diff(before):
    after = _rmg_globals()
    return sorted(key for key in before.keys() | after.keys() if before.get(key) != after.get(key))


def _logger_state():
    return [(logger, logger.handlers, tuple(logger.handlers), logger.level, logger.propagate)
            for logger in (logging.getLogger(), logging.getLogger('rmgpy'))]


def _assert_loggers_unchanged(before):
    for logger, handlers_object, handlers, level, propagate in before:
        assert logger.handlers is handlers_object
        assert tuple(logger.handlers) == handlers
        assert logger.level == level
        assert logger.propagate == propagate


@contextmanager
def _configured_view(view, monkeypatch):
    settings = rmgpy.settings
    values, sources, filename = dict(settings), dict(settings.sources), settings.filename
    sources_object = settings.sources
    try:
        with monkeypatch.context() as caller:
            settings['database.directory'] = str(view)
            yield caller
    finally:
        dict.clear(settings)
        dict.update(settings, values)
        settings.sources = sources_object
        sources_object.clear()
        sources_object.update(sources)
        settings.filename = filename


@pytest.fixture(autouse=True)
def _restore_regression_process_state():
    """A failing regression must not contaminate the next regression."""
    caches = [(name, value, value.copy()) for name, value in vars(quarantine).items()
              if name.endswith(('_WARNED', '_CACHE')) and isinstance(value, (dict, set))]
    loggers = _logger_state()
    warning_filters = list(warnings.filters)
    registries = [(module, '__warningregistry__' in vars(module), vars(module).get('__warningregistry__'))
                  for name, module in sys.modules.copy().items()
                  if name == 'rmgpy' or name.startswith('rmgpy.')]
    registry_contents = [(module, present, original, None if original is None else original.copy())
                         for module, present, original in registries]
    try:
        yield
    finally:
        for module, present, original, contents in registry_contents:
            if not present:
                vars(module).pop('__warningregistry__', None)
            else:
                module.__warningregistry__ = original
                if original is not None:
                    original.clear()
                    original.update(contents)
        for name, original, contents in caches:
            setattr(quarantine, name, original)
            original.clear()
            original.update(contents)
        warnings.filters[:] = warning_filters
        for logger, handlers_object, handlers, level, propagate in loggers:
            for handler in logger.handlers[:]:
                if handler not in handlers:
                    logger.removeHandler(handler)
                    handler.close()
            logger.handlers = handlers_object
            handlers_object[:] = handlers
            logger.setLevel(level)
            logger.propagate = propagate


def _database_view(tmp_path):
    """Copy only the argon library; keep the configured database read-only."""
    try:
        source = rmgpy.settings.require_database_directory()
    except SettingsError as exc:
        pytest.skip(str(exc))
    source = Path(source)
    if not (source / 'kinetics/libraries/PlasmaArgon/reactions.py').is_file():
        pytest.skip('Configured database lacks kinetics/libraries/PlasmaArgon/reactions.py')
    view = tmp_path / 'database'
    view.mkdir()
    for child in source.iterdir():
        if child.name != 'kinetics':
            (view / child.name).symlink_to(child, target_is_directory=child.is_dir())
    kinetics = view / 'kinetics'
    kinetics.mkdir()
    for child in (source / 'kinetics').iterdir():
        if child.name == 'families':
            # Quarantine manifests require real directory components.
            shutil.copytree(child, kinetics / child.name)
        elif child.name != 'libraries':
            (kinetics / child.name).symlink_to(child, target_is_directory=child.is_dir())
    libraries = kinetics / 'libraries'
    libraries.mkdir()
    for child in (source / 'kinetics/libraries').iterdir():
        if child.name == 'PlasmaArgon':
            shutil.copytree(child, libraries / child.name)
        else:
            (libraries / child.name).symlink_to(child, target_is_directory=child.is_dir())
    return view


@pytest.mark.database
@pytest.mark.parametrize('has_handlers', [False, True])
def test_report10_helper_preserves_configuration_refusal_and_dsl_globals(tmp_path, monkeypatch, has_handlers):
    """A real helper call preserves provenance/DSL state and missing-key refusal."""
    settings = rmgpy.settings
    values, sources, filename = dict(settings), dict(settings.sources), settings.filename
    sources_object = settings.sources
    old_database = rmg_data_module.database
    sentinel_rmg = object()
    sentinel_species = {'prior-species': object()}
    sentinel_fragments = {'prior-molecule': {'prior-fragment': 2}}
    with monkeypatch.context() as outer:
        for logger in (logging.getLogger(), logging.getLogger('rmgpy')):
            outer.setattr(logger, 'handlers', [logging.NullHandler()] if has_handlers else [])
            outer.setattr(logger, 'level', logging.WARNING)
            outer.setattr(logger, 'propagate', has_handlers)
        for name, value in vars(quarantine).copy().items():
            if name.endswith(('_WARNED', '_CACHE')) and isinstance(value, (dict, set)):
                replacement = type(value)()
                if has_handlers:
                    replacement.update({'caller-sentinel': None} if isinstance(value, dict)
                                       else {'caller-sentinel'})
                outer.setattr(quarantine, name, replacement)
        outer.setattr(rmg_input, 'rmg', sentinel_rmg)
        outer.setattr(rmg_input, 'species_dict', sentinel_species)
        # read_input_file creates this global even if no earlier input was read.
        outer.setattr(rmg_input, 'mol_to_frag', sentinel_fragments, raising=False)
        output = tmp_path / 'configured'
        output.mkdir()
        try:
            globals_before, loggers_before = _rmg_globals(), _logger_state()
            filters_before = list(warnings.filters)
            with outer.context() as helper_patch:
                _report10_deck(output, helper_patch)
                # Assert before caller teardown can hide an unrestored binding.
                assert dict(settings) == values and settings.filename == filename
                assert dict(settings.sources) == sources
                assert settings.sources is sources_object
                assert rmg_input.mol_to_frag is sentinel_fragments
                assert rmg_input.rmg is sentinel_rmg
                assert rmg_input.species_dict is sentinel_species
                assert rmg_data_module.database is old_database
                changes = _globals_diff(globals_before)
                print('rmgpy globals diff:', changes)
                _assert_loggers_unchanged(loggers_before)
                assert warnings.filters == filters_before
                assert changes == []

            # An existing default directory makes provenance, rather than a
            # nonexistent-path exception, the load-bearing missing-key guard.
            no_database = tmp_path / 'no-database-rmgrc'
            no_database.write_text('# deliberately no database declaration\n')
            settings.load(str(no_database))
            dict.__setitem__(settings, 'database.directory', values['database.directory'])
            assert settings.sources['database.directory'] == settings.DEFAULT_SOURCE
            with pytest.raises(SettingsError, match="no 'database.directory' line"):
                settings.require_database_directory()
            missing_sources = dict(settings.sources)
            missing_output = tmp_path / 'missing-declaration'
            missing_output.mkdir()
            missing_before = _rmg_globals()
            with outer.context() as helper_patch:
                try:
                    _report10_deck(missing_output, helper_patch)
                except pytest.skip.Exception as exc:
                    # The fixed helper refuses to guess a database, so this
                    # negative helper call is intentionally unavailable.
                    assert 'database' in str(exc).lower()
                with pytest.raises(SettingsError, match="no 'database.directory' line"):
                    settings.require_database_directory()
                assert settings.sources == missing_sources
                assert _globals_diff(missing_before) == []
                _assert_loggers_unchanged(loggers_before)
            with pytest.raises(SettingsError, match="no 'database.directory' line"):
                settings.require_database_directory()
            assert settings.sources == missing_sources
            assert settings.sources is sources_object
            assert rmg_input.mol_to_frag is sentinel_fragments
            assert sentinel_fragments == {'prior-molecule': {'prior-fragment': 2}}
        finally:
            # Bypass Settings.__setitem__: teardown must not change provenance.
            dict.clear(settings)
            dict.update(settings, values)
            settings.sources = sources_object
            sources_object.clear()
            sources_object.update(sources)
            settings.filename = filename
            rmg_data_module.database = old_database


@pytest.mark.database
def test_report10_accepts_isomorphic_library_labels(tmp_path, monkeypatch):
    """Equivalent library labels must still execute the scientific assertions."""
    _report10_database()  # Skip an incompatible original before testing a relabelling.
    view = _database_view(tmp_path)
    library = view / 'kinetics/libraries/PlasmaArgon'
    for name in ('dictionary.txt', 'reactions.py'):
        file = library / name
        file.write_text(file.read_text().replace('Arp', 'ArPlus'))
    with _configured_view(view, monkeypatch) as caller:
        sources, filename = dict(rmgpy.settings.sources), rmgpy.settings.filename
        before, loggers = _rmg_globals(), _logger_state()
        output = tmp_path / 'renamed-output'
        output.mkdir()
        try:
            _characterise_report10(output, caller)
        except pytest.skip.Exception as exc:
            pytest.fail('Isomorphic relabelled library was skipped: {0}'.format(exc))
        assert _globals_diff(before) == []
        _assert_loggers_unchanged(loggers)
        assert rmgpy.settings.sources == sources
        assert rmgpy.settings.filename == filename


@pytest.mark.database
@pytest.mark.parametrize('content', ['molecule', 'rate'])
def test_report10_names_incompatible_content(tmp_path, monkeypatch, content):
    """A same-label molecular or rate change must name the precise content."""
    view = _database_view(tmp_path)
    library = view / 'kinetics/libraries/PlasmaArgon'
    if content == 'molecule':
        file = library / 'dictionary.txt'
        changed = file.read_text().replace('1 Ar u0 p4 c0', 'GROUND_STRUCTURE_PLACEHOLDER')
        changed = changed.replace('multiplicity 3\n1 Ar u2 p3 c0', '1 Ar u0 p4 c0')
        changed = changed.replace('GROUND_STRUCTURE_PLACEHOLDER', 'multiplicity 3\n1 Ar u2 p3 c0')
        file.write_text(changed)
        reason = 'PlasmaArgon:86: reactant molecules'
    else:
        file = library / 'reactions.py'
        file.write_text(file.read_text().replace('1.187793909e-13', '1.3065732999e-13'))
        reason = r'PlasmaArgon:88: rate differs.*characterisation must be regenerated'
    with _configured_view(view, monkeypatch) as caller:
        def should_not_initialize(*args, **kwargs):
            pytest.fail('Incompatible preflight allowed RMG.initialize')

        caller.setattr(RMG, 'initialize', should_not_initialize)
        before, loggers = _rmg_globals(), _logger_state()
        output = tmp_path / 'incompatible-output'
        output.mkdir()
        with pytest.raises(pytest.fail.Exception, match=reason):
            _report10_deck(output, caller)
        assert _globals_diff(before) == []
        _assert_loggers_unchanged(loggers)


@pytest.mark.database
def test_report10_accepts_rate_neutral_entry91_t0(tmp_path, monkeypatch):
    """Entry 91's T0 is rate-neutral because its temperature exponent is zero."""
    view = _database_view(tmp_path)
    file = view / 'kinetics/libraries/PlasmaArgon/reactions.py'
    file.write_text(file.read_text().replace("T0 = (11604.51812, 'K'),\n    ),\n    shortDesc = u\"[Ashida1995 via Rehman2016 Table 1] metastable-to-resonance mixing",
                                             "T0 = (1, 'K'),\n    ),\n    shortDesc = u\"[Ashida1995 via Rehman2016 Table 1] metastable-to-resonance mixing"))
    with _configured_view(view, monkeypatch) as caller:
        output = tmp_path / 'neutral-output'
        output.mkdir()
        _characterise_report10(output, caller)


@pytest.mark.database
@pytest.mark.parametrize('kind, replacements, message', [
    ('cross_section', [('1.16e-22', '1.160001e-22')], r'PlasmaArgon:86: cross section sigma differs'),
    ('nan', [('1.16e-22', "float('nan')")], r'PlasmaArgon:86: nonfinite parameters sigma'),
    ('two_rates', [('1.187793909e-13', '1.3065732999e-13'), ('2.0e-13', '2.2e-13')],
     r'PlasmaArgon:88:.*PlasmaArgon:91:'),
    ('inventory', [('index = 93,', 'index = 94,')], r'PlasmaArgon:94: unexpected reaction.*PlasmaArgon:93: missing entry'),
])
def test_report10_contract_names_every_changed_entry(tmp_path, monkeypatch, kind, replacements, message):
    """The contract names all relevant entries instead of silently skipping a deck."""
    view = _database_view(tmp_path)
    file = view / 'kinetics/libraries/PlasmaArgon/reactions.py'
    text = file.read_text()
    for old, new in replacements:
        text = text.replace(old, new, 1)
    file.write_text(text)
    with _configured_view(view, monkeypatch):
        with pytest.raises(pytest.fail.Exception, match=message):
            _report10_database()


@pytest.mark.database
def test_report10_contract_names_multiarrhenius_representation_change(tmp_path, monkeypatch):
    """A one-law wrapper remains a named representation change, not equivalence."""
    view = _database_view(tmp_path)
    file = view / 'kinetics/libraries/PlasmaArgon/reactions.py'
    text = file.read_text()
    start = text.index('    index = 92,')
    end = text.index('\nentry(', start)
    entry = text[start:end]
    entry = entry.replace('kinetics = Arrhenius(', 'kinetics = MultiArrhenius(arrhenius=[Arrhenius(', 1)
    entry = entry.replace("    ),\n    shortDesc", "        )]),\n    shortDesc", 1)
    file.write_text(text[:start] + entry + text[end:])
    with _configured_view(view, monkeypatch):
        with pytest.raises(pytest.fail.Exception, match=r'PlasmaArgon:92: kinetics type changed'):
            _report10_database()


@pytest.mark.database
def test_report10_contract_rejects_superseded_entry88(tmp_path, monkeypatch):
    """A copied current library with the superseded law fails by named entry."""
    view = _database_view(tmp_path)
    file = view / 'kinetics/libraries/PlasmaArgon/reactions.py'
    text = file.read_text()
    start = text.index('    index = 88,')
    end = text.index('\nentry(', start)
    entry = text[start:end]
    old_law = ("        A = (1.187793909e-13, 'm^3/(molecule*s)'),\n"
               "        n = 0.3278287422,\n"
               "        Ea_g = (4.401535934, 'eV/molecule'),\n"
               "        Ea_e = (4.401535934, 'eV/molecule'),")
    superseded_law = ("        A = (6.8e-15, 'm^3/(molecule*s)'),\n"
                      "        n = 0.67,\n"
                      "        Ea_g = (4.20, 'eV/molecule'),\n"
                      "        Ea_e = (4.20, 'eV/molecule'),")
    assert old_law in entry
    file.write_text(text[:start] + entry.replace(old_law, superseded_law, 1) + text[end:])
    with _configured_view(view, monkeypatch):
        with pytest.raises(pytest.fail.Exception, match=r'PlasmaArgon:88: rate differs'):
            _report10_database()


@pytest.mark.database
@pytest.mark.parametrize('stage', ['preflight', 'initialize'])
def test_report10_restores_state_on_error(tmp_path, monkeypatch, stage):
    """The restoration boundary includes preflight and initialization failures."""
    import plasmaEEDFRedTest as red

    def fail(*args, **kwargs):
        quarantine._LEGACY_CALL_SITES_WARNED.add('injected-warning')
        quarantine._DISK_QUARANTINE_CACHE['injected-cache'] = object()
        for logger in (logging.getLogger(), logging.getLogger('rmgpy')):
            logger.addHandler(logging.NullHandler())
            logger.setLevel(logging.DEBUG)
        raise RuntimeError('injected-' + stage)

    with monkeypatch.context() as caller:
        if stage == 'preflight':
            caller.setattr(red, '_report10_database', fail)
        else:
            caller.setattr(RMG, 'initialize', fail)
        before, loggers = _rmg_globals(), _logger_state()
        with pytest.raises(RuntimeError, match='injected-' + stage):
            _report10_deck(tmp_path, caller)
        assert _globals_diff(before) == []
        _assert_loggers_unchanged(loggers)
