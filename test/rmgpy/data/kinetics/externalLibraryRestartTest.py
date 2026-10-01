import os
import json
import pickle
import shutil

import pytest

from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.library import KineticsLibrary, LibraryReaction
from rmgpy.data.base import Entry, DatabaseError
from rmgpy.electron_balance import get_placement_owner
from rmgpy.kinetics.arrhenius import Arrhenius
from rmgpy.reaction import Reaction


def test_external_library_identity_survives_restart_library_load():
    library_path = os.path.abspath(os.path.join(
        os.path.dirname(__file__), "plasma_local_context_data", "plasma-local-context"
    ))
    database = KineticsDatabase()
    database.load_libraries(os.path.dirname(library_path), libraries=[library_path])

    library_path = os.path.realpath(library_path)
    assert database.external_library_labels[library_path] == "plasma-local-context"
    assert "plasma-local-context" in database.libraries

    database.load_libraries(os.path.join(library_path, "missing-restart-dir"))
    assert database.external_library_labels[library_path] == "plasma-local-context"

    database.load_libraries(os.path.dirname(library_path), libraries=[library_path])
    assert "plasma-local-context" in database.libraries


def test_external_library_paths_are_canonical_and_pickle_provenance():
    library_path = os.path.join(os.path.dirname(__file__), "plasma_local_context_data", "plasma-local-context")
    database = KineticsDatabase()
    database.load_libraries('', libraries=[library_path])
    restored = pickle.loads(pickle.dumps(database))
    canonical = os.path.realpath(library_path)
    assert restored.external_library_labels == {canonical: 'plasma-local-context'}


def test_missing_external_library_path_is_not_silently_skipped(tmp_path):
    database = KineticsDatabase()
    missing = str(tmp_path / 'missing-library')
    database.external_library_labels[missing] = 'missing-library'
    with pytest.raises(IOError, match="Couldn't find kinetics library"):
        database.load_libraries('', libraries=[missing])


def test_seed_provenance_keeps_semantic_owner_separate_from_source_path():
    """The seed reader restores a declared owner, not its external filesystem path."""
    library = KineticsLibrary(label='restart', auto_generated=True)
    source = '/an/absolute/external/Cation_R_Recombination'
    reaction = Reaction(reactants=[], products=[], kinetics=Arrhenius(A=(1, 's^-1')))
    library.entries[1] = Entry(
        index=1, label='fixture', item=reaction, data=reaction.kinetics,
        long_desc=('Originally from reaction library: Cation_R_Recombination\n'
                   '[RMG external library provenance v1] '
                   + json.dumps({'label': 'Cation_R_Recombination', 'source': source})
                   + ' [/RMG external library provenance v1]\n'),
    )
    restored = library.get_library_reactions()[0]
    assert restored.library == 'Cation_R_Recombination'
    assert restored.library_source == source
    assert get_placement_owner(restored) == 'Cation_R_Recombination'


def test_load_all_refuses_external_label_collision(tmp_path):
    source = os.path.join(os.path.dirname(__file__), "plasma_local_context_data", "plasma-local-context")
    duplicate = tmp_path / 'restart' / 'plasma-local-context'
    shutil.copytree(source, duplicate)
    database = KineticsDatabase()
    database.load_libraries('', libraries=[source])
    with pytest.raises(DatabaseError, match='already in use'):
        database.load_libraries(str(duplicate.parent))
    assert len(database.libraries) == 1


@pytest.mark.parametrize('alias_kind', ['relative', 'symlink'])
def test_external_library_alias_admitted_to_edge(tmp_path, monkeypatch, alias_kind):
    from rmgpy.data.rmg import RMGDatabase
    from rmgpy.rmg.model import CoreEdgeReactionModel
    import rmgpy.data.rmg as data_rmg

    library_path = tmp_path / 'external' / 'alias-library'
    library_path.mkdir(parents=True)
    (library_path / 'reactions.py').write_text('name = "alias-library"\n')
    (library_path / 'dictionary.txt').write_text('')
    monkeypatch.chdir(tmp_path)
    if alias_kind == 'relative':
        spelling = 'external/alias-library'
    else:
        alias = tmp_path / 'different-basename'
        alias.symlink_to(library_path, target_is_directory=True)
        spelling = str(alias)
    database = RMGDatabase()
    database.kinetics = KineticsDatabase()
    database.kinetics.library_order = [(spelling, 'Reaction Library')]
    database.kinetics.load_libraries('', libraries=[spelling])
    monkeypatch.setattr(data_rmg, 'database', database)
    model = CoreEdgeReactionModel()
    model.add_reaction_library_to_edge(spelling)
    model.add_seed_mechanism_to_core(spelling)
    assert database.kinetics.library_order == [('alias-library', 'Reaction Library')]
    assert database.kinetics.external_library_labels == {str(library_path): 'alias-library'}


def _recorded_restart_library(label, source):
    library = KineticsLibrary(label='restart_edge', auto_generated=True)
    reaction = Reaction(reactants=[], products=[], kinetics=Arrhenius(A=(1, 's^-1')))
    # The legacy phrase is deliberately present as a comment too, so baseline
    # reads the same source and the regression isolates admission behavior.
    record = '[RMG external library provenance v1] ' + json.dumps({'label': label, 'source': source}) + ' [/RMG external library provenance v1]'
    library.entries[1] = Entry(index=1, label='fixture', item=reaction, data=reaction.kinetics,
                              long_desc='Originally from reaction library: ' + label + '\n'
                              + 'External library source: ' + source + '\n' + record + '\n')
    return library


@pytest.mark.parametrize('already_loaded', [False, True])
def test_restart_edge_refuses_missing_recorded_source(tmp_path, monkeypatch, already_loaded):
    from rmgpy.data.rmg import RMGDatabase
    from rmgpy.rmg.model import CoreEdgeReactionModel
    import rmgpy.data.rmg as data_rmg

    label = 'plasma-local-context'
    database = RMGDatabase()
    database.kinetics = KineticsDatabase()
    if already_loaded:
        source = os.path.join(os.path.dirname(__file__), 'plasma_local_context_data', label)
        database.kinetics.load_libraries('', libraries=[source])
    missing = str(tmp_path / 'missing' / label)
    database.kinetics.libraries['restart_edge'] = _recorded_restart_library(label, missing)
    monkeypatch.setattr(data_rmg, 'database', database)
    with pytest.raises(IOError, match='External kinetics library .* recorded by restart seed is missing'):
        CoreEdgeReactionModel().add_reaction_library_to_edge('restart_edge')


def test_restart_edge_refuses_conflicting_loaded_source(tmp_path, monkeypatch):
    from rmgpy.data.rmg import RMGDatabase
    from rmgpy.rmg.model import CoreEdgeReactionModel
    import rmgpy.data.rmg as data_rmg

    label = 'plasma-local-context'
    source = os.path.join(os.path.dirname(__file__), 'plasma_local_context_data', label)
    recorded_source = tmp_path / 'different' / label
    shutil.copytree(source, recorded_source)
    database = RMGDatabase()
    database.kinetics = KineticsDatabase()
    database.kinetics.load_libraries('', libraries=[source])
    database.kinetics.libraries['restart_edge'] = _recorded_restart_library(label, str(recorded_source))
    monkeypatch.setattr(data_rmg, 'database', database)
    with pytest.raises(DatabaseError, match='already in use'):
        CoreEdgeReactionModel().add_reaction_library_to_edge('restart_edge')


def test_plain_source_comment_is_not_restart_provenance():
    library = _recorded_restart_library('owner', '/plain/comment/owner')
    library.entries[1].long_desc = ('Originally from reaction library: owner\n'
                                    'External library source: /plain/comment/owner\n')
    assert not hasattr(library.get_library_reactions()[0], 'library_source')


@pytest.mark.parametrize('conflict', ['source', 'owner', 'duplicate-json-key', 'malformed'])
def test_conflicting_or_invalid_provenance_is_refused(conflict):
    library = _recorded_restart_library('owner', '/external/owner')
    entry = library.entries[1]
    if conflict == 'source':
        entry.long_desc += ('[RMG external library provenance v1] '
                            '{"label": "owner", "source": "/different/owner"} '
                            '[/RMG external library provenance v1]\n')
    elif conflict == 'owner':
        entry.long_desc += 'Originally from reaction library: other-owner\n'
    elif conflict == 'duplicate-json-key':
        entry.long_desc = ('Originally from reaction library: owner\n'
                          '[RMG external library provenance v1] '
                          '{"label": "owner", "source": "/external/owner", "source": "/different/owner"} '
                          '[/RMG external library provenance v1]\n')
    else:
        entry.long_desc += '[RMG external library provenance v1] malformed\n'
    with pytest.raises(DatabaseError, match='provenance'):
        library.get_library_reactions()


def test_reserved_provenance_survives_library_writer_and_reader(tmp_path):
    source = os.path.join(os.path.dirname(__file__), 'plasma_local_context_data', 'plasma-local-context')
    database = KineticsDatabase()
    database.load_libraries('', libraries=[source])
    library = database.libraries['plasma-local-context']
    recorded_source = str(tmp_path / "quoted'\\tabbed" / 'plasma-local-context')
    seed = KineticsLibrary(label='restart', auto_generated=True)
    seed.entries = library.entries
    for entry in seed.entries.values():
        entry.long_desc = _recorded_restart_library('plasma-local-context', recorded_source).entries[1].long_desc
    seed_path = tmp_path / 'seed'
    seed_path.mkdir()
    seed.save(str(seed_path / 'reactions.py'))
    restored = KineticsLibrary(label='restart')
    restored.load(str(seed_path / 'reactions.py'), database.local_context, database.global_context)
    reactions = restored.get_library_reactions()
    assert reactions
    assert all(rxn.library_source == recorded_source and rxn.library == 'plasma-local-context' for rxn in reactions)


def test_external_library_alias_selected_for_output(tmp_path, monkeypatch):
    from rmgpy.data.rmg import RMGDatabase
    from rmgpy.rmg.model import CoreEdgeReactionModel
    import rmgpy.data.rmg as data_rmg

    library_path = tmp_path / 'external' / 'alias-library'
    library_path.mkdir(parents=True)
    (library_path / 'reactions.py').write_text('name = "alias-library"\n')
    (library_path / 'dictionary.txt').write_text('')
    monkeypatch.chdir(tmp_path)
    spelling = 'external/alias-library'
    database = RMGDatabase()
    database.kinetics = KineticsDatabase()
    database.kinetics.load_libraries('', libraries=[spelling])
    monkeypatch.setattr(data_rmg, 'database', database)
    model = CoreEdgeReactionModel()
    reaction = LibraryReaction(library='alias-library', reactants=[], products=[])
    model.edge.reactions.append(reaction)
    model.add_reaction_library_to_output(spelling)
    assert model.output_reaction_list == [reaction]


def test_restart_edge_refuses_source_with_different_label(monkeypatch):
    from rmgpy.data.rmg import RMGDatabase
    from rmgpy.rmg.model import CoreEdgeReactionModel
    import rmgpy.data.rmg as data_rmg

    source = os.path.abspath(os.path.join(os.path.dirname(__file__), 'plasma_local_context_data', 'plasma-local-context'))
    database = RMGDatabase()
    database.kinetics = KineticsDatabase()
    database.kinetics.libraries['restart_edge'] = _recorded_restart_library('different-label', source)
    monkeypatch.setattr(data_rmg, 'database', database)
    with pytest.raises(DatabaseError, match='does not provide the restart label'):
        CoreEdgeReactionModel().add_reaction_library_to_edge('restart_edge')


def test_cached_library_records_explicit_external_source():
    source = os.path.abspath(os.path.join(os.path.dirname(__file__), 'plasma_local_context_data', 'plasma-local-context'))
    database = KineticsDatabase()
    database.load_libraries(os.path.dirname(source), libraries=['plasma-local-context'])
    original = database.libraries['plasma-local-context']
    database.load_libraries('', libraries=[source])
    assert database.libraries['plasma-local-context'] is original
    assert database.external_library_labels == {os.path.realpath(source): 'plasma-local-context'}


def _database_with_semantic_label_path_conflict(tmp_path, monkeypatch):
    label = 'plasma-local-context'
    source = os.path.join(os.path.dirname(__file__), 'plasma_local_context_data', label)
    other_label = 'surface-charge-transfer-bep'
    other = os.path.join(os.path.dirname(__file__), 'surface_charge_transfer_bep_data', other_label)
    database = KineticsDatabase()
    database.load_libraries('', libraries=[source, other])
    database.library_order = [(label, 'Reaction Library'), (other_label, 'Reaction Library')]
    (tmp_path / label).symlink_to(other, target_is_directory=True)
    monkeypatch.chdir(tmp_path)
    return database


@pytest.mark.parametrize('operation', ['resolve', 'order'])
def test_loaded_semantic_label_wins_over_cwd_alias(tmp_path, monkeypatch, operation):
    database = _database_with_semantic_label_path_conflict(tmp_path, monkeypatch)
    label = 'plasma-local-context'
    library = database.libraries[label]
    if operation == 'resolve':
        assert database.resolve_library(label) is library
    else:
        database.load_libraries('', libraries=[])
    assert database.library_order == [(label, 'Reaction Library'),
                                      ('surface-charge-transfer-bep', 'Reaction Library')]


@pytest.mark.parametrize('additive', [False, True])
def test_load_all_rescan_order_and_additive_restart(tmp_path, additive):
    label = 'plasma-local-context'
    source = os.path.join(os.path.dirname(__file__), 'plasma_local_context_data', label)
    database = KineticsDatabase()
    database.load_libraries('', libraries=[source])
    library = database.libraries[label]
    previous_order = [(label, 'Reaction Library')]
    database.library_order = list(previous_order)
    keywords = {'additive': True} if additive else {}
    database.load_libraries(str(tmp_path), **keywords)
    assert database.library_order == (previous_order if additive else [])

    restart = tmp_path / 'restart'
    restart.mkdir()
    (restart / 'reactions.py').write_text('name = "restart"\n')
    (restart / 'dictionary.txt').write_text('')
    database.load_libraries(str(tmp_path), **keywords)
    assert database.library_order == ((previous_order if additive else [])
                                      + [('restart', 'Reaction Library')])
    assert database.libraries[label] is library
    assert database.external_library_labels[os.path.realpath(source)] == label


@pytest.mark.parametrize('request_kind, repeat', [
    pytest.param(kind, repeat, id=kind + '-' + repeat)
    for kind in ('database_name', 'absolute_external', 'cwd_relative_external',
                 'same_source_symlink', 'trailing_slash')
    for repeat in ('once', 'twice_same_call', 'twice_across_calls')
] + [
    pytest.param('different_root', 'refused', id='different-root-refused'),
    pytest.param('label_path_ambiguity', 'refused', id='label-path-ambiguity-refused'),
])
def test_library_request_matrix(tmp_path, monkeypatch, request_kind, repeat):
    label = 'plasma-local-context'
    fixture = os.path.join(os.path.dirname(__file__), 'plasma_local_context_data', label)
    source = tmp_path / 'external' / label
    shutil.copytree(fixture, source)
    source = source.resolve()
    database_root = tmp_path / 'database-libraries'
    database_root.mkdir()
    cwd = tmp_path / 'cwd'
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    path = str(database_root)
    if request_kind in ('database_name', 'different_root', 'label_path_ambiguity'):
        (database_root / label).symlink_to(source, target_is_directory=True)
        spelling = label
    elif request_kind == 'cwd_relative_external':
        monkeypatch.chdir(source.parent)
        spelling = label
    elif request_kind == 'same_source_symlink':
        (cwd / 'library-alias').symlink_to(source, target_is_directory=True)
        spelling = 'library-alias'
    elif request_kind == 'trailing_slash':
        spelling = str(source) + os.sep
        path += os.sep
    else:
        spelling = str(source)

    calls = []
    loaded_objects = []
    original_load = KineticsLibrary.load

    def counted_load(library, filename, *args, **kwargs):
        calls.append(os.path.realpath(filename))
        loaded_objects.append(library)
        return original_load(library, filename, *args, **kwargs)

    monkeypatch.setattr(KineticsLibrary, 'load', counted_load)
    database = KineticsDatabase()
    database.library_order = [(spelling, 'Reaction Library')]
    if repeat == 'twice_same_call':
        database.load_libraries(path, libraries=[spelling, spelling])
    else:
        database.load_libraries(path, libraries=[spelling])
    library = database.libraries[label]
    if repeat == 'twice_across_calls':
        database.load_libraries(path, libraries=[spelling])
        assert database.libraries[label] is library
    elif repeat == 'refused':
        other = os.path.join(os.path.dirname(__file__), 'surface_charge_transfer_bep_data',
                             'surface-charge-transfer-bep')
        if request_kind == 'different_root':
            second_root = tmp_path / 'different-root'
            second_root.mkdir()
            (second_root / label).symlink_to(other, target_is_directory=True)
            path = str(second_root)
        else:
            (cwd / label).symlink_to(other, target_is_directory=True)
        with pytest.raises(DatabaseError):
            database.load_libraries(path, libraries=[label])
        assert database.libraries[label] is library

    assert library.label == label
    assert database.libraries == {label: library}
    assert database.library_sources == {label: str(source)}
    assert database.library_order == [(label, 'Reaction Library')]
    assert calls == [str(source / 'reactions.py')]
    assert library is loaded_objects[0]
    if request_kind not in ('database_name', 'different_root', 'label_path_ambiguity'):
        assert database.external_library_labels == {str(source): label}
