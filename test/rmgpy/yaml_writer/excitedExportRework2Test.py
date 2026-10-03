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

import copy
import multiprocessing
import os
from pathlib import Path
import threading

import pytest
import yaml

from rmgpy.chemkin import load_chemkin_file, save_chemkin, save_species_dictionary
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.exceptions import SpeciesIdentityError
from rmgpy.kinetics import Arrhenius, Lindemann, SurfaceArrhenius, ThirdBody, Troe
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import CoreEdgeReactionModel, ReactionModel
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy import yaml_cantera1, yaml_cantera2, yaml_rms


@pytest.fixture
def state_species():
    species = []
    for index, (smiles, label, electronic, vibrational) in enumerate([
        ('N#N', 'N2', '', -1),
        ('N#N', 'N2(v=1)', '', 1),
        ('N#N', 'N2(A3Su+)', 'A3Su+', -1),
        ('[Ar]', 'Ar(1s5)', '1s5', -1),
    ], 1):
        mol = Molecule(smiles=smiles)
        mol.electronic_state = electronic
        mol.vibrational_level = vibrational
        thermo = NASA(polynomials=[
            NASAPolynomial(coeffs=[3.5, 0, 0, 0, 0, index * 1000, 2], Tmin=(200, 'K'), Tmax=(1000, 'K')),
            NASAPolynomial(coeffs=[3.5, 0, 0, 0, 0, index * 1000, 2], Tmin=(1000, 'K'), Tmax=(3000, 'K')),
        ], Tmin=(200, 'K'), Tmax=(3000, 'K'))
        species.append(Species(index=index, label=label, molecule=[mol], thermo=thermo))
    return species


def falloff(kind=Lindemann, efficiencies=None):
    low = Arrhenius(A=(1e6, 'm^3/(mol*s)'), n=0, Ea=(0, 'J/mol'))
    if kind is ThirdBody:
        return kind(arrheniusLow=low, efficiencies=efficiencies or {})
    high = Arrhenius(A=(1e3, 's^-1'), n=0, Ea=(0, 'J/mol'))
    kwargs = {'alpha': 0.5, 'T3': (100, 'K'), 'T1': (1000, 'K')} if kind is Troe else {}
    return kind(arrheniusHigh=high, arrheniusLow=low, efficiencies=efficiencies or {}, **kwargs)


def library_with(ground, product, collider=None, kinetics=None):
    library = KineticsLibrary(label='state_export')
    reaction = Reaction(reactants=[ground], products=[product], specific_collider=collider)
    collider_label = ' (+{0})'.format(collider.label) if collider else ''
    label = ground.label + collider_label + ' <=> ' + product.label + collider_label
    library.entries[1] = Entry(index=1, label=label, item=reaction,
                               data=kinetics or Arrhenius(A=(1e3, 's^-1'), Ea=(0, 'J/mol')))
    return library


def state_library(ground, product, rate):
    library = library_with(ground, product)
    library.entries[1].data = Arrhenius(A=(rate, 's^-1'), Ea=(0, 'J/mol'))
    return library


def load_library(path):
    return KineticsLibrary().load(str(path), local_context={'Arrhenius': Arrhenius, 'Lindemann': Lindemann})


def loaded_reaction(path):
    library = KineticsLibrary()
    library.load(str(path), local_context={'Arrhenius': Arrhenius, 'Lindemann': Lindemann})
    return next(iter(library.entries.values()))


def assert_library_pair_safe(path, expected_rate, expected_state):
    try:
        entry = loaded_reaction(path)
    except SpeciesIdentityError as error:
        assert str(path) in str(error) and str(path.with_name('dictionary.txt')) in str(error)
    else:
        assert entry.data.A.value_si == expected_rate
        assert entry.item.products[0].molecule[0].vibrational_level == expected_state


def chemkin_model(ground, product, rate):
    model = CoreEdgeReactionModel()
    model.core.species = [ground, product]
    model.core.reactions = [Reaction(reactants=[ground], products=[product],
        kinetics=Arrhenius(A=(rate, 's^-1'), Ea=(0, 'J/mol')))]
    return model


class TestExcitedExportRework2Review:
    def test_review_1_library_crash_between_renames(self, state_species, tmp_path):
        ground, product = state_species[:2]
        ground.label, product.label = 'A', 'B'
        path = tmp_path / 'reactions.py'
        old = copy.deepcopy(product)
        old.molecule[0].vibrational_level = -1
        state_library(ground, old, 1).save(str(path))
        new = state_library(ground, product, 2)
        def crash_save():
            original = os.replace
            def crash(src, dst):
                if Path(dst) == path.with_name('dictionary.txt') and str(src).endswith('.tmp'):
                    os._exit(86)
                return original(src, dst)
            os.replace = crash
            new.save(str(path))
        child = multiprocessing.get_context('fork').Process(target=crash_save)
        child.start()
        child.join(20)
        assert not child.is_alive()
        assert child.exitcode == 86
        assert_library_pair_safe(path, 2, 1)

    def test_review_1_interleaved_successful_library_saves(self, state_species, tmp_path, monkeypatch):
        ground, product = state_species[:2]
        ground.label, product.label = 'A', 'B'
        path = tmp_path / 'reactions.py'
        state_library(ground, ground, 1).save(str(path))
        first = state_library(ground, product, 2)
        second_product = copy.deepcopy(product)
        second_product.molecule[0].vibrational_level = 2
        second = state_library(ground, second_product, 3)
        paused, release = threading.Event(), threading.Event()
        errors = []
        original = os.replace
        def interleave(src, dst):
            if (threading.current_thread().name == 'save-A'
                    and Path(dst) == path.with_name('dictionary.txt') and str(src).endswith('.tmp')):
                paused.set()
                assert release.wait(20)
            return original(src, dst)
        monkeypatch.setattr(os, 'replace', interleave)
        def save_first():
            try:
                first.save(str(path))
            except BaseException as error:
                errors.append(error)
        child = threading.Thread(target=save_first, name='save-A')
        child.start()
        try:
            assert paused.wait(20)
            second.save(str(path))
        finally:
            release.set()
            child.join(20)
        assert not child.is_alive() and not errors
        assert_library_pair_safe(path, 3, 2)

    def test_review_2_chemkin_failed_dictionary_rename(self, state_species, tmp_path, monkeypatch):
        ground, product = state_species[:2]
        ground.label = 'A'
        product.label, product.index = 'N2', 2
        old = copy.deepcopy(product)
        old.label = 'N2(v1)'
        old.molecule[0].vibrational_level = -1
        path, verbose, dictionary = (tmp_path/name for name in ['chem.inp','annotated.inp','dictionary.txt'])
        save_chemkin(chemkin_model(ground, old, 1), str(path), str(verbose), str(dictionary))
        original = os.replace
        def fail(src, dst):
            if Path(dst) == dictionary and str(src).endswith('.tmp'):
                raise OSError('injected dictionary rename')
            return original(src, dst)
        monkeypatch.setattr(os, 'replace', fail)
        with pytest.raises(OSError, match='injected dictionary rename'):
            save_chemkin(chemkin_model(ground, product, 2), str(path), str(verbose), str(dictionary))
        try:
            _, reactions = load_chemkin_file(str(path), str(dictionary))
        except SpeciesIdentityError as error:
            assert str(path) in str(error) and str(dictionary) in str(error)
        else:
            reaction = reactions[0]
            expected = 1 if reaction.kinetics.A.value_si == 2 else -1
            assert reaction.products[0].molecule[0].vibrational_level == expected
        assert not list(tmp_path.glob('*.tmp')) and not list(tmp_path.glob('.*.tmp'))

    @pytest.mark.parametrize('writer', ['cantera1', 'cantera2'])
    def test_review_3_cantera_undeclared_product(self, state_species, tmp_path, writer):
        import cantera as ct
        ground, resolved = state_species[:2]
        ground.index = resolved.index = -1
        ground.label = 'N2(v1)' if writer == 'cantera1' else 'N2'
        resolved.label = 'N2'
        reaction = Reaction(reactants=[ground], products=[resolved],
                            kinetics=Arrhenius(A=(1, 's^-1'), Ea=(0, 'J/mol')))
        path = tmp_path / 'chem.yaml'
        model = ReactionModel(species=[ground], reactions=[reaction])
        try:
            if writer == 'cantera1':
                yaml_cantera1.write_cantera(model.species, model.reactions, model.get_elements(), path=str(path))
            else:
                yaml_cantera2.save_cantera_model(model, str(path))
        except SpeciesIdentityError:
            return
        gas = ct.Solution(str(path), transport_model=None)
        assert gas.n_species == 2, 'Undeclared excited product became a ground self-reaction'

    @pytest.mark.parametrize('writer', ['rms', 'cantera1', 'cantera2'])
    @pytest.mark.parametrize('kind', [ThirdBody, Lindemann, Troe])
    def test_review_4_undeclared_efficiency(self, state_species, writer, kind):
        ground, resolved = state_species[:2]
        reaction = Reaction(reactants=[ground], products=[ground],
                            kinetics=falloff(kind, {resolved.molecule[0]: 4}))
        with pytest.raises(SpeciesIdentityError, match='undeclared|identity'):
            if writer == 'rms':
                yaml_rms.get_mech_dict([ground], [reaction])
            elif writer == 'cantera1':
                yaml_cantera1.get_mech_dict_nonsurface([ground], [reaction])
            else:
                yaml_cantera2.generate_cantera_data([ground], [reaction],
                    elements_in_use=ReactionModel(species=[ground]).get_elements())

    @pytest.mark.parametrize('arrow', ['<=>', '=>'])
    def test_review_5_compact_equation(self, state_species, tmp_path, arrow):
        ground, resolved = state_species[:2]
        ground.label, resolved.label = 'A', 'BC'
        ground_copy = copy.deepcopy(resolved)
        ground_copy.label = 'C'
        ground_copy.molecule[0].vibrational_level = -1
        library = state_library(ground, resolved, 1)
        library.entries[1].label = 'A' + arrow + 'BC'
        library.entries[1].item.reversible = arrow == '<=>'
        # Keep C declared without creating a duplicate of A=>BC.
        other = state_library(ground_copy, ground_copy, 3).entries[1]
        other.index = 2
        library.entries[2] = other
        path = tmp_path / 'reactions.py'
        library.save(str(path))
        entry = loaded_reaction(path)
        assert entry.item.products[0].is_isomorphic(resolved)
        assert entry.item.products[0].label == 'BC'
        assert entry.item.reversible == (arrow == '<=>')

    def test_review_6_old_style_dictionary_refusal(self, state_species, tmp_path):
        path = tmp_path / 'dictionary.txt'
        path.write_text('sentinel')
        with pytest.raises(SpeciesIdentityError, match='old|state|resolved'):
            save_species_dictionary(str(path), [state_species[1]], old_style=True)
        assert path.read_text() == 'sentinel'

    def test_review_7_stale_library_label(self, state_species, tmp_path):
        ground, resolved = state_species[:2]
        ground.label, resolved.label = 'A', 'B'
        library = state_library(ground, resolved, 1)
        library.entries[1].label = 'A <=> A'
        path = tmp_path / 'reactions.py'
        with pytest.raises(SpeciesIdentityError, match='label|equation'):
            library.save(str(path))
        assert not path.exists()

    def test_review_nice_1_new_permissions_follow_umask(self, state_species, tmp_path):
        path = tmp_path / 'reactions.py'
        previous = os.umask(0o027)
        try:
            library_with(*state_species[:2]).save(str(path))
        finally:
            os.umask(previous)
        assert path.stat().st_mode & 0o777 == 0o640
        assert path.with_name('dictionary.txt').stat().st_mode & 0o777 == 0o640

    def test_review_nice_2_same_destination_refusal(self, state_species, tmp_path):
        path = tmp_path / 'dictionary.txt'
        path.write_text('sentinel')
        with pytest.raises(SpeciesIdentityError, match='same|destination|alias'):
            library_with(*state_species[:2]).save(str(path))
        assert path.read_text() == 'sentinel'


CENSUS_WRITERS = ('chemkin-new', 'chemkin-old', 'rms', 'cantera1', 'cantera2', 'library', 'library-old')
CENSUS_POSITIONS = ('reactant', 'product', 'collider', 'ThirdBody', 'Lindemann', 'Troe',
                    'coverage', 'thermo-coverage')
CENSUS_SCENARIOS = ('declared', 'undeclared', 'shared-ground-label')


def census_fixture(state_species, position, scenario):
    """One state-bearing reference; the species inventory varies independently."""
    from rmgpy.kinetics import SurfaceArrhenius
    from rmgpy.quantity import Quantity
    ground, resolved = state_species[:2]
    ground.index = resolved.index = -1
    ground.label, resolved.label = 'ground', 'resolved'
    if position in ('coverage', 'thermo-coverage'):
        ground.molecule = [Molecule().from_adjacency_list('1 X u0 p0 c0')]
        resolved.molecule = [Molecule().from_adjacency_list('vibrationallevel 1\n1 X u0 p0 c0')]
    if scenario == 'shared-ground-label':
        resolved.label = ground.label
    reaction = Reaction(reactants=[ground], products=[ground],
                        kinetics=Arrhenius(A=(1, 's^-1'), Ea=(0, 'J/mol')))
    if position == 'reactant':
        reaction.reactants = [resolved]
    elif position == 'product':
        reaction.products = [resolved]
    elif position == 'collider':
        reaction.specific_collider = resolved
        reaction.kinetics = falloff()
    elif position in ('ThirdBody', 'Lindemann', 'Troe'):
        reaction.kinetics = falloff(globals()[position], {resolved.molecule[0]: 4})
    elif position == 'coverage':
        reaction.kinetics = SurfaceArrhenius(A=(1, 's^-1'), Ea=(0, 'J/mol'), coverage_dependence={
            resolved: {'a': Quantity(0), 'm': Quantity(0), 'E': Quantity(0, 'J/mol')}})
    elif position == 'thermo-coverage':
        ground.thermo.thermo_coverage_dependence = {resolved.molecule[0].to_adjacency_list(): {
            'model': 'polynomial',
            'enthalpy-coefficients': [(1, 'J/mol'), (2, 'J/mol')],
            'entropy-coefficients': [(1, 'J/(mol*K)'), (2, 'J/(mol*K)')]}}
    declarations = [ground] + ([resolved] if scenario != 'undeclared' else [])
    if position in ('coverage', 'thermo-coverage'):
        companion = copy.deepcopy(state_species[3])
        companion.label = 'Ar'
        companion.molecule[0].electronic_state = ''
        declarations.append(companion)
    return declarations, reaction, resolved


def assert_loaded_reference(reaction, position, resolved):
    if position in ('reactant', 'product'):
        references = getattr(reaction, position + 's')
    elif position == 'collider':
        references = [reaction.specific_collider]
    elif position in ('ThirdBody', 'Lindemann', 'Troe'):
        references = list(reaction.kinetics.efficiencies)
    elif position == 'coverage':
        references = list(reaction.kinetics.coverage_dependence)
    else:
        references = [Molecule().from_adjacency_list(adj) for spc in reaction.reactants
                      for adj in spc.thermo.thermo_coverage_dependence]
    assert any(resolved.is_isomorphic(reference) for reference in references)


def exercise_census(writer, declarations, reaction, resolved, position, tmp_path):
    """Use public save paths, then inspect a reloaded reference's full identity."""
    from rmgpy.quantity import Quantity
    if writer.startswith('library'):
        library = library_with(reaction.reactants[0], reaction.products[0],
                               reaction.specific_collider, reaction.kinetics)
        # A library owns its dictionary inventory, including coverage-only keys.
        # Include the other declared identities as auxiliary reactions.
        for index, spc in enumerate(declarations, 2):
            if any(spc is ref for ref in reaction.reactants + reaction.products + [reaction.specific_collider]):
                continue
            library.entries[index] = Entry(index=index, label=spc.label + ' <=> ' + spc.label,
                item=Reaction(reactants=[spc], products=[spc]),
                data=Arrhenius(A=(index, 's^-1'), Ea=(0, 'J/mol')))
        if position == 'thermo-coverage':
            # Libraries do not serialize species thermo. The actual coverage
            # reference must refuse rather than imply it survives the dictionary.
            library.save(str(tmp_path / 'reactions.py')) if writer == 'library' else library.save_old(str(tmp_path))
            raise AssertionError('Library silently omitted a resolved thermodynamic coverage reference')
        if writer == 'library-old':
            library.save_old(str(tmp_path))
            loaded = KineticsLibrary()
            loaded.load_old(str(tmp_path))
        else:
            path = tmp_path / 'reactions.py'
            library.save(str(path))
            loaded = KineticsLibrary()
            loaded.load(str(path), local_context={cls.__name__: cls for cls in
                (Arrhenius, ThirdBody, Lindemann, Troe, SurfaceArrhenius)})
        entry = next(iter(loaded.entries.values()))
        entry.item.kinetics = entry.data
        assert_loaded_reference(entry.item, position, resolved)
        return
    if writer.startswith('chemkin'):
        model = CoreEdgeReactionModel()
        model.core.species, model.core.reactions = declarations, [reaction]
        model.surface_site_density = Quantity(2.5e-5, 'mol/m^2')
        path, verbose, dictionary = (tmp_path/name for name in ['chem.inp','annotated.inp','dictionary.txt'])
        save_chemkin(model, str(path), str(verbose), str(dictionary), old_style_dictionary=writer == 'chemkin-old')
        surface = tmp_path / 'chem-surface.inp'
        _, reactions = load_chemkin_file(str(tmp_path / 'chem-gas.inp') if surface.exists() else str(path), str(dictionary),
                                        surface_path=str(surface) if surface.exists() else None)
        assert_loaded_reference(reactions[0], position, resolved)
        return
    if writer == 'rms':
        path = tmp_path / 'chem.rms'
        yaml_rms.write_rms(declarations, [reaction], path=str(path))
        data = yaml.safe_load(path.read_text())
        by_name = {spc.label: spc for spc in yaml_rms.load_rms_species(str(path))}
        output = data['Reactions'][0]
        names = output[position + 's'] if position in ('reactant', 'product') else list(output['kinetics']['efficiencies'])
        assert any(by_name[name].is_isomorphic(resolved) for name in names)
        return
    path = tmp_path / 'chem.yaml'
    elements = ReactionModel(species=declarations).get_elements()
    if writer == 'cantera1':
        yaml_cantera1.write_cantera(declarations, [reaction], elements, surface_site_density=2.5e-5, path=str(path))
    else:
        data = yaml_cantera2.generate_cantera_data(declarations, [reaction], elements_in_use=elements,
                                                 site_density=2.5e-5)
        path.write_text(yaml.dump(data, Dumper=yaml_cantera2.Dumper, sort_keys=False))
    # Reload the serialized Cantera metadata and schema references. Native
    # Solution loading is independently exercised by review reproduction 3.
    import re
    from rmgpy.util import parse_reaction_equation, get_reaction_collider
    data = yaml.safe_load(path.read_text())
    by_name = {}
    for record in data['species']:
        match = re.search(r'(?:electronicstate|vibrationallevel|multiplicity)\s', record.get('note', ''))
        if match:
            by_name[record['name']] = Species(molecule=[Molecule().from_adjacency_list(record['note'][match.start():])])
    output = next(record for key, records in data.items() if key.endswith('reactions')
                  for record in records if isinstance(record, dict) and 'equation' in record)
    if position in ('reactant', 'product'):
        left, right, _ = parse_reaction_equation(output['equation'])
        names = (left if position == 'reactant' else right).split(' + ')
    elif position == 'collider':
        left, _, _ = parse_reaction_equation(output['equation'])
        collider = get_reaction_collider(left)
        names = [collider[2:-1]] if collider else [left.split(' + ')[-1]]
    elif position in ('ThirdBody', 'Lindemann', 'Troe'):
        names = list(output['efficiencies'])
    elif position == 'coverage':
        names = list(output['coverage-dependencies'])
    else:
        names = [name for record in data['species'] for name in record.get('coverage-dependencies', {})]
    assert any(name in by_name and by_name[name].is_isomorphic(resolved) for name in names)



class TestExcitedExportReferenceCensus:
    @pytest.mark.parametrize('writer', CENSUS_WRITERS)
    @pytest.mark.parametrize('position', CENSUS_POSITIONS)
    @pytest.mark.parametrize('scenario', CENSUS_SCENARIOS)
    def test_every_reference(self, state_species, tmp_path, writer, position, scenario, record_property):
        declarations, reaction, resolved = census_fixture(state_species, position, scenario)
        try:
            exercise_census(writer, declarations, reaction, resolved, position, tmp_path)
        except SpeciesIdentityError as error:
            assert str(error)
            record_property('outcome', 'named refusal: ' + type(error).__name__)
        else:
            record_property('outcome', 'identity/reference roundtrip')


def writer_routing_violations(module, source):
    """Audit writer ASTs: raw names may only allocate declarations or diagnostics.

    Shared contexts own lookup. Primitive formatters allocate identifiers, and
    note strings and errors retain legacy labels, but cannot emit schema names.
    Chemkin's few Cython declarations/types are stripped before the same audit.
    """
    import ast
    import re
    if module == 'chemkin.pyx':
        source = source.replace(' cimport ', ' import ').replace('cpdef ', 'def ')
        source = re.sub(r'^\s*cdef .*$', '', source, flags=re.M)
        source = re.sub(r'<[A-Z]\w*>', '', source)
        source = source.replace('list reaction_list', 'reaction_list').replace('NASAPolynomial poly', 'poly')
        source = source.replace('bint verbose=True', 'verbose=True')
    tree = ast.parse(source)
    parents = {child: parent for parent in ast.walk(tree) for child in ast.iter_child_nodes(parent)}
    violations = []
    selected = {
        'data/base.py': {'render_dictionary'},
        'data/kinetics/library.py': {'save', 'save_old', 'save_entry'},
        'data/kinetics/common.py': {'save_entry', 'library_reaction_equation', 'library_serializable_kinetics'},
    }
    for function in (node for node in ast.walk(tree) if isinstance(node, ast.FunctionDef)):
        if module in selected and function.name not in selected[module]:
            continue
        if module == 'chemkin.pyx' and not function.name.startswith(('write_', 'render_', 'save_', 'validate_', '_chemkin_')):
            continue
        if function.name in ('_cantera_identifier', 'load_rms_species', 'admit_rms_species'):
            continue
        for node in ast.walk(function):
            if isinstance(node, ast.Call):
                callee = node.func
                if ((isinstance(callee, ast.Name) and callee.id == 'get_species_identifier')
                        or (isinstance(callee, ast.Attribute) and callee.attr == 'to_chemkin')
                        or (isinstance(callee, ast.Attribute) and callee.attr == 'index'
                            and ast.unparse(callee.value) in ('spcs', 'species', 'species_list', 'declarations'))):
                    # Legacy ground coverage identifiers are allocated once in
                    # SpeciesReferences; emitting lookups must still use the resolver.
                    ancestor = node
                    allocation = False
                    while ancestor in parents:
                        ancestor = parents[ancestor]
                        if (module == 'data/kinetics/common.py'
                                and function.name == 'library_serializable_kinetics'
                                and isinstance(ancestor, ast.Call)
                                and isinstance(ancestor.func, ast.Name)
                                and ancestor.func.id == 'SpeciesReferences'):
                            allocation = True
                            break
                    if allocation:
                        continue
                    violations.append((function.name, node.lineno, ast.unparse(node)))
            if not isinstance(node, ast.Attribute) or node.attr != 'label' or not isinstance(node.ctx, ast.Load):
                continue
            if isinstance(node.value, ast.Name) and node.value.id in ('entry', 'self'):
                continue  # Entry equations are checked against reaction objects; self.label names a database.
            ancestors = []
            ancestor = node
            while ancestor in parents:
                ancestor = parents[ancestor]
                ancestors.append(ancestor)
                if ancestor is function:
                    break
            if any(isinstance(item, ast.Raise) for item in ancestors):
                continue
            if any(isinstance(item, (ast.Assign, ast.AugAssign, ast.Expr)) and
                   'Specific third body collider:' in ast.unparse(item) for item in ancestors):
                continue  # Non-semantic provenance note, unchanged for ground exports.
            if module == 'yaml_rms.py' and function.name == 'get_mech_dict' and any(
                    isinstance(item, ast.Assign) and ast.unparse(item) == 'names = [x.label for x in spcs]'
                    for item in ancestors):
                continue  # Initial declaration allocation, never reference lookup.
            if module == 'data/kinetics/library.py' and any(isinstance(item, ast.Lambda) for item in ancestors):
                continue  # Sorting dictionary declarations does not look up references.
            if module == 'data/kinetics/common.py' and function.name == 'save_entry' and any(
                    isinstance(item, ast.If) and "collider.label.strip().upper() == 'M'" in ast.unparse(item.test)
                    for item in ancestors):
                continue  # Reserved-name refusal guard.
            violations.append((function.name, node.lineno, ast.unparse(node)))
    return violations


class TestExcitedExportRouting:
    @pytest.mark.parametrize('module', ['chemkin.pyx', 'yaml_rms.py', 'yaml_cantera1.py', 'yaml_cantera2.py',
        'data/base.py', 'data/kinetics/common.py', 'data/kinetics/library.py'])
    def test_writer_has_no_private_name_lookup(self, module):
        import rmgpy
        source = (Path(rmgpy.__file__).parent / module).read_text()
        assert not writer_routing_violations(module, source)

    @pytest.mark.parametrize('module,function', [('yaml_cantera2.py', 'get_label'), ('yaml_cantera1.py', 'species_to_dict'),
        ('yaml_rms.py', 'obj_to_dict'), ('chemkin.pyx', 'write_reaction_string'),
        ('data/kinetics/common.py', 'library_reaction_equation'), ('data/kinetics/library.py', 'save_old'),
        ('data/base.py', 'render_dictionary')])
    def test_static_guard_rejects_private_lookup_mutation(self, module, function):
        source = 'def ' + function + '(spc, spcs):\n    return spcs[spcs.index(spc)].label\n'
        assert writer_routing_violations(module, source)


class TestExcitedExportGeneration:
    @pytest.mark.parametrize('writer', ['library', 'chemkin', 'library-old'])
    @pytest.mark.parametrize('change', ['different', 'kinetics-only', 'dictionary-only', 'malformed'])
    def test_split_generation_names_both_paths_and_tokens(self, state_species, tmp_path, writer, change):
        from rmgpy.exceptions import GenerationMismatchError
        from rmgpy.util import generation_token, strip_generation_marker
        ground, resolved = state_species[:2]
        if writer == 'library':
            kinetics, dictionary = tmp_path/'reactions.py', tmp_path/'dictionary.txt'
            library_with(ground, resolved).save(str(kinetics))
            load = lambda: loaded_reaction(kinetics)
        elif writer == 'library-old':
            kinetics, dictionary = tmp_path/'reactions.txt', tmp_path/'species.txt'
            library_with(ground, ground).save_old(str(tmp_path))
            load = lambda: KineticsLibrary().load_old(str(tmp_path))
        else:
            kinetics, dictionary = tmp_path/'chem.inp', tmp_path/'dictionary.txt'
            save_chemkin(chemkin_model(ground, resolved, 1), str(kinetics), str(tmp_path/'annotated.inp'), str(dictionary))
            load = lambda: load_chemkin_file(str(kinetics), str(dictionary))
        first, second = generation_token(kinetics.read_text()), generation_token(dictionary.read_text())
        assert first == second and first.startswith('v1-sha256:')
        if change == 'different':
            dictionary.write_text(dictionary.read_text().replace(second, 'v1-sha256:' + '0'*64))
        elif change == 'kinetics-only':
            dictionary.write_text(strip_generation_marker(dictionary.read_text()))
        elif change == 'dictionary-only':
            kinetics.write_text(strip_generation_marker(kinetics.read_text()))
        else:
            kinetics.write_text(kinetics.read_text().replace(first, 'v1-sha256:broken'))
        with pytest.raises(GenerationMismatchError) as captured:
            load()
        message = str(captured.value)
        assert str(kinetics) in message and str(dictionary) in message
        assert repr(generation_token(kinetics.read_text())) in message
        assert repr(generation_token(dictionary.read_text())) in message

    @pytest.mark.parametrize('writer', ['library', 'chemkin', 'library-old'])
    def test_token_free_pair_loads_without_changing_identity(self, state_species, tmp_path, writer):
        from rmgpy.util import strip_generation_marker
        ground, resolved = state_species[:2]
        if writer == 'library':
            path = tmp_path/'reactions.py'
            library_with(ground, resolved).save(str(path))
        elif writer == 'library-old':
            library_with(ground, ground).save_old(str(tmp_path))
        else:
            path = tmp_path/'chem.inp'
            save_chemkin(chemkin_model(ground, resolved, 1), str(path), str(tmp_path/'annotated.inp'), str(tmp_path/'dictionary.txt'))
        for path in tmp_path.iterdir():
            path.write_text(strip_generation_marker(path.read_text()))
        if writer == 'library':
            assert loaded_reaction(tmp_path/'reactions.py').item.products[0].is_isomorphic(resolved)
        elif writer == 'library-old':
            # The existing legacy reader rejects save_old's Unit/Reactions:
            # header; generation support must not change token-free parsing.
            from rmgpy.exceptions import ChemkinError
            from rmgpy.chemkin import read_reactions_block
            import io
            with pytest.raises(ChemkinError, match='Invalid reaction block'):
                read_reactions_block(io.StringIO((tmp_path/'reactions.txt').read_text()),
                                     {ground.label: ground})
            with pytest.raises(ChemkinError, match='Invalid reaction block'):
                KineticsLibrary().load_old(str(tmp_path))
        else:
            _, reactions = load_chemkin_file(str(tmp_path/'chem.inp'), str(tmp_path/'dictionary.txt'))
            assert reactions[0].products[0].is_isomorphic(resolved)

    @pytest.mark.parametrize('arrow', ['<=>', '=>', '='])
    @pytest.mark.parametrize('spacing', ['', ' '])
    def test_all_arrow_forms_roundtrip(self, state_species, tmp_path, arrow, spacing):
        ground, resolved = state_species[:2]
        ground.label, resolved.label = 'A', 'BC'
        library = library_with(ground, resolved)
        library.entries[1].item.reversible = arrow != '=>'
        library.entries[1].label = 'A' + spacing + arrow + spacing + 'BC'
        path = tmp_path/'reactions.py'
        library.save(str(path))
        # Exercise the loader's compact grammar too, rather than only save's canonical equation.
        from rmgpy.util import strip_generation_marker
        for output in tmp_path.iterdir():
            output.write_text(strip_generation_marker(output.read_text()))
        path.write_text(path.read_text().replace('A <=> BC' if arrow != '=>' else 'A => BC', library.entries[1].label))
        entry = loaded_reaction(path)
        assert entry.item.products[0].is_isomorphic(resolved)
        assert entry.item.reversible == (arrow != '=>')

    @pytest.mark.parametrize('alias', ['same', 'symlink', 'hardlink'])
    def test_pair_destination_alias_refuses_before_replacement(self, state_species, tmp_path, alias):
        path = tmp_path/'reactions.py'
        dictionary = tmp_path/'dictionary.txt'
        path.write_text('sentinel')
        if alias == 'same':
            path = dictionary
            path.write_text('sentinel')
        elif alias == 'symlink':
            dictionary.symlink_to(path)
        else:
            os.link(path, dictionary)
        with pytest.raises(SpeciesIdentityError, match='alias|same'):
            library_with(*state_species[:2]).save(str(path))
        assert path.read_text() == 'sentinel'

    def test_preserves_existing_modes_and_content_digest_is_deterministic(self, state_species, tmp_path):
        from rmgpy.util import generation_token
        path, dictionary = tmp_path/'reactions.py', tmp_path/'dictionary.txt'
        path.write_text('old kinetics'); dictionary.write_text('old dictionary')
        path.chmod(0o754); dictionary.chmod(0o640)
        library = library_with(*state_species[:2])
        library.save(str(path))
        token = generation_token(path.read_text())
        assert path.stat().st_mode & 0o777 == 0o754
        assert dictionary.stat().st_mode & 0o777 == 0o640
        library.save(str(path))
        assert generation_token(path.read_text()) == token

    @pytest.mark.parametrize('writer', ['library', 'chemkin'])
    def test_loader_parses_the_validated_snapshots(self, state_species, tmp_path, monkeypatch, writer):
        from rmgpy import util
        ground, resolved = state_species[:2]
        ground.label, resolved.label = 'A', 'B'
        later = copy.deepcopy(resolved)
        later.molecule[0].vibrational_level = 2
        if writer == 'library':
            path = tmp_path/'reactions.py'
            state_library(ground, resolved, 1).save(str(path))
            original = util.read_generation_pair
            def change_after_validation(*args):
                snapshots = original(*args)
                state_library(ground, later, 2).save(str(path))
                return snapshots
            monkeypatch.setattr(util, 'read_generation_pair', change_after_validation)
            entry = loaded_reaction(path)
            assert entry.data.A.value_si == 1 and entry.item.products[0].is_isomorphic(resolved)
        else:
            path, dictionary = tmp_path/'chem.inp', tmp_path/'dictionary.txt'
            verbose = tmp_path/'annotated.inp'
            save_chemkin(chemkin_model(ground, resolved, 1), str(path), str(verbose), str(dictionary))
            original = util.read_generation_files
            def change_after_validation(*args):
                snapshots = original(*args)
                save_chemkin(chemkin_model(ground, later, 2), str(path), str(verbose), str(dictionary))
                return snapshots
            monkeypatch.setattr(util, 'read_generation_files', change_after_validation)
            _, reactions = load_chemkin_file(str(path), str(dictionary))
            assert reactions[0].kinetics.A.value_si == 1 and reactions[0].products[0].is_isomorphic(resolved)

    @pytest.mark.parametrize('failure', ['crash', 'interleave'])
    def test_chemkin_split_generation_never_silently_changes_state(self, state_species, tmp_path, monkeypatch, failure):
        ground, resolved = state_species[:2]
        ground.label, resolved.label = 'A', 'B'
        later = copy.deepcopy(resolved)
        later.molecule[0].vibrational_level = 2
        path, dictionary, verbose = (tmp_path/name for name in ['chem.inp', 'dictionary.txt', 'annotated.inp'])
        def save(product, rate):
            save_chemkin(chemkin_model(ground, product, rate), str(path), str(verbose), str(dictionary))
        save(ground, 1)
        original = os.replace
        if failure == 'crash':
            def crash_save():
                def crash(src, dst):
                    if Path(dst) == dictionary and str(src).endswith('.tmp'):
                        os._exit(86)
                    return original(src, dst)
                os.replace = crash
                save(resolved, 2)
            child = multiprocessing.get_context('fork').Process(target=crash_save)
            child.start(); child.join(20)
            assert not child.is_alive() and child.exitcode == 86
            expected_rate, expected = 2, resolved
        else:
            paused, release = threading.Event(), threading.Event()
            errors = []
            def interleave(src, dst):
                if threading.current_thread().name == 'save-A' and Path(dst) == dictionary and str(src).endswith('.tmp'):
                    paused.set()
                    assert release.wait(20)
                return original(src, dst)
            monkeypatch.setattr(os, 'replace', interleave)
            def save_first():
                try:
                    save(resolved, 2)
                except BaseException as error:
                    errors.append(error)
            child = threading.Thread(target=save_first, name='save-A')
            child.start()
            try:
                assert paused.wait(20)
                save(later, 3)
            finally:
                release.set(); child.join(20)
            assert not child.is_alive() and not errors
            expected_rate, expected = 3, later
        try:
            _, reactions = load_chemkin_file(str(path), str(dictionary))
        except SpeciesIdentityError as error:
            assert str(path) in str(error) and str(dictionary) in str(error)
        else:
            assert reactions[0].kinetics.A.value_si == expected_rate
            assert reactions[0].products[0].is_isomorphic(expected)

    def test_database_library_save_keeps_the_generation_pair(self, state_species, tmp_path):
        from rmgpy.data.kinetics.database import KineticsDatabase
        database = KineticsDatabase()
        database.libraries['state'] = library_with(*state_species[:2])
        database.save_libraries(str(tmp_path))
        entry = loaded_reaction(tmp_path/'state/reactions.py')
        assert entry.item.products[0].is_isomorphic(state_species[1])
