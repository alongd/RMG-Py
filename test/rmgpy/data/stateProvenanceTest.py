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

"""Indirect data provenance and consistent resolved-state training refusal."""
from copy import deepcopy
import pickle
import pytest

from rmgpy.data.base import Entry
from rmgpy.data.kinetics.family import KineticsFamily
from rmgpy.data.kinetics.groups import KineticsGroups
from rmgpy.data.kinetics.rules import KineticsRules
from rmgpy.data.thermo import ThermoDatabase, ThermoGroups
from rmgpy.data.transport import TransportDatabase, TransportGroups, CriticalPointGroupContribution
from rmgpy.data.solvation import SolvationDatabase, SoluteGroups, SoluteData
from rmgpy.kinetics import ArrheniusEP
from rmgpy.molecule import Molecule, Group
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import ThermoData

STATES = [('A', -1), ('', 1), ('A', 1)]
WILDCARD = 'electronicstate x\nvibrationallevel x\n'


def datum(kind, value=12):
    if kind == 'thermo':
        return ThermoData(Tdata=([300,400,500,600,800,1000,1500], 'K'),
                          Cpdata=([1]*7, 'J/(mol*K)'), H298=(value, 'J/mol'), S298=(12, 'J/(mol*K)'))
    if kind == 'solute':
        return SoluteData(S=value, B=0, E=0, L=0, A=0)
    return CriticalPointGroupContribution(Tc=value, Pc=0, Vc=0, Tb=0, structureIndex=0)


def tree(kind, structure='1 * C u0', root_data=None, child_data=None):
    cls = {'thermo': ThermoGroups, 'solute': SoluteGroups, 'transport': TransportGroups, 'kinetics': KineticsGroups}[kind]
    db = cls(label=kind)
    root = Entry(label='Any', item=Group().from_adjacency_list(WILDCARD + structure), data=root_data)
    child = Entry(label='Ground', item=Group().from_adjacency_list(structure), data=child_data, parent=root, nodal_distance=1)
    root.children = [child]
    db.entries = {e.label: e for e in (root, child)}
    db.top = [root]
    db.generic_nodes = []
    return db, root, child


def molecule(state, smiles='C'):
    mol = Molecule(smiles=smiles, electronic_state=state[0], vibrational_level=state[1])
    mol.atoms[0].label = '*'
    return mol


def refusal(call, name='StateProvenanceError'):
    with pytest.raises(Exception) as err:
        call()
    assert type(err.value).__name__ == name, str(err.value)


def group_lookup(kind, db, mol, remove=False):
    atoms = {'*': mol.atoms[0]}
    if kind == 'thermo':
        owner = ThermoDatabase()
        if remove:
            return owner._remove_group_thermo_data(datum(kind, 0), db, mol, atoms)
        return owner._add_group_thermo_data(None, db, mol, atoms)[0]
    if kind == 'solute':
        owner = SolvationDatabase()
        if remove:
            return owner._remove_group_solute_data(datum(kind, 0), db, mol, atoms)
        return owner._add_group_solute_data(None, db, mol, atoms)
    return TransportDatabase()._add_critical_point_contribution(datum(kind, 0), db, mol, atoms)


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('kind,remove', [('thermo',False), ('thermo',True), ('solute',False), ('solute',True), ('transport',False)])
def test_alias_checks_target(kind, remove, state):
    db, root, child = tree(kind, root_data='Ground', child_data=datum(kind, 12345 if kind == 'thermo' else 12))
    mol = molecule(state)
    assert db.descend_tree(mol, {'*': mol.atoms[0]}) is root
    assert not db.match_node_to_structure(child, mol, {'*': mol.atoms[0]})
    refusal(lambda: group_lookup(kind, db, mol, remove), name='ExcitedSpeciesThermoError' if kind in ('thermo', 'solute') else 'StateProvenanceError')


@pytest.mark.parametrize('kind', ['thermo','solute','transport'])
@pytest.mark.parametrize('state', STATES)
def test_matched_alias_and_ancestor_data_allowed(kind, state):
    db, root, child = tree(kind, root_data='Ground', child_data=datum(kind))
    child.item.electronic_state = ['x']; child.item.vibrational_level = ['x']
    mol = molecule(state)
    if kind == 'transport':
        group_lookup(kind, db, mol)
    else:
        refusal(lambda: group_lookup(kind, db, mol), 'ExcitedSpeciesThermoError')
    root.data = datum(kind); child.data = None
    if kind == 'transport':
        group_lookup(kind, db, mol)
    else:
        refusal(lambda: group_lookup(kind, db, mol), 'ExcitedSpeciesThermoError')


@pytest.mark.parametrize('kind', ['thermo','solute','transport'])
@pytest.mark.parametrize('state', STATES)
def test_unmatched_ancestor_refused(kind, state):
    db, root, child = tree(kind, root_data=datum(kind))
    # Starting at an explicitly supplied descendant can skip matching ancestors.
    root.item.electronic_state = []; root.item.vibrational_level = []
    child.item.electronic_state = ['x']; child.item.vibrational_level = ['x']
    db.top = [child]
    refusal(lambda: group_lookup(kind, db, molecule(state)), name='ExcitedSpeciesThermoError' if kind in ('thermo', 'solute') else 'StateProvenanceError')


@pytest.mark.parametrize('kind', ['thermo','solute'])
@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('warm', [False, True])
def test_ring_average_provenance(kind, state, warm):
    mol = molecule(state, 'C1CC1')
    structure = '1 * C u0 {2,S} {3,S}\n2 C u0 {1,S} {3,S}\n3 C u0 {1,S} {2,S}'
    db, root, child = tree(kind, structure, child_data=datum(kind))
    owner = ThermoDatabase() if kind == 'thermo' else SolvationDatabase()
    lookup = getattr(owner, '_add_ring_correction_' + ('thermo' if kind == 'thermo' else 'solute') + '_data_from_tree')
    ring = mol.get_smallest_set_of_smallest_rings()[0]
    if warm:
        ordinary = molecule(('', -1), 'C1CC1')
        # Ensure averaging is cached on the wildcard node, not its data-bearing child.
        child.item = Group().from_adjacency_list('1 * C u0 {2,S}\n2 O u0 {1,S}')
        lookup(None, db, ordinary, ordinary.get_smallest_set_of_smallest_rings()[0])
    refusal(lambda: lookup(None, db, mol, ring), name='ExcitedSpeciesThermoError')


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('kind', ['thermo','solute'])
def test_ring_alias_target_checks(kind, state):
    mol = molecule(state, 'C1CC1')
    structure = '1 * C u0 {2,S} {3,S}\n2 C u0 {1,S} {3,S}\n3 C u0 {1,S} {2,S}'
    db, _, _ = tree(kind, structure, root_data='Ground', child_data=datum(kind))
    owner = ThermoDatabase() if kind == 'thermo' else SolvationDatabase()
    lookup = getattr(owner, '_add_ring_correction_' + ('thermo' if kind == 'thermo' else 'solute') + '_data_from_tree')
    refusal(lambda: lookup(None, db, mol, mol.get_smallest_set_of_smallest_rings()[0]), name='ExcitedSpeciesThermoError')


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('kind', ['thermo','solute'])
def test_polycyclic_decomposition_refused(kind, state):
    owner = ThermoDatabase() if kind == 'thermo' else SolvationDatabase()
    groups = {}
    for label in ['group','ring','polycyclic']:
        db, root, child = tree(kind, '1 * R u0', root_data=datum(kind,10), child_data=datum(kind,99))
        db.label=label; groups[label]=db
    groups['group'].top[0].data = datum(kind,1)
    cls = ThermoGroups if kind == 'thermo' else SoluteGroups
    for label in ['other','longDistanceInteraction_cyclic','longDistanceInteraction_noncyclic']:
        groups[label] = cls(label=label)
    owner.groups = groups
    mol = molecule(state, 'c1ccc2ccccc2c1')
    refusal(lambda: getattr(owner, 'compute_group_additivity_' + kind)(mol), name='ExcitedSpeciesThermoError')


def kinetic_fixture(state, exact=False):
    db, root, child = tree('kinetics')
    rules=KineticsRules(label='probe/rules')
    rules.entries={'Ground':[Entry(label='Ground', item=[child], data=ArrheniusEP(A=(123,'s^-1'), n=0, alpha=0, E0=(0,'J/mol')), rank=1)]}
    if exact:
        rules.entries['Any']=[Entry(label='Any', item=[root], data=ArrheniusEP(A=(321,'s^-1'), n=0, alpha=0, E0=(0,'J/mol')), rank=1)]
    family=KineticsFamily(label='probe',allow_excited_reactants=True);family.groups=db;family.rules=rules
    return family, root, child, Reaction(reactants=[molecule(state)], products=[])


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('verbose', [False,True])
def test_rate_average_provenance(state, verbose):
    family, root, child, rxn = kinetic_fixture(state)
    family.rules.fill_rules_by_averaging_up([root], {}, verbose=verbose)
    template=family.get_reaction_template(rxn)
    refusal(lambda: family.get_kinetics_for_template(template))
    restored=pickle.loads(pickle.dumps(family.rules))
    refusal(lambda: restored.estimate_kinetics(template))


@pytest.mark.parametrize('state', STATES)
def test_exact_rate_allowed(state):
    family, root, child, rxn=kinetic_fixture(state, exact=True)
    family.rules.fill_rules_by_averaging_up([root], {}, verbose=True)
    k,_=family.get_kinetics_for_template(family.get_reaction_template(rxn))
    assert k.A.value_si == 321


def test_unresolved_wildcard_average_unchanged():
    family, root, child, rxn=kinetic_fixture(('',-1))
    child.item=Group().from_adjacency_list('1 * O u0')
    family.rules.fill_rules_by_averaging_up([root], {}, verbose=False)
    k,_=family.get_kinetics_for_template(family.get_reaction_template(rxn))
    assert k.A.value_si == 123


@pytest.mark.parametrize('state', STATES)
def test_rate_ancestor_provenance(state):
    family, root, child, rxn=kinetic_fixture(state)
    child.item.electronic_state=['x'];child.item.vibrational_level=['x']
    root.item.electronic_state=[];root.item.vibrational_level=[]
    family.groups.top=[child]
    family.rules.entries={'Any':[Entry(label='Any',item=[root],data=family.rules.entries['Ground'][0].data,rank=1)]}
    refusal(lambda: family.get_kinetics_for_template(family.get_reaction_template(rxn)))


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('route', ['thermo_hbi','solute_hbi','solute_halogen','transport_radical','transport_fallback','surface_desorption','adsorption_groups','poly_thermo','poly_solute'])
def test_derived_molecule_routes_refuse(route,state):
    mol=molecule(state, '[CH3]' if 'hbi' in route or 'radical' in route else 'CCl' if 'halogen' in route else 'C1CC2CCC1C2')
    td=ThermoDatabase();sd=SolvationDatabase();tr=TransportDatabase()
    calls={
        'thermo_hbi':lambda:td.estimate_radical_thermo_via_hbi(mol, lambda sp:datum('thermo')),
        'solute_hbi':lambda:sd.estimate_radical_solute_data_via_hbi(mol, lambda sp:datum('solute')),
        'solute_halogen':lambda:sd.estimate_halogen_solute_data(mol, lambda sp:datum('solute')),
        'transport_radical':lambda:tr.estimate_critical_properties_via_group_additivity(mol),
        'transport_fallback':lambda:tr.get_transport_properties_via_lennard_jones_parameters(Species(molecule=[mol])),
        'surface_desorption':lambda:td.get_thermo_data_for_surface_species(Species(molecule=[mol])),
        'adsorption_groups':lambda:td._add_adsorption_correction(datum('thermo'),None,mol,[]),
        'poly_thermo':lambda:td._add_polycyclic_correction_thermo_data(datum('thermo'),mol,mol.get_disparate_cycles()[1][0]),
        'poly_solute':lambda:sd._add_polycyclic_correction_solute_data(datum('solute'),mol,mol.get_disparate_cycles()[1][0]),
    }
    if route=='adsorption_groups':td.adsorption_groups='adsorptionPt111'
    refusal(calls[route])


@pytest.mark.parametrize('state', STATES)
def test_thermo_copy_and_removal_provenance(state):
    db, root, child=tree('thermo', child_data=datum('thermo'))
    db.copy_data(child,root)
    refusal(lambda:group_lookup('thermo',db,molecule(state)), name='ExcitedSpeciesThermoError')
    root.data='Ground';db.remove_group(child)
    refusal(lambda:group_lookup('thermo',db,molecule(state)), name='ExcitedSpeciesThermoError')


def training_fixture(state):
    mols=[molecule(state, '[CH3]') for _ in range(2)]
    for i,mol in enumerate(mols,1):mol.atoms[0].label='*'+str(i)
    rxn=Reaction(reactants=[Species(molecule=[m]) for m in mols], products=[Species(smiles='CC')])
    family=KineticsFamily(label='Training',allow_excited_reactants=True)
    roots=[Entry(index=i,label='Root'+str(i),item=Group().from_adjacency_list(WILDCARD+'1 *'+str(i)+' C u1')) for i in [1,2]]
    for e in roots:
        e.item.electronic_state=[state[0]] if state[0] else []
        e.item.vibrational_level=[state[1]] if state[1]>=0 else []
    family.groups=KineticsGroups();family.groups.top=roots;family.groups.entries={e.label:e for e in roots}
    family.forward_template=Reaction(reactants=roots);family.rules=KineticsRules()
    return family,roots,rxn


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('route', ['get_reaction_matches','rxns_match_node','generate_tree','save_training_reactions','get_rxn_batches','prune_tree','make_tree_nodes','make_bm_rules_from_template_rxn_map','cross_validate','eval_ext','get_extension_edge','extend_node'])
def test_training_entry_points_named_refusal(state,route):
    f, roots, rxn=training_fixture(state)
    mapping={roots[0].label:[rxn]}
    calls={
        'get_reaction_matches':lambda:f.get_reaction_matches(rxns=[rxn],estimate_thermo=False),
        'rxns_match_node':lambda:f.rxns_match_node(roots[0],[rxn]),
        'generate_tree':lambda:f.generate_tree(rxns=[rxn]),
        'save_training_reactions':lambda:f.save_training_reactions([rxn]),
        'get_rxn_batches':lambda:f.get_rxn_batches([rxn]),
        'prune_tree':lambda:f.prune_tree([], [rxn]),
        'make_tree_nodes':lambda:f.make_tree_nodes(template_rxn_map=mapping),
        'make_bm_rules_from_template_rxn_map':lambda:f.make_bm_rules_from_template_rxn_map(mapping),
        'cross_validate':lambda:f.cross_validate(template_rxn_map=mapping),
        'eval_ext':lambda:f.eval_ext(roots[0],roots[0].item,'test',mapping),
        'get_extension_edge':lambda:f.get_extension_edge(roots[0],mapping,None,1000),
        'extend_node':lambda:f.extend_node(roots[0],mapping),
        '_split_reactions':lambda:f._split_reactions([rxn],roots[0].item),
    }
    refusal(calls[route], 'ResolvedStateTrainingError')


@pytest.mark.parametrize('state', STATES)
def test_clean_tree_keeps_constraints(state):
    f, roots, rxn=training_fixture(state)
    f.clean_tree_groups()
    merged=f.groups.top[0].item
    assert merged.electronic_state == ([state[0]] if state[0] else [])
    assert merged.vibrational_level == ([state[1]] if state[1]>=0 else [])
    assert f.reaction_matches(rxn,merged)


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('route', ['get_training_set', 'add_rules_from_training', 'cross_validate_old', 'regularize', 'make_tree'])
def test_training_depository_entry_points_named_refusal(state, route):
    from rmgpy.data.kinetics.depository import KineticsDepository
    family, roots, rxn = training_fixture(state)
    dep = KineticsDepository(label='Training/training')
    dep.entries = {1: Entry(index=1, item=rxn)}
    family.depositories = [dep]
    refusal(lambda: getattr(family, route)(), 'ResolvedStateTrainingError')


@pytest.mark.parametrize('state', STATES)
def test_auto_generated_rates_refuse_unverified_state(state):
    family, root, child, rxn = kinetic_fixture(state, exact=True)
    family.rules.auto_generated = True
    refusal(lambda: family.get_kinetics_for_template(family.get_reaction_template(rxn)))


@pytest.mark.parametrize('state', STATES)
def test_multistep_alias_and_cycle_refusal(state):
    db, root, child = tree('thermo', root_data='Bridge', child_data=datum('thermo'))
    bridge = Entry(label='Bridge',item=root.item.copy(deep=True),data='Ground')
    db.entries['Bridge'] = bridge
    refusal(lambda: group_lookup('thermo',db,molecule(state)), name='ExcitedSpeciesThermoError')
    bridge.data = 'Any'
    refusal(lambda: group_lookup('thermo',db,molecule(state)), name='ExcitedSpeciesThermoError')


@pytest.mark.parametrize('field,left,right', [('electronic_state',['A'],['B']),('vibrational_level',[1],[2])])
def test_tree_merge_refuses_incompatible_states(field,left,right):
    first=Group().from_adjacency_list('1 *1 C u1')
    second=Group().from_adjacency_list('1 *2 C u1')
    setattr(first,field,left);setattr(second,field,right)
    refusal(lambda:first.merge_groups(second), 'StateConstraintMergeError')


@pytest.mark.parametrize('state', STATES)
def test_sidt_group_tree_refuses_unverified_root_data(state):
    from pysidt import MultiTargetSingleEvalSubgraphIsomorphicDecisionTree
    from pysidt.sidt import Node
    from types import SimpleNamespace
    owner=ThermoDatabase();owner.adsorption_groups='SIDT'
    node=Node(name='Root',group=Group().from_adjacency_list('1 * R ux'),
              rule=SimpleNamespace(value=[12345]*9,uncertainty=[0]*9))
    tree=MultiTargetSingleEvalSubgraphIsomorphicDecisionTree(nodes={'Root':node})
    owner.sidts={'Pt111_monodentate_adsorption_corrections':tree}
    mol=molecule(state,'[CH2]=*')
    refusal(lambda:owner._add_adsorption_correction(datum('thermo'),None,mol,mol.get_surface_sites()))


@pytest.mark.parametrize('kind', ['thermo','solute'])
@pytest.mark.parametrize('legacy', [False,True])
def test_text_save_refuses_lost_average_provenance(kind,legacy,tmp_path):
    from io import StringIO
    db,root,child=tree(kind,child_data=datum(kind))
    owner=ThermoDatabase() if kind=='thermo' else SolvationDatabase()
    success,data=getattr(owner,'_average_children_'+kind)(root,db)
    assert success
    root.data=data
    if legacy:
        refusal(lambda:db.save_old_library(str(tmp_path/'derived.txt')))
    else:
        refusal(lambda:db.save_entry(StringIO(),root))
    # Existing unresolved-only text output remains available.
    root.item.electronic_state=[];root.item.vibrational_level=[]
    if legacy and kind == 'solute':
        with pytest.raises(NotImplementedError):
            db.save_old_library(str(tmp_path/'ordinary.txt'))
    elif legacy:
        db.save_old_library(str(tmp_path/'ordinary.txt'))
    else:
        db.save_entry(StringIO(),root)


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('input_kind', ['map', 'reactions'])
def test_regularize_explicit_inputs_named_refusal(state, input_kind):
    family, roots, rxn = training_fixture(state)
    kwargs = {'template_rxn_map': {roots[0].label: [rxn]}} if input_kind == 'map' else {'rxns': [rxn]}
    refusal(lambda: family.regularize(**kwargs), 'ResolvedStateTrainingError')


@pytest.mark.parametrize('state', STATES)
@pytest.mark.parametrize('kind', ['thermo', 'solute'])
@pytest.mark.parametrize('logic', ['OR{Pattern}', 'OR{Bridge}'])
@pytest.mark.parametrize('writer', ['save', 'save_entry', 'save_old_library'])
def test_logical_text_save_preserves_provenance(kind, state, logic, writer, tmp_path):
    from io import StringIO
    from rmgpy.data.base import make_logic_node
    db, root, child = tree(kind, child_data=datum(kind))
    pattern = Entry(label='Pattern', item=root.item.copy(deep=True))
    pattern.item.electronic_state = [state[0]] if state[0] else []
    pattern.item.vibrational_level = [state[1]] if state[1] >= 0 else []
    pattern.parent = root
    root.children.append(pattern)
    db.entries[pattern.label] = pattern
    if logic == 'OR{Bridge}':
        db.entries['Bridge'] = Entry(label='Bridge', item=make_logic_node('AND{Pattern}'))
    root.item = make_logic_node(logic)
    owner = ThermoDatabase() if kind == 'thermo' else SolvationDatabase()
    success, root.data = getattr(owner, '_average_children_' + kind)(root, db)
    assert success
    refusal(lambda: group_lookup(kind, db, molecule(state)), name='ExcitedSpeciesThermoError')
    def save():
        target = StringIO() if writer == 'save_entry' else str(tmp_path / 'derived.py')
        getattr(db, writer)(target, root) if writer == 'save_entry' else getattr(db, writer)(target)
    refusal(save)
    # A logical domain confined to unresolved species can still be written.
    pattern.item.electronic_state = []; pattern.item.vibrational_level = []
    if writer == 'save_old_library' and kind == 'solute':
        with pytest.raises(NotImplementedError):
            save()
    else:
        save()
