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

"""Shared validation of sourced collision ownership and cascade yields.

The generated yield document owns the numbers. Consumers validate a byte
snapshot and never renormalize it or turn an unresolved family into chemistry.
"""

import ast
import copy
import hashlib
import math
from pathlib import Path
import re

import yaml

from rmgpy.tools.eedf.schema import FingerprintMismatch, SpecError


YIELD_SUM_ATOL = 1e-12
ARGON_PRODUCTS = ('Ar(1s5)', 'Ar(1s4)', 'Ar(1s3)', 'Ar(1s2)')
FAMILY_SCHEMA = 'eedf-branch-family/1'
ARGON_STATES = {'Ar(1s5)': ('3P2', 3), 'Ar(1s4)': ('3P1', 3),
                'Ar(1s3)': ('3P0', 3), 'Ar(1s2)': ('1P1', 1)}


class _UniqueLoader(yaml.SafeLoader):
    """Refuse duplicate source keys instead of accepting the last value."""

    def construct_mapping(self, node, deep=False):
        result = {}
        for key_node, value_node in node.value:
            key = self.construct_object(key_node, deep=deep)
            if type(key) not in (str, int) or key in result:
                raise SpecError('duplicate or invalid cascade source key: ' + repr(key))
            result[key] = self.construct_object(value_node, deep=deep)
        return result


def validate_branch_yield(value):
    """Accept a finite plain numeric fraction, including an explicit zero."""
    if type(value) not in (int, float) or not 0 <= value <= 1 or not math.isfinite(value):
        raise SpecError('branch yield must be finite and between zero and one')
    return value


def validate_branch_metadata(product, value, digest, side):
    """Validate an atomic marker triple without imposing a family sum."""
    if product is None and value is None and digest is None:
        return
    if not isinstance(product, str) or product not in ARGON_PRODUCTS:
        raise SpecError('branch product must name a sourced argon 1s state')
    validate_branch_yield(value)
    if not isinstance(digest, str) or not re.fullmatch('[0-9a-f]{64}', digest):
        raise SpecError('branch yield_table_sha256 must be a source digest')
    if side != 'ine':
        raise SpecError('cascade branch cannot use the superelastic side')


def validate_yield_vector(values):
    """Require all four sourced products and unity within an absolute 1e-12.

    This is the sole engine normalization rule. Return a detached copy with
    unchanged values; absent, null and zero yields have different meanings.
    """
    if not isinstance(values, dict) or set(values) != set(ARGON_PRODUCTS):
        raise SpecError('cascade yield products must contain all four argon 1s states')
    for value in values.values():
        validate_branch_yield(value)
    if abs(math.fsum(values.values()) - 1.) > YIELD_SUM_ATOL:
        raise SpecError('cascade yields must sum to one within 1e-12')
    return dict(values)


def _read_source(path, expected):
    """Hash and return the same bytes, even if the source path is replaced."""
    try:
        data = Path(path).read_bytes()
    except OSError as error:
        raise SpecError('cannot read cascade source: ' + str(error)) from error
    digest = hashlib.sha256(data).hexdigest()
    if expected is not None and expected != digest:
        raise FingerprintMismatch('cascade source hash: ' + str(path))
    return data, digest


def _validate_source_row(row):
    """Check source completeness without inferring missing transition data."""
    if not isinstance(row, dict):
        raise SpecError('invalid cascade source channel')
    if (type(row.get('id')) is not int or row['id'] < 1
            or not isinstance(row.get('process'), str) or not row['process'].strip()):
        raise SpecError('invalid cascade collision identity')
    threshold = row.get('threshold_eV')
    if type(threshold) not in (int, float) or not math.isfinite(threshold) or threshold < 0:
        raise SpecError('invalid cascade threshold')
    status = row.get('status')
    if status == 'unresolved':
        if 'product_yields' not in row or row['product_yields'] is not None:
            raise SpecError('unresolved cascade must have null yields')
    elif status == 'complete_for_E1':
        validate_yield_vector(row.get('product_yields'))
        transitions = row.get('transition_ids')
        if (row.get('missing_transitions') != [] or not isinstance(transitions, list)
                or not transitions or any(not isinstance(t, str) or not t for t in transitions)
                or len(set(transitions)) != len(transitions)):
            raise SpecError('complete cascade has missing or invalid transition provenance')
    else:
        raise SpecError('unsupported cascade source status')


def load_cascade_yields(path, radiative_source, *, yield_table_sha256=None, source_sha256=None):
    """Read a generated E1 yield document and verify its radiative source.

    Optional pins constrain the external bytes. Both files are read once;
    the returned snapshot contains their digests, never machine-local paths.
    Unresolved records may be carried as evidence, but cannot be requested
    as usable families through :func:`require_cascade_yields`.
    """
    data, digest = _read_source(path, yield_table_sha256)
    _, raw_digest = _read_source(radiative_source, source_sha256)
    try:
        document = yaml.load(data, Loader=_UniqueLoader)
    except (yaml.YAMLError, UnicodeError) as error:
        raise SpecError('malformed cascade yield YAML') from error
    if (not isinstance(document, dict)
            or set(document) != {'schema_version', 'radiative_model', 'source_sha256', 'channels'}
            or type(document['schema_version']) is not int or document['schema_version'] != 1
            or document['radiative_model'] != 'E1'):
        raise SpecError('cascade yield document schema')
    source_hash = document['source_sha256']
    if not isinstance(source_hash, str) or not re.fullmatch('[0-9a-f]{64}', source_hash):
        raise SpecError('invalid cascade radiative source hash')
    if source_hash != raw_digest:
        raise FingerprintMismatch('cascade radiative source')
    rows = document['channels']
    if not isinstance(rows, list) or not rows:
        raise SpecError('empty cascade source channels')
    ids, processes = set(), set()
    for row in rows:
        _validate_source_row(row)
        if row['id'] in ids or row['process'] in processes:
            raise SpecError('duplicate cascade collision identity')
        ids.add(row['id'])
        processes.add(row['process'])
    snapshot = copy.deepcopy(document)
    snapshot['yield_table_sha256'] = digest
    return snapshot


def require_cascade_yields(snapshot, channel_id, *, process=None, threshold_eV=None):
    """Return one complete sourced row, cross-checking its physical collision."""
    if type(channel_id) is not int or channel_id < 1:
        raise SpecError('invalid cascade channel id')
    rows = [row for row in snapshot['channels'] if row['id'] == channel_id]
    if len(rows) != 1:
        raise SpecError('missing or duplicate cascade source family: ' + str(channel_id))
    row = rows[0]
    _validate_source_row(row)
    if row['status'] != 'complete_for_E1':
        raise SpecError('unresolved cascade yields: ' + str(channel_id))
    if ((process is not None and row['process'] != process)
            or (threshold_eV is not None and row['threshold_eV'] != threshold_eV)):
        raise SpecError('cascade collision identity differs from sourced row: ' + str(channel_id))
    return copy.deepcopy(row)


def is_branch_family(identity):
    """Recognize the ownership union without accepting unknown schema versions."""
    return isinstance(identity, dict) and 'schema' in identity


def marker_fields(representation):
    """Read marker literals without executing code from an artifact."""
    try:
        call = ast.parse(representation, mode='eval').body
        if (not isinstance(call, ast.Call) or not isinstance(call.func, ast.Name)
                or call.func.id != 'EEDFChannel' or call.args
                or any(k.arg is None for k in call.keywords)
                or len({k.arg for k in call.keywords}) != len(call.keywords)):
            raise ValueError('not an EEDFChannel keyword constructor')
        fields = {item.arg: ast.literal_eval(item.value) for item in call.keywords}
    except (TypeError, SyntaxError, ValueError, AttributeError) as error:
        raise SpecError('unreadable branch marker metadata') from error
    required = {'process', 'collision_set', 'side', 'branch_product',
                'branch_yield', 'yield_table_sha256'}
    if not required <= set(fields) or set(fields) - required - {'Tmin', 'Tmax', 'Pmin', 'Pmax', 'comment'}:
        raise SpecError('branch marker fields')
    validate_branch_metadata(fields['branch_product'], fields['branch_yield'],
                             fields['yield_table_sha256'], fields['side'])
    return fields


def validate_family(family, *, process=None, threshold_eV=None):
    """Validate the persisted ownership envelope, returning a detached snapshot.

    The source digests and complete vector are frozen into the table identity.
    Runtime does not discover a newer database file or reinterpret its yields.
    """
    required = {'schema', 'channel_id', 'process', 'threshold_eV', 'collision_set',
                'radiative_model', 'status', 'yield_table_sha256', 'source_sha256',
                'sum_atol', 'product_yields', 'transition_ids', 'branches'}
    if not isinstance(family, dict) or set(family) != required or family['schema'] != FAMILY_SCHEMA:
        raise SpecError('branch family schema/fields')
    if (family['radiative_model'] != 'E1' or family['status'] != 'complete_for_E1'
            or type(family['sum_atol']) not in (int, float) or family['sum_atol'] != YIELD_SUM_ATOL
            or not isinstance(family['collision_set'], str) or not family['collision_set'].strip()):
        raise SpecError('branch family source policy')
    _validate_source_row(dict(id=family['channel_id'], process=family['process'],
                              threshold_eV=family['threshold_eV'], status=family['status'],
                              product_yields=family['product_yields'], missing_transitions=[],
                              transition_ids=family['transition_ids']))
    for name in ('yield_table_sha256', 'source_sha256'):
        if not isinstance(family[name], str) or not re.fullmatch('[0-9a-f]{64}', family[name]):
            raise SpecError('branch family source digest')
    if ((process is not None and family['process'] != process)
            or (threshold_eV is not None and family['threshold_eV'] != threshold_eV)):
        raise SpecError('branch family collision identity')
    branches = family['branches']
    if not isinstance(branches, list) or len(branches) != len(ARGON_PRODUCTS):
        raise SpecError('branch family must contain all four branches, including zero yields')
    products, entries = [], set()
    for owner in branches:
        if not isinstance(owner, dict) or set(owner) != {'product', 'library', 'index', 'repr'}:
            raise SpecError('branch owner fields')
        if (not isinstance(owner['product'], str) or owner['product'] not in ARGON_PRODUCTS
                or not isinstance(owner['library'], str) or not owner['library'].strip()
                or type(owner['index']) is not int or owner['index'] < 0):
            raise SpecError('branch owner identity')
        key = (owner['library'], owner['index'])
        if key in entries:
            raise SpecError('duplicate branch source entry')
        entries.add(key)
        products.append(owner['product'])
        fields = marker_fields(owner['repr'])
        if (fields['process'] != family['process'] or fields['collision_set'] != family['collision_set']
                or fields['branch_product'] != owner['product']
                or fields['branch_yield'] != family['product_yields'][owner['product']]
                or fields['yield_table_sha256'] != family['yield_table_sha256']):
            raise FingerprintMismatch('branch marker differs from sourced family')
    if products != list(ARGON_PRODUCTS):
        raise SpecError('branch family products must be unique and in source order')
    return copy.deepcopy(family)


def iter_reaction_owners(identity):
    """Enumerate source identities while keeping one envelope per collision."""
    if identity is None:
        return iter(())
    if is_branch_family(identity):
        family = validate_family(identity)
        return iter({key: owner[key] for key in ('library', 'index', 'repr')}
                    for owner in family['branches'])
    if not isinstance(identity, dict) or set(identity) != {'library', 'index', 'repr'}:
        raise SpecError('mapped reaction identity')
    return iter((copy.deepcopy(identity),))


def _argon_identity(species):
    """Require a neutral single atom with a consistent resolved state."""
    molecules = species.molecule
    if not molecules:
        raise SpecError('cascade argon species requires molecular identity')
    identities = set()
    for molecule in molecules:
        if (len(molecule.atoms) != 1 or molecule.atoms[0].element.symbol != 'Ar'
                or molecule.get_net_charge() != 0 or molecule.vibrational_level != -1):
            raise SpecError('cascade requires neutral monatomic argon')
        identities.add((molecule.electronic_state, molecule.multiplicity))
    if len(identities) != 1:
        raise SpecError('inconsistent cascade argon molecular states')
    return identities.pop()


def reaction_product_state(reaction, process=None):
    """Check collision stoichiometry and return the actual Paschen product.

    Display aliases (notably Ars for 1s5) carry no identity authority.
    """
    if reaction.reversible:
        raise SpecError('cascade branch must be irreversible')
    heavy = []
    for side in (reaction.reactants, reaction.products):
        if len(side) != 2 or sum(species.is_electron() for species in side) != 1:
            raise SpecError('cascade branch requires one electron and one argon on each side')
        heavy.append(next(species for species in side if not species.is_electron()))
    if process is None:
        process = getattr(reaction.kinetics, 'process', None)
    source = re.match(r'^e \+ Ar\((1S0|3P2|3P1|3P0|1P1)\) -> ', process or '')
    if source is None:
        raise SpecError('cascade process requires a resolved argon source state')
    actual_source = _argon_identity(heavy[0])
    if source[1] == '1S0':
        if actual_source not in (('', 1), ('1S0', 1)):
            raise SpecError('cascade source reactant must be ground-state argon')
    elif actual_source != (source[1], 1 if source[1] == '1P1' else 3):
        raise SpecError('cascade source differs from molecular reactant state')
    state = _argon_identity(heavy[1])
    for product, identity in ARGON_STATES.items():
        if state == identity:
            return product
    raise SpecError('cascade product is not a sourced argon 1s state')


def make_branch_family(snapshot, channel_id, process, threshold_eV, collision_set, owners):
    """Cross-check library reactions against the source and freeze the family.

    Owners are (library name, entry index, reaction, kinetics marker) tuples;
    kinetics live in entry.data before a library reaction is instantiated.
    """
    row = require_cascade_yields(snapshot, channel_id, process=process, threshold_eV=threshold_eV)
    branches = []
    for library, index, reaction, marker in owners:
        product = reaction_product_state(reaction, process)
        if (not getattr(marker, 'uses_eedf', False)
                or marker.branch_product != product):
            raise SpecError('branch marker product differs from molecular product')
        branches.append(dict(product=product, library=library, index=index, repr=repr(marker)))
    branches.sort(key=lambda branch: ARGON_PRODUCTS.index(branch['product']))
    return validate_family(dict(
        schema=FAMILY_SCHEMA, channel_id=channel_id, process=process, threshold_eV=threshold_eV,
        collision_set=collision_set, radiative_model=snapshot['radiative_model'], status=row['status'],
        yield_table_sha256=snapshot['yield_table_sha256'], source_sha256=snapshot['source_sha256'],
        sum_atol=YIELD_SUM_ATOL, product_yields=row['product_yields'],
        transition_ids=row['transition_ids'], branches=branches))


def validate_bound_branch(family, reaction):
    """Match a runtime reaction to a previously validated family snapshot."""
    product = reaction_product_state(reaction, family['process'])
    marker = reaction.kinetics
    entry = getattr(reaction, 'entry', None)
    actual = dict(product=product, library=getattr(reaction, 'library', None),
                  index=getattr(entry, 'index', None), repr=repr(marker))
    expected = next(owner for owner in family['branches'] if owner['product'] == product)
    if actual != expected:
        raise FingerprintMismatch('cascade branch reaction identity differs from sourced family')
    return product


def is_declared_branch_owner(reaction, library):
    """Allow admission only for a marker explicitly named by a complete family.

    The declaration travels in entry.long_desc through library serialization.
    Table binding still checks the sourced digest and full frozen family.
    """
    from rmgpy.kinetics import EEDFChannel
    marker = reaction.kinetics
    entry = getattr(reaction, 'entry', None)
    if (not isinstance(marker, EEDFChannel) or marker.branch_product is None
            or entry is None or getattr(reaction, 'library', None) != library.label
            or library.entries.get(entry.index) is not entry
            or repr(marker) != repr(entry.data) or reaction.reversible):
        return False
    declarations = re.findall(r'^EEDF branch family = (.+)$', entry.long_desc, re.M)
    if len(declarations) != 1:
        return False
    try:
        indices = ast.literal_eval(declarations[0])
        if (not isinstance(indices, tuple) or len(indices) != 4
                or any(type(index) is not int for index in indices)
                or len(set(indices)) != 4 or entry.index not in indices):
            return False
        yields = {}
        for index in indices:
            sibling = library.entries[index]
            kin = sibling.data
            if (not isinstance(kin, EEDFChannel) or sibling.item.reversible
                    or kin.process != marker.process
                    or kin.collision_set != marker.collision_set
                    or kin.yield_table_sha256 != marker.yield_table_sha256
                    or kin.side != 'ine'
                    or re.findall(r'^EEDF branch family = (.+)$', sibling.long_desc, re.M) != declarations
                    or kin.branch_product in yields
                    or reaction_product_state(sibling.item, kin.process) != kin.branch_product):
                return False
            yields[kin.branch_product] = kin.branch_yield
        validate_yield_vector(yields)
        validate_branch_metadata(marker.branch_product, marker.branch_yield,
                                marker.yield_table_sha256, marker.side)
        return reaction_product_state(reaction, marker.process) == marker.branch_product
    except (ValueError, SyntaxError, KeyError, TypeError):
        return False


def check_branch_promotion(core, edge, proposed):
    """Preflight a complete promotion before changing any model membership.

    Initial seed admission is incremental. Once any family member is in the
    edge, admission to core must move the entire four-product family at once.
    Reactor binding separately refuses an incomplete initial seed family.
    """
    def key(reaction):
        marker = reaction.kinetics
        if getattr(marker, 'branch_product', None) is None:
            return None
        if not isinstance(marker.process, str) or not isinstance(marker.collision_set, str):
            raise SpecError('cascade promotion requires a string collision identity')
        return (marker.process, marker.collision_set, marker.yield_table_sha256)

    touched = {key(reaction) for reaction in proposed} - {None}
    if not touched:
        return
    edge_keys = {key(reaction) for reaction in edge} - {None}
    for family_key in touched & edge_keys:
        members = []
        for reaction in list(core) + list(edge) + list(proposed):
            if key(reaction) == family_key and reaction not in members:
                members.append(reaction)
        products = [reaction.kinetics.branch_product for reaction in members]
        if len(products) != len(ARGON_PRODUCTS) or set(products) != set(ARGON_PRODUCTS):
            raise SpecError('cascade promotion requires a complete four-product family')
        if any(reaction in edge and reaction not in proposed for reaction in members):
            raise SpecError('split cascade family promotion between core and edge is unsupported')
