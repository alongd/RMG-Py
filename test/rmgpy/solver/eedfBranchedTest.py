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

"""Literal sourced product rates through the public artifact/provider boundary."""

import copy
import importlib.util
import json
import math
from pathlib import Path

import pytest

from rmgpy.solver.eedf_provider import EEDFProvider
from rmgpy.tools.eedf.schema import _file_hash
from plasmaEEDFCoreTest import artifact, species_and_reaction


_spec = importlib.util.spec_from_file_location(
    'cascade_test_source', Path(__file__).parents[1] / 'tools/eedf/branchYieldTest.py')
_source = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_source)
source_files = _source.source_files
branch_fixture = _source.branch_fixture
branch_reactions = _source.branch_reactions
family_from_reactions = _source.family_from_reactions
ID6 = _source.ID6


@pytest.fixture
def family_table(tmp_path, branch_fixture):
    _, reactions, family = branch_fixture
    path = artifact(tmp_path, reactions[0].kinetics, composition=False,
                    include_ar4s_envelope=False, family=family)
    model = json.loads((path / 'manifest.json').read_text())['row_inputs']
    return path, model, {i: (0, 'ine', product) for i, product in enumerate(ID6)}


def provider_for(table, mapping=None):
    path, model, original = table
    return EEDFProvider(path, model, artifact_sha256=_file_hash(path / 'table.h5'),
                         reaction_map=original if mapping is None else mapping, branch='branch_0')


def test_id6_literal_provider_rates(family_table):
    provider = provider_for(family_table)
    row = provider.row(0., {})
    rates = [provider.reaction_rate(i, row) for i in range(4)]
    assert rates == pytest.approx([
        7.4204946996466445e-16, 2.120141342756184e-16,
        3.8476639183352967e-17, 7.45975657636435e-18], rel=2e-14, abs=0.)
    assert math.fsum(rates) == pytest.approx(1e-15, rel=2e-14, abs=0.)
    for i, weight in enumerate(ID6.values()):
        assert provider.reaction_yield(i) == weight
        assert rates[i] == float(row.k_ine[0]) * weight


@pytest.mark.parametrize('problem', ['missing', 'extra', 'ordinary', 'sup', 'unknown', 'empty'])
def test_provider_requires_complete_unique_family(family_table, problem):
    mapping = dict(family_table[2])
    if problem == 'missing':
        mapping.pop(3)
    elif problem == 'extra':
        mapping[4] = mapping[0]
    elif problem == 'ordinary':
        mapping[0] = (0, 'ine')
    elif problem == 'sup':
        mapping[0] = (0, 'sup', 'Ar(1s5)')
    elif problem == 'unknown':
        mapping[0] = (0, 'ine', 'Ar(4p)')
    else:
        mapping.clear()
    with pytest.raises(ValueError, match='cascade|double-count'):
        provider_for(family_table, mapping)


def test_binding_owns_family_snapshot(family_table):
    provider = provider_for(family_table)
    before = provider.reaction_rate(0, provider.row(.5, {}))
    family_table[1]['reactions'][0]['product_yields']['Ar(1s5)'] = 0.
    family_table[2][0] = (0, 'ine', 'Ar(1s2)')
    assert provider.reaction_rate(0, provider.row(.5, {})) == before
    with pytest.raises(TypeError):
        provider.runtime_metadata['channel_map'][0]['reaction']['product_yields']['Ar(1s5)'] = 0.


def test_legacy_single_owner_bit_identical_and_duplicate_refused(tmp_path):
    _, reaction = species_and_reaction()
    path = artifact(tmp_path, reaction.kinetics, composition=False, include_ar4s_envelope=False)
    model = json.loads((path / 'manifest.json').read_text())['row_inputs']
    provider = provider_for((path, model, {0: (0, 'ine')}))
    for u in [0., .5, 1., 2.3, 3.]:
        row = provider.row(u, {})
        assert provider.reaction_rate(0, row).hex() == float(row.k_ine[0]).hex()
        assert provider.reaction_yield(0) == 1.
    with pytest.raises(ValueError, match='double-count one tabulated rate'):
        provider_for((path, model, {0: (0, 'ine'), 1: (0, 'ine')}))


def test_explicit_zero_provider_rate(tmp_path, branch_fixture):
    snapshot, _, _ = branch_fixture
    snapshot['channels'][0]['product_yields'] = dict(zip(ID6, [1., 0., 0., 0.]))
    reactions = branch_reactions(snapshot)
    family = family_from_reactions(snapshot, reactions)
    path = artifact(tmp_path, reactions[0].kinetics, composition=False,
                    include_ar4s_envelope=False, family=family)
    model = json.loads((path / 'manifest.json').read_text())['row_inputs']
    provider = provider_for((path, model, {i: (0, 'ine', product) for i, product in enumerate(ID6)}))
    row = provider.row(1., {})
    assert [provider.reaction_rate(i, row) for i in [1, 2, 3]] == [0., 0., 0.]
