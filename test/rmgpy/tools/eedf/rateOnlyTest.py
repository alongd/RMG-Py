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

"""Rate-only chemistry must not become a second electron-energy operator."""

import copy
import pickle
import os

import numpy as np
import pytest

from rmgpy.kinetics import EEDFChannel
from rmgpy.tools.eedf.schema import SpecError
from rmgpy.tools.eedf.rate_only import (
    energy_schema, validate_energy_schema, validate_rate_only_marker,
    read_dissociation_source, SOURCE_SHA256, validate_spec_sources,
    validate_artifact,
)
from rmgpy.tools.eedf.moments import _distribution_moments


def channels():
    ordinary = dict(description='e + N2(X,v=0) -> e + N2(sing-sup), Excitation',
                    kind='excitation', threshold_eV=13., classification='B',
                    reaction=None, cross_section={'energy_eV': [13., 20.], 'sigma_m2': [0., 1e-20]})
    rate = dict(description='e + N2(X,v=0) -> e + N(4S) + N(4S), Excitation',
                kind='excitation', threshold_eV=10., classification='A', role='rate_only',
                source_sha256='1' * 64,
                reaction={'library': 'PlasmaAtomicNitrogen', 'index': 18,
                          'repr': "EEDFChannel(process='e + N2(X,v=0) -> e + N(4S) + N(4S), Excitation', collision_set='test', role='rate_only')"},
                cross_section={'energy_eV': [10., 20.], 'sigma_m2': [0., 1e-20]})
    return [ordinary, rate]


def test_rates_do_not_contribute_electron_power():
    row = dict(energy_eV=np.arange(.5, 30.), energy_edges_eV=np.arange(31.),
               f0=np.exp(-np.arange(.5, 30.) / 3.), target_fractions=[1., 1.],
               product_fractions=[0., 0.], gas_temperature_K=300.)
    result = _distribution_moments(row, channels())
    assert result['k_ine'][1] > 0.
    assert result['channel_power'][0] > 0.
    assert result['channel_power'][1] == 0.


def test_energy_schema_refuses_missing_chemical_metadata():
    data = channels()
    schema = energy_schema(data)
    validate_energy_schema(data, schema)
    schema['chemical_rates'].clear()
    with pytest.raises(SpecError, match='energy schema'):
        validate_energy_schema(data, schema)


def test_energy_schema_requires_complete_nitrogen_coverage():
    data = channels()
    data.insert(0, dict(data[0], description='e + N2(X,v=0) -> e + N2(A3Su+), Excitation'))
    schema = energy_schema(data)
    schema['coverage']['channels'].pop()
    with pytest.raises(SpecError, match='energy schema'):
        validate_energy_schema(data, schema)


def test_energy_schema_requires_sing_sup_and_unique_owners():
    data = channels()
    with pytest.raises(SpecError, match='sing-sup'):
        energy_schema(data[1:])
    with pytest.raises(SpecError, match='duplicate'):
        energy_schema(data + [copy.deepcopy(data[1])])


def test_rate_sources_do_not_change_energy_operator_fingerprint():
    data = channels()
    first = energy_schema(data)
    data[1]['source_sha256'] = '2' * 64
    data[1]['cross_section']['sigma_m2'][1] *= 2.
    second = energy_schema(data)
    assert first['energy_operator_sha256'] == second['energy_operator_sha256']
    assert first['chemical_rates'] != second['chemical_rates']


def test_rate_only_marker_round_trip_and_identity():
    marker = EEDFChannel(process='SS', collision_set='test', role='rate_only')
    assert marker.role == 'rate_only'
    assert pickle.loads(pickle.dumps(marker)).is_identical_to(marker)
    assert eval(repr(marker), {'EEDFChannel': EEDFChannel}).is_identical_to(marker)
    assert not EEDFChannel(process='SS', collision_set='test').is_identical_to(marker)
    with pytest.raises(ValueError, match='rate.only'):
        marker.side = 'sup'


@pytest.mark.parametrize('side,product', [('sup', None), ('ine', 'Ar(1s5)')])
def test_rate_only_rejects_superelastic_or_cascade(side, product):
    with pytest.raises(ValueError, match='rate.only'):
        validate_rate_only_marker('rate_only', side, product)


@pytest.fixture
def physical_source():
    source = os.environ.get('EEDF_DISSOCIATION_SOURCE')
    if not source:
        pytest.skip('pinned dissociation source required')
    return source


def test_energy_dependent_branch_split(physical_source):
    sections = read_dissociation_source(physical_source)
    assert [s['threshold_eV'] for s in sections] == [10., 12.]
    assert [len(s['cross_section']['energy_eV']) for s in sections] == [7, 24]
    values = [np.interp([14., 18., 20.], s['cross_section']['energy_eV'],
                        s['cross_section']['sigma_m2']) for s in sections]
    assert values[0] / (values[0] + values[1]) == pytest.approx(
        [.578815, .2229309872605626, .007118262386286167], rel=1e-12)


def test_refuse_source_in_collision_inputs():
    with pytest.raises(SpecError, match='normal collision'):
        validate_spec_sources({'input_files': {'extra': {'sha256': SOURCE_SHA256}}}, [])


def test_refuse_altered_rate_source(tmp_path, physical_source):
    from pathlib import Path
    from rmgpy.tools.eedf.schema import FingerprintMismatch
    path = tmp_path / 'extra.txt'
    path.write_bytes(Path(physical_source).read_bytes() + b'\n')
    with pytest.raises(FingerprintMismatch, match='source'):
        read_dissociation_source(path)


def test_refuse_chemical_columns_with_stripped_role_and_schema():
    data = channels()
    del data[-1]['role']
    with pytest.raises(SpecError, match='energy schema'):
        validate_artifact(data, {})


def test_full_spec_source_guards(physical_source):
    data = channels()
    first = data[-1]
    data = data[:1]
    for index, section in enumerate(read_dissociation_source(physical_source), 18):
        entry = dict(first, **section, source_sha256=SOURCE_SHA256)
        entry['reaction'] = dict(first['reaction'], index=index,
                                 repr=f"EEDFChannel(process={section['description']!r}, collision_set='test', role='rate_only')")
        data.append(entry)
    spec = {'input_files': {}, 'rate_only': {'sources': {SOURCE_SHA256: physical_source},
                                            'energy_schema': energy_schema(data)}}
    validate_spec_sources(spec, data)
    data[-1]['threshold_eV'] += 1.
    spec['rate_only']['energy_schema'] = energy_schema(data)
    with pytest.raises(ValueError, match='source'):
        validate_spec_sources(spec, data)
