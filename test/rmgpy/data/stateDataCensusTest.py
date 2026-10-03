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

"""Static coverage of every database function and explicit data-route verdicts."""

import json
from pathlib import Path
import runpy

import pytest

ROOT = Path(__file__).resolve().parents[3]
CENSUS = runpy.run_path(str(ROOT / 'scripts/state_data_census.py'))
MANIFEST = ROOT / 'test/rmgpy/data/stateDataCensus.json'


def test_generated_census_is_classified():
    sites = CENSUS['collect_sites'](ROOT)
    classifications = json.loads(MANIFEST.read_text())['sites']
    CENSUS['check_classifications'](sites, classifications)
    for row in classifications.values():
        for evidence in row['evidence']:
            assert (ROOT / evidence).is_file(), evidence


@pytest.mark.parametrize('mutation', ['remove_site', 'remove_verdict', 'source_changed', 'new_site'])
def test_census_rejects_missing_or_stale_classification(mutation):
    sites = CENSUS['collect_sites'](ROOT)
    classifications = json.loads(MANIFEST.read_text())['sites']
    site = 'rmgpy/data/base.py:Database.descend_tree'
    if mutation == 'remove_site':
        del classifications[site]
    elif mutation == 'remove_verdict':
        del classifications[site]['verdict']
    elif mutation == 'source_changed':
        classifications[site]['source_sha256'] = 'stale'
    else:
        sites['rmgpy/data/new.py:opaque_lookup'] = sites[site]
    with pytest.raises(ValueError):
        CENSUS['check_classifications'](sites, classifications)


def test_ast_census_includes_aliases_wrappers_and_nested_callbacks(tmp_path):
    directory = tmp_path / 'rmgpy/data/kinetics'
    directory.mkdir(parents=True)
    (directory / 'opaque.py').write_text('''
from package import Entry as Record

def wrapper(db):
    return opaque(db)

async def opaque(db):
    def callback(obj):
        return getattr(obj, "data")
    return callback(Record()), lambda obj: obj.data

class Routes:
    def selection(self, db):
        return db.descend_tree()
''')
    sites = CENSUS['collect_sites'](tmp_path)
    assert {key.split(':')[1] for key in sites} == {
        'wrapper', 'opaque', 'opaque.callback', 'opaque.<lambda1>', 'Routes.selection'}


@pytest.mark.parametrize('directory', [
    'rmgpy/qm', 'rmgpy/ml', 'rmgpy/rmg', 'rmgpy/statmech',
    'rmgpy/pdep', 'rmgpy/thermo', 'rmgpy/kinetics', 'arkane/ess', 'arkane',
])
def test_new_supplier_is_discovered_and_requires_a_verdict(directory, tmp_path):
    source = tmp_path / directory / 'new_supplier.py'
    source.parent.mkdir(parents=True, exist_ok=True)
    source.write_text('def supply(species):\n    return species.thermo\n')
    sites = CENSUS['collect_sites'](tmp_path)
    assert str(source.relative_to(tmp_path)) + ':supply' in sites
    with pytest.raises(ValueError, match='unclassified new site'):
        CENSUS['check_classifications'](sites, {})


def test_cython_supplier_properties_callbacks_and_declarations(tmp_path):
    path = tmp_path / 'rmgpy/pdep/supplier.pyx'
    path.parent.mkdir(parents=True)
    path.write_text('cdef class Supplier:\n    property value:\n        def __get__(self):\n            return 1\n    cpdef double supply(self, double temperature):\n        return temperature\n    def wrapper(self):\n        def callback():\n            return 2\n        return callback(), lambda x: x\n')
    sites = CENSUS['collect_sites'](tmp_path)
    assert {key.split(':')[1] for key in sites} == {
        'Supplier.value.__get__', 'Supplier.supply', 'Supplier.wrapper',
        'Supplier.wrapper.callback', 'Supplier.wrapper.<lambda1>'}
    path.with_suffix('.pxd').write_text('cdef class Supplier:\n    cpdef double supply(self, double temperature)\n')
    changed = CENSUS['collect_sites'](tmp_path)
    assert all(sites[key]['source_sha256'] != changed[key]['source_sha256'] for key in sites)


def test_inventory_covers_independent_assignment_boundaries():
    sites = CENSUS['collect_sites'](ROOT)
    for key in [
        'rmgpy/species.py:Species.generate_energy_transfer_model',
        'rmgpy/qm/main.py:QMCalculator.get_thermo_data',
        'rmgpy/ml/estimator.py:MLEstimator.get_thermo_data',
        'rmgpy/rmg/model.py:CoreEdgeReactionModel.generate_thermo',
        'rmgpy/pdep/collision.pyx:SingleExponentialDown.get_alpha',
        'rmgpy/thermo/thermoengine.py:generate_thermo_data',
        'rmgpy/statmech/conformer.pyx:Conformer.get_heat_capacity',
        'arkane/statmech.py:StatMechJob.load',
    ]:
        assert key in sites, key
