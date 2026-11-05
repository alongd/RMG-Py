"""Argon population realization and pinned-solver applicability."""
import copy
import pytest
from rmgpy.tools.eedf.argon_population import (
    AXIS, NODES, CONVENTION, populations, validate_population_axis, check_zero_target_nodes)
from rmgpy.tools.eedf.schema import SpecError


@pytest.mark.parametrize('fraction', NODES + [3.7e-8, 2.3e-6])
def test_statistical_populations_preserve_requested_ratio(fraction):
    result = populations(fraction)
    ground = result.pop('Ar(1S0)')
    assert sum(result.values()) / ground == pytest.approx(fraction, rel=5.e-16, abs=1.e-25)
    assert ground + sum(result.values()) == pytest.approx(1., abs=1.e-15)
    if fraction:
        assert [v / sum(result.values()) for v in result.values()] == pytest.approx([5/12,3/12,1/12,3/12])


@pytest.mark.parametrize('mutation', ['declaration', 'axis', 'old_bound', 'weights'])
def test_population_axis_refuses_missing_or_inconsistent_convention(mutation):
    spec = {'axes': {AXIS: list(NODES)}, 'argon_4s_population': copy.deepcopy(CONVENTION)}
    validate_population_axis(spec)
    if mutation == 'declaration':
        del spec['argon_4s_population']
    elif mutation == 'axis':
        del spec['axes'][AXIS]
    elif mutation == 'old_bound':
        spec['excited_fraction_max'] = 1.e-5
    else:
        spec['argon_4s_population']['weights'] = [1,1,1,1]
    with pytest.raises(SpecError):
        validate_population_axis(spec)


def test_production_node_check_covers_all_nodes_and_refuses_zero_target(monkeypatch):
    import rmgpy.tools.eedf.integrity as integrity
    seen = []
    def resolve(spec, point):
        seen.append(point[AXIS])
        return {}, {}, {'Ar': 1.}, populations(point[AXIS])
    monkeypatch.setattr(integrity, '_resolved_properties', resolve)
    channel = {'description': 'e + Ar(3P2) <-> e + Ar(3P1), Excitation',
               'kind': 'excitation', 'target_fraction': 'state:Ar(3P2)',
               'product_fraction': 'state:Ar(3P1)'}
    spec = {'axes': {'u': [1., 2.], AXIS: NODES}}
    assert check_zero_target_nodes(spec, [channel]) == 5
    assert seen == NODES
    def invalid(spec, point):
        return {}, {}, {'Ar': 1.}, {'Ar(3P2)': 0., 'Ar(3P1)': 1.e-5}
    monkeypatch.setattr(integrity, '_resolved_properties', invalid)
    with pytest.raises(SpecError, match='zero-population target with populated product'):
        check_zero_target_nodes(spec, [channel])


def test_scan_agreement_uses_physical_quantities_not_relative_error_of_residual():
    import numpy as np
    from rmgpy.tools.eedf.validation import _rows_agree
    a = {key: np.array([1.e-15]) for key in
         ('k_ine', 'k_sup', 'channel_power', 'attachment_energy_eV', 'f0')}
    a.update(swarm={'mobility_N': 1.e24}, power_groups={
        'field': 1.e-16, 'balance': 1.e-28, 'relative_balance_percent': 1.e-9},
        converged=True)
    b = copy.deepcopy(a)
    b['power_groups']['relative_balance_percent'] += 2.3e-11
    tolerance = {'rtol': 1.e-6, 'atol': 1.e-24}
    assert _rows_agree(a, b, tolerance)
    b['swarm']['mobility_N'] *= 1.00001
    assert not _rows_agree(a, b, tolerance)
    b = copy.deepcopy(a)
    b['power_groups']['field'] *= 1.00001
    assert not _rows_agree(a, b, tolerance)


def test_population_refinement_preserves_frozen_coordinate_transition():
    spec = {'axes': {AXIS: sorted(NODES + [1.e-7 / 3, 2.e-7 / 3, 4.e-7])},
            'argon_4s_population': copy.deepcopy(CONVENTION)}
    validate_population_axis(spec)
    from rmgpy.tools.eedf.argon_population import coordinate
    assert coordinate([1.e-7 / 3, 1.e-7]).tolist() == pytest.approx([1/3, 1.])
    spec['axes'][AXIS].remove(1.e-7)
    with pytest.raises(SpecError, match='five base nodes'):
        validate_population_axis(spec)
