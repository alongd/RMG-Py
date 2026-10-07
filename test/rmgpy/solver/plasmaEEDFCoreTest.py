#!/usr/bin/env python3

import copy
import json

import numpy as np
import pytest
from scipy.optimize import brentq

import rmgpy.constants as constants
import rmgpy.data.rmg as rmg_data_module
import rmgpy.rmg.input as rmg_input
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.data.base import Entry
from rmgpy.exceptions import PlasmaStateError
from rmgpy.kinetics import EEDFChannel, TwoTemperaturePlasma
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.solver.eedf_provider import (
    DEVELOPMENT_UNQUALIFIED_STATUS,
    DevelopmentWallBudgetExceeded,
    development_unqualified_table_route,
)
from rmgpy.rmg.main import RMG
from rmgpy.species import Species
from rmgpy.thermo import ThermoData
from rmgpy.tools.eedf.artifact import write_artifact
from rmgpy.tools.eedf.schema import EEDFError, content_hash, file_hash, interpolant_identity


TGAS = 298.15
PRESSURE = 666.61
LAMBDA = 0.02
MU0 = 1.535e-4
VOLUME = 1.e-4


@pytest.fixture(autouse=True)
def isolate_database():
    saved = rmg_data_module.database
    rmg_data_module.database = None
    try:
        yield
    finally:
        rmg_data_module.database = saved


def thermo(energy_eV=0.):
    return ThermoData(Tdata=([300., 1000.], 'K'), Cpdata=([20., 20.], 'J/(mol*K)'),
                      H298=(energy_eV * constants.e * constants.Na, 'J/mol'),
                      S298=(0., 'J/(mol*K)'))


def species_and_reaction():
    electron = Species(label='e-').from_adjacency_list('1 e u1 p0 c-1')
    ar = Species(label='Ar').from_adjacency_list('1 Ar u0 p4 c0')
    ion = Species(label='Ar+').from_adjacency_list('multiplicity 2\n1 Ar u1 p3 c+1')
    excited = Species(label='Ars').from_adjacency_list('1 Ar u2 p3 c0')
    ar.thermo, ion.thermo, excited.thermo = thermo(), thermo(15.76), thermo(11.55)
    marker = EEDFChannel(process='e + Ar -> e + e + Ar+', collision_set='fixture', side='ine')
    reaction = LibraryReaction(reactants=[electron, ar], products=[ion, electron, electron],
                               reversible=False, kinetics=marker, library='Toy',
                               entry=Entry(index=1))
    reaction.index = 1
    return [electron, ar, ion, excited], reaction


def artifact(tmp_path, marker, composition=True, cold_mean=False, accepted=True,
             reaction_owned=True, gases=None, pressure=PRESSURE,
             epsilon_k_eV=None, include_ar4s_envelope=True,
             variable_epsilon=False):
    gases = list(gases or ['Ar'])
    axes = {'u': [0., 1., 2., 3.]}
    if composition:
        axes['x(Ars)'] = [0., .2]
    shape = tuple(len(values) for values in axes.values())
    rows = np.empty(shape, dtype=object)
    for index in np.ndindex(shape):
        x_index = index[1] if composition else 0
        scale = (1., 2., 4., 8.)[index[0]] * (1. + x_index)
        mean_energy = ((1.5 * constants.kB * TGAS / constants.e + .01 * index[0])
                       if cold_mean else 1. + index[0] + x_index)
        mobility_n = 1.e24 * scale
        epsilon_k = ((1. + axes['u'][index[0]]) ** 2 if variable_epsilon
                     else (2. if epsilon_k_eV is None else epsilon_k_eV))
        diffusion_n = epsilon_k * mobility_n
        rows[index] = {
            'u': axes['u'][index[0]], 'EN_Td': np.exp(axes['u'][index[0]]),
            'swarm': {'mean_energy_eV': mean_energy,
                      'mobility_N': mobility_n, 'diffusion_N': diffusion_n},
            'k_ine': np.array([scale * 1.e-15, scale * 2.e-14]),
            'k_sup': np.array([.2 * scale * 1.e-15, 0.]),
            'channel_power': np.array([15.76 * .8 * scale * 1.e-15,
                                       2.e-16 * scale]),
            'power_groups': {'field': (1.e24 * scale) *
                                      (np.exp(axes['u'][index[0]]) * 1.e-21) ** 2,
                             'elastic_gain': 1.e-16 * scale,
                             'elastic_loss': -3.e-16 * scale},
            'attachment_energy_eV': np.array([0., 0.]),
            'target_fractions': np.array([1., 1.]),
            'product_fractions': np.array([1., 0.]),
            'rate_floors': np.array([1.e-40, 1.e-40]),
            'below_floor': np.zeros((2, 2), bool),
            'f0': np.array([1.]), 'energy_eV': np.array([.5]),
            'energy_edges_eV': np.array([0., 1.]),
            'gas_fractions': {'Ar': 1.}, 'state_populations': {'Ar': 1.},
            'composition': ({'x(Ars)': axes['x(Ars)'][x_index]} if composition else {}),
            'converged': True, 'iteration_count': -1, 'setup': 'fixture',
        }
    identity = {'library': 'Toy', 'index': 1, 'repr': repr(marker)}
    model = {
        'axes': axes, 'Tg_K': TGAS, 'P_Pa': pressure,
        'arm': {'gases': gases, 'additive': None, 'feed_partial_pressure_Pa': 0.},
        'cross_sections': {'Ar': 'sha'}, 'channel_map_sha256': 'map',
        'reactions': ([identity] if reaction_owned else []),
        'solver_options': {'growth': 'temporal'},
        'interpolant': interpolant_identity(), 'state_populations': {'Ar': 1.},
    }
    envelopes = {
        'Tg_K': {'min': TGAS, 'max': TGAS, 'reference': TGAS},
        'P_Pa': {'min': pressure, 'max': pressure, 'reference': pressure},
        'electron_density': {'min': 0., 'max': 1.e30, 'reference': 1.e16},
    }
    if include_ar4s_envelope:
        envelopes['Ar4s_total'] = {'min': 0., 'max': .15, 'reference': .1}
    manifest = {
        'accepted': accepted,
        'held_out_verdicts': [{'point': {'u': .5}, 'branch_id': 'branch_0',
                               'passed': accepted, 'checks': [{'passed': accepted}]}],
        'branch_certification': {'branch_0': 'fixture'}, 'schema_version': 1,
        'row_inputs': model, 'fingerprint': content_hash(model), 'axes': axes,
        'envelopes': envelopes,
        'tolerances': {'fingerprint_rtol': 1.e-9, 'G1': {'rtol': 1.e-6, 'atol': 1.e-24},
                       'u_condition': {'min_abs_denergy_du': 1.e-6,
                                       'min_scaled_denergy_du': 1.e-6}},
        'floors': {'rate_absolute': 1.e-40, 'absolute_power_share': 1.e-4},
        'repositories': {},
        'channel_map': [{'description': marker.process, 'kind': 'ionization',
                         'classification': ('A' if reaction_owned else 'B'),
                         'reaction': (identity if reaction_owned else None)},
                        {'description': 'e + Ar -> e + Ar, Elastic', 'kind': 'elastic',
                         'classification': 'B', 'reaction': None, 'mass_ratio': 1.e-5}],
        'energy_eV': [.5], 'energy_edges_eV': [0., 1.],
    }
    path = write_artifact(tmp_path, manifest, {'branch_0': rows}, [])
    return path


def build_reactor(tmp_path, composition=True, extra_reaction=None, empirical_laws=None,
                  reaction_index=1, entry_index=1, power_w=.5, cold_mean=False,
                  extra_species=None, electron_energies=None, accepted=True,
                  development=False, excited_label=None):
    core, reaction = species_and_reaction()
    if excited_label is not None:
        core[3].label = excited_label
    path = artifact(tmp_path, reaction.kinetics, composition=composition,
                    cold_mean=cold_mean, accepted=accepted)
    reaction.index = reaction_index
    reaction.entry.index = entry_index
    config = {
        'provider': 'loki-table', 'table': (str(path), file_hash(path / 'table.h5')),
        'branch': 'branch_0', 'initial_reduced_field': (np.e, 'Td'),
        'empirical_laws': dict(empirical_laws or {}),
    }
    initial = {core[0]: 1.e-6, core[2]: 1.e-6, core[3]: .1 if composition else 0.,
               core[1]: .899998 if composition else .999998}
    if extra_species is not None:
        core.append(extra_species)
        initial[extra_species] = 0.
    reactions = [reaction] + ([] if extra_reaction is None else [extra_reaction(core)])
    energy_balance = {'absorbed_power': (power_w, 'W'),
                      'chamber_volume': (VOLUME, 'm^3'),
                      'sheath': 'floating_wall'}
    if electron_energies is not None:
        energy_balance['electron_energies'] = dict(electron_energies)
    reactor = PlasmaReactor(
        (TGAS, 'K'), (PRESSURE, 'Pa'), initial, (TGAS, 'K'), termination=[],
        diffusion_length=(LAMBDA, 'm'), ion_reduced_mobility=(MU0, 'm^2/(V*s)'),
        wall_single_bath_approximation=extra_species is not None,
        wall_neutralization_products={'Ar+': 'Ar'},
        electron_energy_balance=energy_balance,
        electron_kinetics=config, thermo_source_assertions={'Ar+': 'ion'})
    if development:
        with development_unqualified_table_route():
            reactor.initialize_model(core, reactions, [], [])
    else:
        reactor.initialize_model(core, reactions, [], [])
    return reactor, core, reactions


def test_development_progress_is_timed_and_production_is_untouched(
        tmp_path, monkeypatch, capsys):
    import rmgpy.solver.plasma as plasma_module

    now = [100.]
    monkeypatch.setattr(plasma_module.time, 'monotonic', lambda: now[0])
    development, _, _ = build_reactor(
        tmp_path / 'development', accepted=False, development=True)
    production, _, _ = build_reactor(tmp_path / 'production')

    development.configure_development_run(
        progress_interval_seconds=30., wall_budget_seconds=90.)
    with pytest.raises(PlasmaStateError, match='development-only'):
        production.configure_development_run(
            progress_interval_seconds=30., wall_budget_seconds=90.)

    zeros = np.zeros_like(development.y0)
    development.residual(0., development.y0.copy(), zeros)
    assert capsys.readouterr().err == ''
    now[0] += 30.
    development.residual(1.e-9, development.y0.copy(), zeros)
    line = capsys.readouterr().err.strip()
    assert line.startswith('EEDF_DEVELOPMENT_PROGRESS ')
    progress = json.loads(line.split(' ', 1)[1])
    assert set(progress) >= {
        'wall_seconds', 't_s', 'dt_s', 'step_count',
        'residual_evaluations', 'jacobian_evaluations',
        'u', 'EN_Td', 'electron_density_m^-3', 'Ar4s_total',
        'max_residual', 'max_residual_row',
    }
    assert progress['t_s'] == pytest.approx(1.e-9)
    assert progress['residual_evaluations'] == 2
    assert progress['scientific_status'] == (
        'DEVELOPMENT ONLY — TABLE QUALIFICATION FAILED')
    assert progress['artifact_blockers']
    assert progress['export_allowed'] is False
    assert progress['qualification_allowed'] is False
    assert development.development_last_progress == progress

    production.residual(0., production.y0.copy(), np.zeros_like(production.y0))
    assert capsys.readouterr().err == ''


def test_development_wall_budget_forces_a_last_progress_line(
        tmp_path, monkeypatch, capsys):
    import rmgpy.solver.plasma as plasma_module

    reactor, _, _ = build_reactor(
        tmp_path, accepted=False, development=True)
    accepted_t = reactor.t
    accepted_y = np.array(reactor.y)
    published_budget = reactor.energy_budget
    now = [10.]

    def advancing_clock():
        result = now[0]
        now[0] += 10.
        return result

    monkeypatch.setattr(plasma_module.time, 'monotonic', advancing_clock)
    reactor.configure_development_run(
        progress_interval_seconds=1.e6, wall_budget_seconds=45.)

    with pytest.raises(DevelopmentWallBudgetExceeded):
        reactor.step(1.e-8)

    line = capsys.readouterr().err.strip()
    assert line.startswith('EEDF_DEVELOPMENT_PROGRESS ')
    progress = json.loads(line.split(' ', 1)[1])
    assert progress['wall_seconds'] >= 45.
    assert progress['scientific_status'] == (
        'DEVELOPMENT ONLY — TABLE QUALIFICATION FAILED')
    assert progress['artifact_blockers']
    assert reactor.development_wall_budget_hit is True
    assert reactor.development_last_progress == progress
    assert reactor.t == accepted_t
    np.testing.assert_array_equal(reactor.y, accepted_y)
    assert reactor.energy_budget is published_budget


def test_native_residual_callback_reraises_original_before_publication(
        tmp_path, monkeypatch):
    reactor, _, _ = build_reactor(tmp_path)
    reactor.step(1.e-8)
    accepted_t = reactor.t
    accepted_y = np.array(reactor.y)
    published_budget = reactor.energy_budget
    failure_times = []

    class ResidualCallbackError(ValueError):
        pass

    def raising_row(*args, **kwargs):
        failure_times.append(reactor.t)
        raise ResidualCallbackError('residual callback failed')

    monkeypatch.setattr(reactor.eedf_provider, 'row', raising_row)

    with pytest.raises(ResidualCallbackError, match='residual callback failed'):
        reactor.step(1.e-8)

    assert failure_times
    assert reactor.t == accepted_t
    assert reactor.t <= failure_times[0]
    np.testing.assert_array_equal(reactor.y, accepted_y)
    assert reactor.energy_budget is published_budget
    with pytest.raises(ResidualCallbackError, match='residual callback failed'):
        reactor.step(1.e-8)


def test_native_jacobian_callback_reraises_original_before_publication(
        tmp_path, monkeypatch):
    import rmgpy.solver.plasma as plasma_module

    reactor, _, _ = build_reactor(tmp_path)
    accepted_t = reactor.t
    accepted_y = np.array(reactor.y)
    published_budget = reactor.energy_budget
    original_zeros = np.zeros
    failure_times = []

    class JacobianCallbackError(ArithmeticError):
        pass

    def raising_jacobian_allocation(shape, *args, **kwargs):
        if isinstance(shape, tuple) and shape == (reactor.neq, reactor.neq):
            failure_times.append(reactor.t)
            raise JacobianCallbackError('jacobian callback failed')
        return original_zeros(shape, *args, **kwargs)

    monkeypatch.setattr(plasma_module.np, 'zeros', raising_jacobian_allocation)

    with pytest.raises(JacobianCallbackError, match='jacobian callback failed'):
        reactor.step(1.e-8)

    assert failure_times
    assert reactor.t == accepted_t
    assert reactor.t <= failure_times[0]
    np.testing.assert_array_equal(reactor.y, accepted_y)
    assert reactor.energy_budget is published_budget


def test_residual_and_progress_hot_paths_do_not_copy_the_manifest(
        tmp_path, monkeypatch, capsys):
    import rmgpy.solver.eedf_provider as provider_module
    import rmgpy.solver.plasma as plasma_module

    production, _, _ = build_reactor(tmp_path / 'production')
    development, _, _ = build_reactor(
        tmp_path / 'development', accepted=False, development=True)
    original_deepcopy = copy.deepcopy
    copies = []

    def counted_deepcopy(value, *args, **kwargs):
        if (isinstance(value, dict) and 'artifact_sha256' in value and
                'channel_map' in value and 'branches' in value):
            copies.append(value)
        return original_deepcopy(value, *args, **kwargs)

    monkeypatch.setattr(provider_module.copy, 'deepcopy', counted_deepcopy)
    production.residual(
        production.t, production.y.copy(), np.zeros_like(production.y))
    assert copies == []

    now = [100.]
    monkeypatch.setattr(plasma_module.time, 'monotonic', lambda: now[0])
    development.configure_development_run(progress_interval_seconds=30.)
    development.residual(
        development.t, development.y.copy(), np.zeros_like(development.y))
    copies.clear()
    now[0] += 30.
    development.residual(
        development.t, development.y.copy(), np.zeros_like(development.y))
    assert capsys.readouterr().err.startswith('EEDF_DEVELOPMENT_PROGRESS ')
    assert copies == []


def test_development_trajectory_compacts_without_losing_endpoint_or_extrema(
        tmp_path, monkeypatch):
    import rmgpy.solver.plasma as plasma_module

    reactor, _, _ = build_reactor(
        tmp_path, accepted=False, development=True)
    reactor.configure_development_run(progress_interval_seconds=30.)
    monkeypatch.setattr(plasma_module, 'DEVELOPMENT_TRAJECTORY_RECORD_LIMIT', 4)
    reduced_fields = [1.0, 1.01, 1.02, 1.03, 1.04, 1.2, 1.05]

    for step, u in enumerate(reduced_fields):
        state = reactor.y0.copy()
        state[reactor.te_index] = u
        reactor._latch_energy_budget(state, float(step))

    assert [row['time_s'] for row in reactor.development_trajectory] == [0., 4., 6.]
    assert reactor.development_trajectory[0]['u'] == reduced_fields[0]
    assert reactor.development_trajectory[-1]['u'] == reduced_fields[-1]
    assert reactor.development_transient_extrema['EN_Td'] == pytest.approx({
        'min': np.exp(min(reduced_fields)),
        'max': np.exp(max(reduced_fields)),
    })


def test_table_owned_source_cannot_also_run_declared_temperature_law(tmp_path):
    def duplicate_temperature_law(core):
        return LibraryReaction(
            reactants=[core[0], core[1]],
            products=[core[2], core[0], core[0]],
            reversible=False,
            library='Toy',
            entry=Entry(index=1),
            kinetics=TwoTemperaturePlasma(
                A=(1.e-15, 'm^3/(molecule*s)'), n=0.,
                Ea_g=(0., 'J/mol'), Ea_e=(0., 'J/mol')),
        )

    with pytest.raises(
            PlasmaStateError,
            match='table-owned EEDF channel source Toy:1.*electron-temperature law'):
        build_reactor(
            tmp_path,
            extra_reaction=duplicate_temperature_law,
            empirical_laws={
                'Toy:1': {
                    'evaluate_at': 'Te_eff', 'class': 'A',
                    'basis': 'regression', 'sensitivity': True,
                },
            },
        )


def test_ar4s_envelope_uses_species_identity_after_label_alias(tmp_path):
    reactor, core, _ = build_reactor(
        tmp_path, composition=False, excited_label='ArStar')
    state = np.array(reactor.y0)
    excited = reactor.species_index[core[3]]
    ground = reactor.species_index[core[1]]
    state[excited] = .18
    state[ground] -= .18

    with pytest.raises(PlasmaStateError, match='Ar4s_total'):
        reactor._check_accepted_plasma_domain(state)


def build_eedf_electronegative_reactor(
        tmp_path, monkeypatch, variable_epsilon=False, accepted=True,
        development=False):
    import plasmaElectronegativeWallTest as en
    import rmgpy.solver.plasma as plasma_module

    _, marker_reaction = species_and_reaction()
    pressure = 5. * 101325. / 760.
    path = artifact(
        tmp_path, marker_reaction.kinetics, composition=False,
        reaction_owned=False, gases=['Ar', 'Cl'], pressure=pressure,
        epsilon_k_eV=4., include_ar4s_envelope=False,
        variable_epsilon=variable_epsilon, accepted=accepted)
    config = {
        'provider': 'loki-table',
        'table': (str(path), file_hash(path / 'table.h5')),
        'branch': 'branch_0',
        'initial_reduced_field': (np.e, 'Td'),
        'empirical_laws': {},
    }
    radial_inputs = []
    reference_inputs = []
    actual_radial_source = plasma_module.radial_source_frequency

    def capture_radial(*args):
        radial_inputs.append(args[1])
        return actual_radial_source(*args)

    def capture_reference(context):
        reference_inputs.append(context['electron_temperature'])
        return en.fixture_reference(context)

    monkeypatch.setattr(plasma_module, 'radial_source_frequency', capture_radial)
    kwargs = {
        'electron_kinetics': config,
        'electron_energy_balance': {
            'absorbed_power': (0., 'W'),
            'chamber_volume': (VOLUME, 'm^3'),
            'sheath': 'floating_wall',
        },
        'electronegative_wall_qualification': en.qualification(
            full_profile_reference=capture_reference,
            geometry_reference=capture_reference),
    }
    if development:
        with development_unqualified_table_route():
            reactor, species, reactions = en.model(**kwargs)
    else:
        reactor, species, reactions = en.model(**kwargs)
    return reactor, species, reactions, radial_inputs, reference_inputs


def test_eedf_en_wall_closures_receive_row_epsilon_k(tmp_path, monkeypatch):
    reactor, _, _, radial_inputs, reference_inputs = \
        build_eedf_electronegative_reactor(tmp_path, monkeypatch)
    epsilon_k_temperature = 4. * constants.Na * constants.e / constants.R
    te_eff_eV = reactor.Te.value_si * constants.R / (constants.Na * constants.e)

    assert te_eff_eV == pytest.approx(4. / 3.)
    assert radial_inputs
    assert reference_inputs
    assert radial_inputs == pytest.approx([epsilon_k_temperature] * len(radial_inputs))
    assert reference_inputs == pytest.approx(
        [epsilon_k_temperature] * len(reference_inputs))


def actual_eedf_wall_residual(reactor, candidate):
    """Isolate the wall delta from the integrated residual and its energy ledger."""
    try:
        reactor.has_wall = True
        wall_on = reactor.residual(
            0., candidate, np.zeros_like(candidate))[0].copy()
        terms_on = dict(reactor.electron_energy_terms)
        reactor.has_wall = False
        reactor.wall_loss_rates.fill(0.)
        wall_off = reactor.residual(
            0., candidate, np.zeros_like(candidate))[0].copy()
        terms_off = dict(reactor.electron_energy_terms)
        result = wall_on - wall_off
        denominator = (candidate[reactor.electron_index] * constants.Na *
                       constants.e * reactor.eedf_row.denergy_du)
        power_delta = terms_on['P_abs'] - terms_off['P_abs']
        loss_delta = sum(
            terms_on[key] - terms_off[key]
            for key in ('Q_inelastic', 'Q_elastic', 'Q_wall_electron',
                        'Q_wall_ion')
        )
        storage_delta = constants.Na * constants.e * (
            terms_on['mean_energy_eV'] * terms_on['dNe_dt'] -
            terms_off['mean_energy_eV'] * terms_off['dNe_dt'])
        composition_delta = (
            terms_on['composition_energy_derivative'] -
            terms_off['composition_energy_derivative'])
        result[reactor.te_index] = (
            power_delta - loss_delta - storage_delta - composition_delta
        ) / denominator
        return result
    finally:
        reactor.has_wall = True


@pytest.mark.parametrize('variable_epsilon', [False, True])
def test_public_eedf_en_wall_jacobian_matches_wall_residual(
        tmp_path, monkeypatch, variable_epsilon):
    reactor, _, _, _, _ = build_eedf_electronegative_reactor(
        tmp_path, monkeypatch, variable_epsilon=variable_epsilon)
    state = np.array(reactor.y0)
    analytic = np.asarray(reactor.compute_electronegative_wall_jacobian(state))

    numerical = np.zeros_like(analytic)
    for column in range(len(state)):
        step = 1.e-4 * max(abs(state[column]), reactor.atol_array[column])
        plus, minus = state.copy(), state.copy()
        plus[column] += step
        minus[column] -= step
        numerical[:, column] = (
            actual_eedf_wall_residual(reactor, plus) -
            actual_eedf_wall_residual(reactor, minus)
        ) / (2. * step)
    reactor.residual(0., state, np.zeros_like(state))

    assert numerical[reactor.te_index, reactor.te_index] != 0.
    np.testing.assert_allclose(analytic, numerical, rtol=2.e-3, atol=2.e-6)


@pytest.mark.parametrize('endpoint', [0, -1])
def test_public_eedf_en_wall_jacobian_uses_inward_difference_at_u_endpoint(
        tmp_path, monkeypatch, endpoint):
    reactor, _, _, _, _ = build_eedf_electronegative_reactor(
        tmp_path, monkeypatch, variable_epsilon=True)
    state = np.array(reactor.y0)
    state[reactor.te_index] = reactor.eedf_provider.axis('u')[endpoint]
    jacobian = np.asarray(reactor.compute_electronegative_wall_jacobian(state))
    u_axis = reactor.eedf_provider.axis('u')
    step = 6.e-6 * max(
        abs(state[reactor.te_index]), float(u_axis[-1] - u_axis[0]))
    inside = state.copy()
    direction = 1. if endpoint == 0 else -1.
    inside[reactor.te_index] += direction * step
    expected = direction * (
        actual_eedf_wall_residual(reactor, inside) -
        actual_eedf_wall_residual(reactor, state)
    ) / step

    np.testing.assert_allclose(
        jacobian[:, reactor.te_index], expected, rtol=2.e-3, atol=2.e-6)


def test_development_en_wall_manifest_carries_unqualified_notice(
        tmp_path, monkeypatch):
    reactor, _, _, _, _ = build_eedf_electronegative_reactor(
        tmp_path, monkeypatch, accepted=False, development=True)
    notice = reactor.eedf_provider.development_notice
    manifest = reactor.electronegative_wall_manifest()

    assert notice is not None
    assert {key: manifest[key] for key in notice} == notice
    assert manifest['closure_scientific_status'] == (
        'FALSIFIED SIMPLIFIED CLOSURE')
    assert reactor.electronegative_wall_diagnostics[
        'closure_scientific_status'] == 'FALSIFIED SIMPLIFIED CLOSURE'
    assert {key: reactor.electronegative_wall_diagnostics[key]
            for key in notice} == notice


def test_development_run_labels_provider_diagnostics_budget_and_manifest(tmp_path):
    with pytest.raises(PlasmaStateError, match='held-out qualification'):
        build_reactor(tmp_path / 'production', accepted=False)

    reactor, core, reactions = build_reactor(
        tmp_path / 'development', accepted=False, development=True)
    expected = {
        'scientific_status': DEVELOPMENT_UNQUALIFIED_STATUS,
        'artifact_blockers': [
            'artifact manifest is not accepted',
            '1 of 1 held-out verdicts failed',
            'branch_0: fixture',
        ],
        'export_allowed': False,
        'qualification_allowed': False,
    }
    assert {key: reactor.eedf_provider.manifest[key] for key in expected} == expected
    assert {key: reactor.electron_energy_terms[key] for key in expected} == expected
    assert {key: reactor.energy_budget[key] for key in expected} == expected
    assert {key: reactor.eedf_run_manifest()[key] for key in expected} == expected

    # The process-local route is not serialized: the reconstructed reactor
    # returns to the production loader and refuses the same artifact.
    restored = copy.deepcopy(reactor)
    restored_core = list(restored.initial_mole_fractions)
    restored_by_label = {species.label: species for species in restored_core}
    restored_reactions = [LibraryReaction(
        reactants=[restored_by_label['e-'], restored_by_label['Ar']],
        products=[restored_by_label['Ar+'], restored_by_label['e-'],
                  restored_by_label['e-']],
        reversible=False, kinetics=copy.deepcopy(reactions[0].kinetics),
        library='Toy', entry=Entry(index=1))]
    with pytest.raises(PlasmaStateError, match='held-out qualification'):
        restored.initialize_model(restored_core, restored_reactions, [], [])


def test_input_deck_naming_unaccepted_artifact_still_refuses(tmp_path):
    core, reaction = species_and_reaction()
    path = artifact(tmp_path / 'table', reaction.kinetics, composition=False,
                    accepted=False)
    digest = file_hash(path / 'table.h5')
    input_path = tmp_path / 'input.py'
    deck = """
species(label='e-', reactive=True, structure=adjacencyList('1 e u1 p0 c-1'))
species(label='Ar', reactive=True, structure=adjacencyList('1 Ar u0 p4 c0'))
species(label='Arp', reactive=True,
        structure=adjacencyList('multiplicity 2\\n1 Ar u1 p3 c+1'))
species(label='Ars', reactive=True,
        structure=adjacencyList('multiplicity 3\\n1 Ar u2 p3 c0'))
plasmaReactor(
    temperature=(298.15, 'K'), pressure=(666.61, 'Pa'),
    initialMoleFractions={'e-': 1e-6, 'Ar': .999998, 'Arp': 1e-6, 'Ars': 0.0},
    chamberGeometry={'diffusionLength': (.02, 'm')},
    ionReducedMobility=(1.535e-4, 'm^2/(V*s)'),
    wallNeutralizationProducts={'Arp': 'Ar'},
    thermoSourceAssertions={'Arp': 'ion'},
    electronKinetics=__CONFIG__,
    electronEnergyBalance={'absorbedPower': (.5, 'W'),
                           'chamberVolume': (1e-4, 'm^3'),
                           'sheath': 'floatingWall'},
    terminationTime=(1, 's'))
"""
    input_path.write_text(deck.replace('__CONFIG__', repr({
        'provider': 'loki-table', 'table': (str(path / 'table.h5'), digest),
        'branch': 'branch_0', 'initialReducedField': (np.e, 'Td'),
        'empiricalLaws': {},
    })))
    rmg = RMG()
    rmg_input.read_input_file(str(input_path), rmg)
    reactor = rmg.reaction_systems[0]
    by_label = {species.label: species for species in reactor.initial_mole_fractions}
    by_label['Ar'].thermo = thermo()
    by_label['Arp'].thermo = thermo(15.76)
    by_label['Ars'].thermo = thermo(11.55)
    marker = EEDFChannel(process=reaction.kinetics.process,
                         collision_set='fixture', side='ine')
    deck_reaction = LibraryReaction(
        reactants=[by_label['e-'], by_label['Ar']],
        products=[by_label['Arp'], by_label['e-'], by_label['e-']],
        reversible=False, kinetics=marker, library='Toy', entry=Entry(index=1))
    with pytest.raises(PlasmaStateError, match='held-out qualification'):
        reactor.initialize_model(list(by_label.values()), [deck_reaction], [], [])


def test_u_state_atol_restart_and_mapped_rate(tmp_path):
    reactor, _, _ = build_reactor(tmp_path)
    assert reactor.y0[reactor.te_index] == pytest.approx(1.)
    assert reactor.atol_array[reactor.te_index] == 1.e-8
    assert reactor.kf[0] == pytest.approx(constants.Na * reactor.eedf_row.k_ine[0])
    args = reactor.__reduce__()[1]
    assert args[-1] == reactor.electron_kinetics
    assert copy.deepcopy(reactor).electron_kinetics == reactor.electron_kinetics
    manifest = reactor.eedf_run_manifest()
    assert manifest['wall_sheath_basis'] == 'Maxwellian form at epsilon_k'
    assert manifest['branch'] == 'branch_0'


def test_eedf_requires_the_energy_balance_at_construction():
    config = {'provider': 'loki-table', 'table': ('unused', '0' * 64),
              'branch': 'branch_0', 'initial_reduced_field': (1., 'Td'),
              'empirical_laws': {}}
    with pytest.raises(PlasmaStateError, match='electron_energy_balance is required'):
        PlasmaReactor((TGAS, 'K'), (PRESSURE, 'Pa'), {}, (TGAS, 'K'),
                      electron_kinetics=config)


def test_full_coordinate_cache_eos_and_transport_move_together(tmp_path):
    reactor, core, _ = build_reactor(tmp_path)
    y = np.array(reactor.y0)
    first = reactor.eedf_provider.cache_identity
    rate = reactor.kf[0]
    volume = reactor.compute_volume(y)
    nu = reactor.compute_nu_wall(y, volume)
    reactor.residual(0., y, np.zeros_like(y))
    power = reactor.electron_energy_terms['Q_inelastic']
    excited = reactor.species_index[core[3]]
    ground = reactor.species_index[core[1]]
    y[excited] += .04
    y[ground] -= .04
    reactor.residual(0., y, np.zeros_like(y))
    assert reactor.eedf_provider.cache_identity != first
    assert reactor.kf[0] != rate
    assert reactor.compute_volume(y) != volume
    assert reactor.compute_nu_wall(y, reactor.compute_volume(y)) != nu
    assert reactor.electron_energy_terms['Q_inelastic'] != power


def test_elastic_row_power_is_a_positive_loss_and_steady_state_gates_run(tmp_path):
    reactor, _, _ = build_reactor(tmp_path)
    y = np.array(reactor.y0)
    reactor.residual(0., y, np.zeros_like(y))
    assert reactor.electron_energy_terms['Q_elastic'] > 0.
    reactor.energy_budget = {'A6a_relative': 0.0, 'A6a_steady_relative': 0.0,
                             'A6b_relative': 0.0,
                             'A6b_tolerance': 1.e-6}
    reactor.validate_steady_state()
    reactor.energy_budget['A6a_steady_relative'] = 0.02
    with pytest.raises(PlasmaStateError, match='A6a'):
        reactor.validate_steady_state()


def test_terminal_state_qualifies_once_and_records_the_outcome(tmp_path, monkeypatch):
    reactor, _, _ = build_reactor(tmp_path)
    reactor.energy_budget = {'A6a_relative': 0., 'A6a_steady_relative': 0.,
                             'A6b_relative': 0., 'A6b_tolerance': 1.e-6}
    calls = []

    def qualify(y, t, *, runner=None):
        calls.append((np.array(y), t))
        return {'status': 'PASS', 'record': 'fixture'}

    monkeypatch.setattr(reactor.eedf_provider, 'qualify', qualify)
    reactor.validate_terminal_state()
    reactor.validate_terminal_state()

    assert len(calls) == 1
    assert reactor.energy_terminal['qualification'] == 'PASS'
    assert reactor.energy_terminal['qualification_record']['record'] == 'fixture'


def test_terminal_qualification_failure_marks_the_state_before_refusing(tmp_path, monkeypatch):
    reactor, _, _ = build_reactor(tmp_path)
    reactor.energy_budget = {'A6a_relative': 0., 'A6a_steady_relative': 0.,
                             'A6b_relative': 0., 'A6b_tolerance': 1.e-6}

    calls = []

    def refuse(y, t, *, runner=None):
        calls.append((y, t))
        raise PlasmaStateError('EEDF qualification failed: A3')

    monkeypatch.setattr(reactor.eedf_provider, 'qualify', refuse)
    with pytest.raises(PlasmaStateError, match='EEDF qualification failed: A3'):
        reactor.validate_terminal_state()
    assert reactor.energy_terminal['qualification'] == 'FAILED'
    with pytest.raises(PlasmaStateError, match='latched FAILED'):
        reactor.validate_terminal_state()
    assert len(calls) == 1


def test_terminal_power_refusal_replaces_an_earlier_pass_and_blocks_export(
        tmp_path):
    reactor, _, _ = build_reactor(tmp_path)
    provider = reactor.eedf_provider
    record_path = tmp_path / 'eedf_qualification.json'
    state = {
        'u': 1., 'composition': {'x(Ars)': 0.},
        'record_path': str(record_path),
        'gas_fractions': {'Ar': 1.},
        'state_populations': {'Ar': 1.},
        'state_statistical_weights': {'Ar': 1.},
    }
    row = provider.row(state['u'], state['composition'])
    state['energy_budget'] = {
            'P_abs': float(row.power_groups['field']),
            'joule_power': float(row.power_groups['field']),
            'A6b_relative': 0.,
            'A6b_tolerance': 1.e-6,
    }
    provider.bind_qualification_state(lambda y, t: state)
    provider._qualification_status = 'PASS'
    reactor.energy_terminal = {
        'termination': 'accepted', 'qualification': 'PASS',
        'qualification_record': {'status': 'PASS'},
    }
    provider.require_qualified('export')
    reactor.energy_budget = {
        'A6a_relative': .02, 'A6a_steady_relative': .02,
        'A6b_relative': 0., 'A6b_tolerance': 1.e-6,
    }

    with pytest.raises(PlasmaStateError, match='A6a'):
        reactor.validate_terminal_state()

    assert reactor.energy_terminal['qualification'] == 'FAILED'
    assert reactor.energy_terminal['qualification_record']['status'] == 'FAIL'
    assert json.loads(record_path.read_text())['status'] == 'FAIL'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_extinction_power_refusal_runs_inside_terminal_transaction(tmp_path):
    reactor, _, _ = build_reactor(tmp_path, cold_mean=True)
    provider = reactor.eedf_provider
    record_path = tmp_path / 'eedf_qualification.json'
    state = {
        'u': 0., 'composition': {'x(Ars)': 0.},
        'record_path': str(record_path),
        'gas_fractions': {'Ar': 1.},
        'state_populations': {'Ar': 1.},
        'state_statistical_weights': {'Ar': 1.},
    }
    row = provider.row(state['u'], state['composition'])
    state['energy_budget'] = {
        'P_abs': float(row.power_groups['field']),
        'joule_power': float(row.power_groups['field']),
        'A6b_relative': 0.,
        'A6b_tolerance': 1.e-6,
    }
    provider.bind_qualification_state(lambda y, t: state)
    provider._qualification_status = 'PASS'
    provider.require_qualified('export')

    y = np.array(reactor.y0)
    y[reactor.te_index] = 0.
    reactor._check_accepted_plasma_domain(y)
    reactor.energy_budget.update(
        discharge_state='extinct', nu_ionisation=0., nu_source=0., nu_loss=1.e99,
        A6a_relative=.02, A6a_steady_relative=.02,
        A6b_relative=0., A6b_tolerance=1.e-6)
    reactor.energy_was_self_sustained = True
    reactor._update_terminal_state(y, 0.)
    reactor._update_terminal_state(y, reactor.extinction_persistence_time(y))

    assert reactor.energy_terminal['termination'] == 'extinct'
    with pytest.raises(PlasmaStateError, match='A6a'):
        reactor.validate_terminal_state()

    assert reactor.energy_terminal['qualification'] == 'FAILED'
    assert reactor.energy_terminal['qualification_record']['status'] == 'FAIL'
    assert json.loads(record_path.read_text())['status'] == 'FAIL'
    with pytest.raises(EEDFError, match='FAILED.*export'):
        provider.require_qualified('export')


def test_development_diagnostic_accepts_an_explicit_legacy_state_map(
        tmp_path, monkeypatch):
    reactor, core, _ = build_reactor(
        tmp_path, accepted=False, development=True)
    state_map = []
    for index, species in enumerate(core):
        molecule = species.molecule[0]
        if molecule.get_net_charge() != 0:
            continue
        state_map.append({
            'formula': molecule.get_formula(),
            'electronic_state': str(molecule.electronic_state or ''),
            'vibrational_level': molecule.vibrational_level,
            'multiplicity': molecule.multiplicity,
            'loki_state': 'Ar(fixture-{0})'.format(index),
            'statistical_weight': float(molecule.multiplicity),
        })

    captured = {}

    def diagnostic(state, t, *, runner=None):
        captured.update(state)
        return {'status': 'DEVELOPMENT — not a qualification'}

    monkeypatch.setattr(reactor.eedf_provider, 'development_diagnostic', diagnostic)
    record = reactor.development_eedf_diagnostic(engine_state_map=state_map)

    assert record['status'] == 'DEVELOPMENT — not a qualification'
    assert captured['state_populations']
    assert captured['state_statistical_weights']


def test_qualification_state_keeps_a_mapped_zero_population_gas_absent(
        tmp_path):
    helium = Species(label='He').from_adjacency_list('1 He u0 p1 c0')
    helium.thermo = thermo()
    reactor, core, _ = build_reactor(tmp_path, extra_species=helium)
    state_map = []
    for index, species in enumerate(core):
        molecule = species.molecule[0]
        if molecule.get_net_charge() != 0:
            continue
        formula = molecule.get_formula()
        state_map.append({
            'formula': formula,
            'electronic_state': str(molecule.electronic_state or ''),
            'vibrational_level': molecule.vibrational_level,
            'multiplicity': molecule.multiplicity,
            'loki_state': '{0}(fixture-{1})'.format(formula, index),
            'statistical_weight': float(molecule.multiplicity),
        })

    state = reactor._eedf_qualification_state(
        reactor.y0, 0., engine_state_map=state_map)

    assert state['gas_fractions']['He'] == 0.
    assert state['state_populations'][
        next(entry['loki_state'] for entry in state_map
             if entry['formula'] == 'He')] == 0.


def test_unowned_temperature_law_is_refused(tmp_path):
    def extra(core):
        kinetics = TwoTemperaturePlasma(A=(1.e-18, 'm^3/(molecule*s)'), n=0.,
                                        Ea_g=(0., 'J/mol'), Ea_e=(0., 'J/mol'))
        return LibraryReaction(reactants=[core[0], core[2]], products=[core[1]],
                               reversible=False, kinetics=kinetics, library='Empirical')
    with pytest.raises(PlasmaStateError, match='unowned electron-temperature law'):
        build_reactor(tmp_path, extra_reaction=extra)


def test_source_entry_identity_survives_model_renumbering_and_mismatch_is_refused(tmp_path):
    reactor, _, _ = build_reactor(tmp_path / 'renumbered', reaction_index=29)
    assert reactor.eedf_reaction_map == {0: (0, 'ine')}
    with pytest.raises(PlasmaStateError, match='reaction identity mismatch'):
        build_reactor(tmp_path / 'mismatch', entry_index=2)


def test_separate_superelastic_entry_maps_the_same_collision_sup_column(tmp_path):
    def superelastic(core):
        marker = EEDFChannel(process='e + Ar -> e + e + Ar+',
                             collision_set='fixture', side='sup')
        return LibraryReaction(
            reactants=[core[0], core[2]], products=[core[0], core[1]],
            reversible=False, kinetics=marker, library='Toy', entry=Entry(index=2))

    reactor, _, _ = build_reactor(tmp_path, extra_reaction=superelastic)
    assert reactor.eedf_reaction_map == {0: (0, 'ine'), 1: (0, 'sup')}
    reactor.residual(0., reactor.y0, np.zeros_like(reactor.y0))
    terms = reactor.electron_energy_terms
    assert terms['Q_inelastic_by_reaction_gross'][0] > 0.
    assert terms['Q_superelastic_by_reaction'][1] > 0.
    assert terms['Q_inelastic_by_reaction'][0] > 0.
    assert terms['Q_inelastic_by_reaction'][1] < 0.


def test_duplicate_superelastic_binding_is_refused(tmp_path):
    def duplicates(core):
        marker = EEDFChannel(process='e + Ar -> e + e + Ar+',
                             collision_set='fixture', side='sup')
        first = LibraryReaction(reactants=[core[0], core[2]], products=[core[0], core[1]],
                                reversible=False, kinetics=marker, library='Toy',
                                entry=Entry(index=2))
        second = LibraryReaction(reactants=[core[0], core[2]], products=[core[0], core[1]],
                                 reversible=False, kinetics=copy.deepcopy(marker), library='Toy',
                                 entry=Entry(index=3))
        return [first, second]

    core, reaction = species_and_reaction()
    path = artifact(tmp_path, reaction.kinetics)
    config = {'provider': 'loki-table', 'table': (str(path), file_hash(path / 'table.h5')),
              'branch': 'branch_0', 'initial_reduced_field': (np.e, 'Td'),
              'empirical_laws': {}}
    initial = {core[0]: 1.e-6, core[2]: 1.e-6, core[3]: .1, core[1]: .899998}
    reactor = PlasmaReactor(
        (TGAS, 'K'), (PRESSURE, 'Pa'), initial, (TGAS, 'K'), termination=[],
        diffusion_length=(LAMBDA, 'm'), ion_reduced_mobility=(MU0, 'm^2/(V*s)'),
        wall_neutralization_products={'Ar+': 'Ar'},
        electron_energy_balance={'absorbed_power': (.5, 'W'),
                                 'chamber_volume': (VOLUME, 'm^3'),
                                 'sheath': 'floating_wall'},
        electron_kinetics=config, thermo_source_assertions={'Ar+': 'ion'})
    with pytest.raises(PlasmaStateError, match='double-count'):
        reactor.initialize_model(core, [reaction] + duplicates(core), [], [])


def test_declared_heavy_particle_electron_energy_is_kept_and_table_owned_is_refused(tmp_path):
    def empirical(core):
        kinetics = TwoTemperaturePlasma(A=(1.e-18, 'm^3/(molecule*s)'), n=0.,
                                        Ea_g=(0., 'J/mol'), Ea_e=(0., 'J/mol'))
        return LibraryReaction(reactants=[core[0], core[3]], products=[core[0], core[1]],
                               reversible=False, kinetics=kinetics, library='Empirical',
                               entry=Entry(index=90))

    law = {'Empirical:90': {'evaluate_at': 'Te_eff', 'class': 'C',
                            'basis': 'fixture', 'sensitivity': True}}
    reactor, _, _ = build_reactor(
        tmp_path / 'heavy', extra_reaction=empirical, empirical_laws=law,
        electron_energies={'Empirical:90': (-7.34, 'eV')})
    reactor.residual(0., reactor.y0, np.zeros_like(reactor.y0))
    assert reactor.energy_reaction_keys[1] == 'Empirical:90'
    assert reactor.electron_energy_terms['Q_heavy_particle_electron_by_reaction'][1] < 0.

    with pytest.raises(PlasmaStateError, match='table-owned'):
        build_reactor(tmp_path / 'owned', electron_energies={'Toy:1': (15.76, 'eV')})


def test_new_neutral_additive_is_refused_at_an_accepted_state(tmp_path):
    helium = Species(label='He').from_adjacency_list('1 He u0 p1 c0')
    helium.thermo = thermo()
    reactor, core, _ = build_reactor(tmp_path, extra_species=helium)
    y = np.array(reactor.y0)
    y[reactor.species_index[helium]] = .01
    y[reactor.species_index[next(spc for spc in core if spc.label == 'Ar')]] -= .01
    with pytest.raises(PlasmaStateError, match='additive identity'):
        reactor._check_accepted_plasma_domain(y)


def test_eedf_steady_state_hooks_reattach_u_not_effective_temperature(tmp_path):
    reactor, _, _ = build_reactor(tmp_path)
    species_state = np.array(reactor.y0[:reactor.num_core_species])
    assert reactor.steady_state_relaxation_time(1., species_state) > 0.
    assert reactor.eedf_row.u == pytest.approx(reactor.y[reactor.te_index])


def test_trial_clamps_but_accepted_domain_refuses(tmp_path):
    reactor, _, _ = build_reactor(tmp_path)
    trial = np.array(reactor.y0)
    trial[reactor.te_index] = -100.
    reactor.residual(0., trial, np.zeros_like(trial))
    assert reactor.eedf_row.u == 0.
    with pytest.raises(PlasmaStateError, match='outside the table domain'):
        reactor._check_energy_state(trial)

    envelope = np.array(reactor.y0)
    excited = reactor.species_index[next(spc for spc in reactor.initial_mole_fractions
                                         if spc.label == 'Ars')]
    ground = reactor.species_index[next(spc for spc in reactor.initial_mole_fractions
                                        if spc.label == 'Ar')]
    envelope[excited] = .18
    envelope[ground] -= .08
    with pytest.raises(PlasmaStateError, match='Ar4s_total'):
        reactor._check_energy_state(envelope)


def test_energy_chain_rule_and_central_difference_jacobian(tmp_path):
    reactor, _, _ = build_reactor(tmp_path)
    y = np.array(reactor.y0)
    dydt = np.zeros_like(y)
    delta = reactor.residual(0., y, dydt)[0]
    terms = reactor.electron_energy_terms
    ne = y[reactor.electron_index]
    stored = (constants.Na * constants.e *
              (terms['mean_energy_eV'] * delta[reactor.electron_index] +
               ne * reactor.eedf_row.denergy_du * delta[reactor.te_index]) +
              terms['composition_energy_derivative'])
    losses = sum(terms[name] for name in
                 ('Q_inelastic', 'Q_elastic', 'Q_wall_electron', 'Q_wall_ion', 'Q_flow'))
    assert stored == pytest.approx(terms['P_abs'] - losses, rel=2.e-8)

    pd = np.asarray(reactor.jacobian(0., y, dydt, 0.))
    for column in (reactor.species_index[next(sp for sp in reactor.initial_mole_fractions
                                              if sp.label == 'Ars')], reactor.te_index):
        h = 5.e-6 * max(abs(y[column]), reactor.atol_array[column])
        yp, ym = y.copy(), y.copy()
        yp[column] += h
        ym[column] -= h
        fd = ((reactor.residual(0., yp, dydt)[0] - reactor.residual(0., ym, dydt)[0]) /
              (2. * h))
        scale = np.max(np.abs(fd)) or 1.
        assert np.allclose(pd[:, column], fd, rtol=3.e-4, atol=2.e-6 * scale)


def test_terminal_a6_gates_pass_and_fail_by_name(tmp_path):
    reactor, _, _ = build_reactor(tmp_path)
    reactor.energy_budget.update(A6a_relative=.01, A6b_relative=5.e-7,
                                 A6b_tolerance=1.e-6)
    reactor._check_eedf_power_gates()
    reactor.energy_budget['A6a_relative'] = .0100001
    with pytest.raises(PlasmaStateError, match='A6a.*tolerance'):
        reactor._check_eedf_power_gates()
    reactor.energy_budget['A6a_relative'] = 0.
    reactor.energy_budget['A6b_relative'] = 1.0001e-6
    with pytest.raises(PlasmaStateError, match='A6b.*tolerance'):
        reactor._check_eedf_power_gates()


def test_only_development_route_continues_past_recorded_a6b_failure(tmp_path):
    production, _, _ = build_reactor(tmp_path / 'production')
    development, _, _ = build_reactor(
        tmp_path / 'development', accepted=False, development=True)
    failed = {
        'A6a_relative': 0.0,
        'A6a_steady_relative': 0.0,
        'A6b_relative': 3.1240246954502026e-6,
        'A6b_tolerance': 1.0294419182226285e-6,
        'A6b_numerator': 0.002454087231,
        'A6b_denominator': 785.288329393846,
    }
    production.energy_budget.update(failed)
    development.energy_budget.update(failed)

    with pytest.raises(PlasmaStateError, match='A6b.*tolerance'):
        production._check_eedf_power_gates(steady=True)

    development.configure_development_run(progress_interval_seconds=30.)
    development._check_eedf_power_gates(steady=True)

    record = development.development_a6b_failure
    assert record['outcome'] == 'FAIL'
    assert record['value'] == failed['A6b_relative']
    assert record['tolerance'] == failed['A6b_tolerance']
    assert record['numerator'] == failed['A6b_numerator']
    assert record['denominator'] == failed['A6b_denominator']
    assert record['artifact_sha256'] == development.electron_kinetics['table'][1]
    assert record['scientific_status'] == (
        'DEVELOPMENT ONLY \u2014 TABLE QUALIFICATION FAILED')
    assert record['export_allowed'] is False
    assert record['qualification_allowed'] is False


@pytest.mark.parametrize('field,value', [
    ('A6b_relative', 'missing'),
    ('A6b_relative', np.nan),
    ('A6b_relative', np.inf),
    ('A6b_relative', -1.e-6),
    ('A6b_tolerance', 'missing'),
    ('A6b_tolerance', np.nan),
    ('A6b_tolerance', np.inf),
    ('A6b_tolerance', 0.0),
    ('A6b_tolerance', -1.e-6),
    ('A6b_numerator', 'missing'),
    ('A6b_numerator', np.nan),
    ('A6b_numerator', np.inf),
    ('A6b_numerator', -1.0),
    ('A6b_denominator', 'missing'),
    ('A6b_denominator', np.nan),
    ('A6b_denominator', np.inf),
    ('A6b_denominator', 0.0),
    ('A6b_denominator', -1.0),
])
def test_development_route_refuses_invalid_a6b_diagnostics(tmp_path, field, value):
    reactor, _, _ = build_reactor(tmp_path, accepted=False, development=True)
    diagnostics = {
        'A6a_relative': 0.0,
        'A6a_steady_relative': 0.0,
        'A6b_relative': 3.e-6,
        'A6b_tolerance': 1.e-6,
        'A6b_numerator': 0.002,
        'A6b_denominator': 700.0,
    }
    if value != 'missing':
        diagnostics[field] = value
    reactor.energy_budget.update(diagnostics)
    if value == 'missing':
        reactor.energy_budget.pop(field, None)
    reactor.configure_development_run(progress_interval_seconds=30.)

    with pytest.raises(PlasmaStateError, match='A6b.*tolerance'):
        reactor._check_eedf_power_gates(steady=True)
    assert reactor.development_a6b_failure is None


def test_extinction_uses_row_effective_temperature_and_elastic_frequency(tmp_path):
    reactor, _, _ = build_reactor(tmp_path, cold_mean=True)
    y = np.array(reactor.y0)
    y[reactor.te_index] = 0.0
    reactor._check_accepted_plasma_domain(y)
    reactor.energy_budget.update(
        discharge_state='extinct', nu_ionisation=0., nu_source=0., nu_loss=1.e99)
    reactor.energy_was_self_sustained = True

    volume = reactor.compute_volume(y)
    neutral_density = sum(y[j] for j in range(reactor.num_core_species)
                          if reactor.neutral_heavy_mask[j]) * constants.Na / volume
    elastic = reactor.eedf_provider.manifest['channel_map'][1]
    relaxation_rate = (2. * elastic['mass_ratio'] *
                       reactor.eedf_row.target_fractions[1] *
                       reactor.eedf_row.k_ine[1] * neutral_density)
    assert reactor.extinction_persistence_time(y) == pytest.approx(2. / relaxation_rate)

    reactor._update_terminal_state(y, 0.)
    assert reactor.Te.value_si == pytest.approx(TGAS)
    assert reactor.energy_extinct_since == 0.


def test_prescribed_power_moves_the_fixed_composition_u_root(tmp_path):
    """C1 section 9: the smallest steady energy-coordinate harness.

    Species amounts are held fixed; only the EEDF energy coordinate is solved.
    This isolates the prescribed-power selection of u from chemical transients.
    """
    roots = []
    for name, power in (('half', 1.e4), ('base', 2.e4), ('double', 4.e4)):
        reactor, _, _ = build_reactor(tmp_path / name, power_w=power)
        y = np.array(reactor.y0)
        dydt = np.zeros_like(y)

        def energy_residual(u):
            state = y.copy()
            state[reactor.te_index] = u
            return reactor.residual(0., state, dydt)[0][reactor.te_index]

        nodes = reactor.eedf_provider.axis('u')
        values = [energy_residual(float(u)) for u in nodes]
        bracket = next((pair for pair in zip(nodes[:-1], nodes[1:], values[:-1], values[1:])
                        if pair[2] * pair[3] <= 0.), None)
        assert bracket is not None
        roots.append(brentq(energy_residual, bracket[0], bracket[1]))
    assert roots[0] < roots[1] < roots[2]
