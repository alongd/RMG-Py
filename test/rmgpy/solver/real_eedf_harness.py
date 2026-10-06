"""Real LoKI-B fixture and pure-Ar EEDF run harness support.

This module owns no scientific inputs.  It reads the spec named by
``EEDF_REAL_SPEC``, verifies that every external file exists, and lets the
existing generator qualify the artifact.  It never edits acceptance,
fingerprints, branch metadata, classifications, or manifests.
"""

import copy
import json
import math
from pathlib import Path

import h5py
import numpy as np
from scipy.interpolate import PchipInterpolator

from rmgpy import constants
from rmgpy.kinetics import EEDFChannel
from rmgpy.solver.eedf_provider import unqualified_artifact_blockers

from rmgpy.tools.eedf.generate import generate
from rmgpy.tools.eedf.schema import file_hash, load_spec


DEVELOPMENT_COLLISION_SET = 'IST-Lisbon-Ar-2021'
DEVELOPMENT_BINDINGS = {
    86: 'e + Ar(1S0) -> e + e + Ar(+,gnd), Ionization',
    87: 'e + Ar(1S0) -> e + Ar(3P2), Excitation',
}
DEVELOPMENT_RETAINED_ENTRIES = {86, 87, 88, 90, 91, 92, 93}
DEVELOPMENT_SEED_RULING = (
    '/home/alon/runs/owner-rulings/2026-10-06-dev-deck-initial-ne.md')
DEVELOPMENT_ABSORBED_POWER_W = 0.5
DEVELOPMENT_CHAMBER_VOLUME_M3 = 2.356194490192345e-3
DEVELOPMENT_INITIAL_FIELD_TD = 17.0
DEVELOPMENT_SWEEP_TOLERANCES = {
    'electron_density_m^-3': {'rtol': 2.e-4, 'atol': 0.0},
    'EN_Td': {'rtol': 2.e-5, 'atol': 0.0},
    'mean_energy_eV': {'rtol': 2.e-5, 'atol': 0.0},
    'composition.Ar': {'rtol': 0.0, 'atol': 2.e-8},
    'composition.Arp': {'rtol': 2.e-4, 'atol': 1.e-16},
    'wall_flux_mol_s.e-': {'rtol': 2.e-4, 'atol': 1.e-18},
    'wall_flux_mol_s.Arp': {'rtol': 2.e-4, 'atol': 1.e-18},
    'electron_power_partition_W.Q_inelastic': {'rtol': 2.e-4, 'atol': 1.e-10},
    'electron_power_partition_W.Q_elastic': {'rtol': 2.e-4, 'atol': 1.e-10},
    'electron_power_partition_W.Q_wall_electron': {'rtol': 2.e-4, 'atol': 1.e-10},
    'electron_power_partition_W.Q_wall_ion': {'rtol': 2.e-4, 'atol': 1.e-10},
    'A6a_relative': {'rtol': 0.0, 'atol': 1.e-8},
    'A6b_relative': {'rtol': 2.e-4, 'atol': 1.e-12},
}


def _interpolate_initial_row(artifact, branch, field_td):
    """Read the same shape-preserving PCHIP-in-u row used by the provider."""
    with h5py.File(Path(artifact) / 'table.h5', 'r') as table:
        rows = table['branches'][branch]
        u_axis = np.asarray(rows['u'], dtype=float)
        u = math.log(field_td)
        if not u_axis[0] <= u <= u_axis[-1]:
            raise ValueError('initial reduced field is outside the table u axis')
        upper = int(np.searchsorted(u_axis, u, side='right'))
        upper = min(max(upper, 1), len(u_axis) - 1)
        lower = upper - 1

        def interpolate(path):
            return float(PchipInterpolator(
                u_axis, np.asarray(rows[path], dtype=float),
                extrapolate=False)(u))

        return {
            'u': u,
            'bracketing_row_indices': [lower, upper],
            'bracketing_EN_Td': [float(math.exp(u_axis[lower])),
                                  float(math.exp(u_axis[upper]))],
            'field_power_coefficient_eV_m3_s^-1': interpolate('power_groups/field'),
            'mobility_N': interpolate('swarm/mobility_N'),
            'mean_energy_eV': interpolate('swarm/mean_energy_eV'),
        }


def power_consistent_initialization(artifact, manifest, multiplier=1.0):
    """Derive the balanced pure-Ar seed from power, gas density and the table row.

    The reference density equates chamber power density to the table row's
    field-power coefficient.  The returned mole fraction then solves the
    reactor's two-temperature EOS exactly, rather than rounding to 1.7e-10.
    """
    multiplier = float(multiplier)
    if not math.isfinite(multiplier) or multiplier <= 0.0:
        raise ValueError('seed multiplier must be finite and positive')
    branch = manifest['branches'][0]
    pressure = float(manifest['row_inputs']['P_Pa'])
    gas_temperature = float(manifest['row_inputs']['Tg_K'])
    row = _interpolate_initial_row(
        artifact, branch, DEVELOPMENT_INITIAL_FIELD_TD)
    gas_density = pressure / (constants.kB * gas_temperature)
    power_density = DEVELOPMENT_ABSORBED_POWER_W / DEVELOPMENT_CHAMBER_VOLUME_M3
    loss_scale = row['field_power_coefficient_eV_m3_s^-1']
    transport_field = (row['mobility_N'] *
                       (DEVELOPMENT_INITIAL_FIELD_TD * 1.e-21) ** 2)
    table_a6b_numerator = abs(loss_scale - transport_field)
    table_a6b_denominator = max(abs(loss_scale), abs(transport_field))
    table_a6b_relative = table_a6b_numerator / table_a6b_denominator
    table_a6b_tolerance = float(manifest['tolerances']['G1']['rtol'])
    reference_density = power_density / (gas_density * constants.e * loss_scale)
    electron_density = multiplier * reference_density
    electron_temperature = (
        2.0 * row['mean_energy_eV'] * constants.e / (3.0 * constants.kB))
    denominator = pressure - electron_density * constants.kB * (
        electron_temperature - gas_temperature)
    if denominator <= 0.0:
        raise ValueError('requested seed is outside the two-temperature EOS domain')
    electron_fraction = electron_density * constants.kB * gas_temperature / denominator
    argon_fraction = 1.0 - 2.0 * electron_fraction
    if argon_fraction <= 0.0:
        raise ValueError('requested balanced seed leaves no neutral argon')
    heavy_moles = argon_fraction + electron_fraction
    reference_volume = heavy_moles * constants.R * gas_temperature / pressure
    reactor_power = (DEVELOPMENT_ABSORBED_POWER_W * reference_volume /
                     DEVELOPMENT_CHAMBER_VOLUME_M3)
    eos_volume = constants.R * (
        heavy_moles * gas_temperature + electron_fraction * electron_temperature) / pressure
    exact_density = electron_fraction * constants.Na / eos_volume

    old_fraction = 1.0e-6
    old_argon = 1.0 - 2.0 * old_fraction
    old_heavy = old_argon + old_fraction
    old_volume = constants.R * (
        old_heavy * gas_temperature + old_fraction * electron_temperature) / pressure
    old_neutral_density = old_argon * constants.Na / old_volume
    old_field_power = (old_fraction * constants.Na * old_neutral_density *
                       constants.e * loss_scale)
    return {
        'policy': 'power-consistent quasineutral pure-Ar development seed',
        'provenance': DEVELOPMENT_SEED_RULING,
        'artifact_sha256': manifest['artifact_sha256'],
        'branch': branch,
        'prescribed_absorbed_power_W': DEVELOPMENT_ABSORBED_POWER_W,
        'chamber_volume_m3': DEVELOPMENT_CHAMBER_VOLUME_M3,
        'power_density_W_m^-3': power_density,
        'pressure_Pa': pressure,
        'gas_temperature_K': gas_temperature,
        'gas_density_m^-3': gas_density,
        'initial_reduced_field_Td': DEVELOPMENT_INITIAL_FIELD_TD,
        'table_u': row['u'],
        'table_bracketing_row_indices': row['bracketing_row_indices'],
        'table_bracketing_EN_Td': row['bracketing_EN_Td'],
        'table_field_power_coefficient_eV_m3_s^-1': loss_scale,
        'table_transport_field_coefficient_eV_m3_s^-1': transport_field,
        'table_A6b_numerator': table_a6b_numerator,
        'table_A6b_denominator': table_a6b_denominator,
        'table_A6b_relative': table_a6b_relative,
        'table_A6b_tolerance': table_a6b_tolerance,
        'table_A6b_outcome': ('PASS' if table_a6b_relative <= table_a6b_tolerance
                              else 'FAIL'),
        'table_mean_energy_eV': row['mean_energy_eV'],
        'table_effective_electron_temperature_K': electron_temperature,
        'reference_electron_density_m^-3': reference_density,
        'seed_multiplier': multiplier,
        'electron_density_m^-3': exact_density,
        'electron_mole_fraction': electron_fraction,
        'argon_ion_mole_fraction': electron_fraction,
        'argon_mole_fraction': argon_fraction,
        'initial_charge_number_density_m^-3': 0.0,
        'reactor_reference_volume_m3': reference_volume,
        'reactor_power_scaling': 'P_reactor = P_chamber * V_ref / V_chamber',
        'reactor_absorbed_power_W': reactor_power,
        'over_ionised_fixture_mole_fraction': old_fraction,
        'over_ionised_field_power_W': old_field_power,
    }


def _nested_value(record, dotted):
    value = record
    for part in dotted.split('.'):
        value = value[part]
    return float(value)


def compare_development_endpoints(arms):
    """Compare endpoints, and report convergence and transient differences."""
    if not arms:
        raise ValueError('at least one sweep arm is required')
    reference = next((arm for arm in arms if arm.get('arm') == '1x'), arms[0])
    metrics = {}
    all_passed = True
    for name, tolerance in DEVELOPMENT_SWEEP_TOLERANCES.items():
        expected = _nested_value(reference, name)
        values = {arm['arm']: _nested_value(arm, name) for arm in arms}
        checks = {}
        for arm, value in values.items():
            limit = tolerance['atol'] + tolerance['rtol'] * max(abs(expected), abs(value))
            checks[arm] = {
                'absolute_difference': abs(value - expected),
                'limit': limit,
                'passed': abs(value - expected) <= limit,
            }
            all_passed = all_passed and checks[arm]['passed']
        metrics[name] = {'reference': expected, 'values': values, 'checks': checks}
    convergence = {
        arm['arm']: {
            'value_s': float(arm['convergence_time_s']),
            'difference_from_reference_s': float(
                arm['convergence_time_s'] - reference['convergence_time_s']),
            'ratio_to_reference': float(
                arm['convergence_time_s'] / reference['convergence_time_s']),
        }
        for arm in arms
    }
    transient_names = sorted({
        name for arm in arms for name in arm['transient_extrema']})
    transients = {
        name: {
            arm['arm']: copy.deepcopy(arm['transient_extrema'][name])
            for arm in arms
        }
        for name in transient_names
    }
    return {
        'tolerances': copy.deepcopy(DEVELOPMENT_SWEEP_TOLERANCES),
        'reference_arm': reference['arm'],
        'metrics': metrics,
        'all_endpoints_agree': all_passed,
        'convergence_time_comparison': {
            'agreement_required': False,
            'reason': 'initial-condition transients need not agree',
            'arms': convergence,
        },
        'transient_extrema_comparison': {
            'agreement_required': False,
            'reason': 'initial-condition transients need not agree',
            'metrics': transients,
        },
    }


def external_input_paths(spec):
    """Return every external file pinned by a loaded generation spec."""
    paths = [spec['binary']['path'], spec['cmake_cache']['path'], spec['channel_map']['path']]
    paths.extend(item['path'] for item in spec['input_files'].values())
    paths.extend(spec['shared_objects'])
    return tuple(Path(path) for path in paths)


def load_real_spec(spec_path):
    """Load ``spec_path`` and raise with the exact missing external inputs."""
    path = Path(spec_path)
    if not path.is_file():
        raise FileNotFoundError('EEDF_REAL_SPEC does not exist: {0}'.format(path))
    spec = load_spec(path)
    missing = [str(item) for item in external_input_paths(spec) if not item.is_file()]
    if missing:
        raise FileNotFoundError('missing EEDF external inputs: {0}'.format(', '.join(missing)))
    return spec


def generate_real_artifact(spec_path, work_directory, command_line):
    """Generate the full, unmodified-grid artifact under ``work_directory``."""
    spec = copy.deepcopy(load_real_spec(spec_path))
    work_directory = Path(work_directory)
    spec['output_root'] = str(work_directory / 'tables')
    spec['scratch_root'] = str(work_directory / 'loki')
    artifact = Path(generate(spec, command_line=command_line))
    manifest = json.loads((artifact / 'manifest.json').read_text())
    digest = file_hash(artifact / 'table.h5')
    if digest != manifest['artifact_sha256']:
        raise RuntimeError('generated table hash differs from its manifest')
    return artifact, manifest


def reactor_readiness(manifest):
    """Name real-artifact blockers that prevent an honest reactor run."""
    blockers = list(unqualified_artifact_blockers(manifest))
    mapped = [row for row in manifest['channel_map'] if row.get('reaction') is not None]
    if not mapped:
        classes = {}
        for row in manifest['channel_map']:
            label = row.get('classification')
            classes[label] = classes.get(label, 0) + 1
        blockers.append(
            'artifact has no reaction-owned channels; channel classifications are {0!r}'.format(classes))
    return blockers


def pure_argon_deck(artifact, manifest, seed_multiplier=1.0,
                    restart_state=None, over_ionised=False):
    """Render the dispatch operating point with no machine path in repo data."""
    branch = manifest['branches'][0]
    if restart_state is not None and over_ionised:
        raise ValueError('restart and over-ionised fixture are mutually exclusive')
    if restart_state is not None:
        initial = dict(restart_state['species_amounts_mol'])
        field_td = float(math.exp(restart_state['u']))
    elif over_ionised:
        initial = {'Ar': 0.999998, 'Arp': 1.0e-6, 'e-': 1.0e-6, 'Ars': 0.0}
        field_td = DEVELOPMENT_INITIAL_FIELD_TD
    else:
        initialization = power_consistent_initialization(
            artifact, manifest, multiplier=seed_multiplier)
        initial = {
            'Ar': initialization['argon_mole_fraction'],
            'Arp': initialization['argon_ion_mole_fraction'],
            'e-': initialization['electron_mole_fraction'],
            'Ars': 0.0,
        }
        field_td = DEVELOPMENT_INITIAL_FIELD_TD
    if set(initial) != {'Ar', 'Arp', 'e-', 'Ars'}:
        raise ValueError('initial state must contain exactly Ar, Arp, e-, and Ars')
    declaration = {
        'provider': 'loki-table',
        'table': (str(Path(artifact) / 'table.h5'), manifest['artifact_sha256']),
        'branch': branch,
        'initialReducedField': (field_td, 'Td'),
        'empiricalLaws': {
            'PlasmaArgon:88': {
                'evaluate_at': 'Te_eff',
                'class': 'B',
                'basis': 'I-318: Ar(4s) ionisation cross section absent from LoKI set',
                'sensitivity': True,
            },
            'PlasmaArgon:91': {
                'evaluate_at': 'Te_eff',
                'class': 'C',
                'basis': 'I-318: retained metastable-mixing on/off sensitivity',
                'sensitivity': True,
            },
        },
    }
    deck = """database(
    thermoLibraries=['primaryThermoLibrary', 'PlasmaThermo', 'electrocatThermo',
                     'PlasmaExcitedNeutralThermo'],
    reactionLibraries=['PlasmaArgon', 'PlasmaRadiativeRecombination'],
    seedMechanisms=[], kineticsDepositories=['training'],
    kineticsFamilies=['Plasma_Electron_Attachment'], kineticsEstimator='rate rules',
)

species(label='e-', reactive=True, structure=adjacencyList("1 e u1 p0 c-1"))
species(label='Ar', reactive=True, structure=adjacencyList("1 Ar u0 p4 c0"))
species(label='Arp', reactive=True, structure=adjacencyList("multiplicity 2\\n1 Ar u1 p3 c+1"))
species(label='Ars', reactive=True, structure=adjacencyList("multiplicity 3\\n1 Ar u2 p3 c0"))

plasmaReactor(
    temperature=(298.15, 'K'),
    # Exact pressure fingerprint used by the pinned spec: 5 torr = 666.61184 Pa.
    pressure=(666.61184, 'Pa'),
    initialMoleFractions=__INITIAL_STATE__,
    chamberGeometry={'diffusionLength': (20.3, 'mm')},
    ionReducedMobility=(1.535e-4, 'm^2/(V*s)'), wallRecycling=1.0,
    wallNeutralizationProducts={'Arp': 'Ar'},
    wallNeutralDiffusion={'Ars': {'product': 'Ar',
                                  'diffusivity': (47.0, 'cm^2*torr/s')}},
    ionisationSource=(6.6e4, 'm^-3/s'),
    electronKinetics=__ELECTRON_KINETICS__,
    electronEnergyBalance={'absorbedPower': (0.5, 'W'),
                           'sheath': 'floatingWall',
                           'chamberVolume': (2.356194490192345e-3, 'm^3')},
    terminationSteadyState=1e-6, terminationTime=(600.0, 's'),
)

simulator(atol=1e-16, rtol=1e-8)
model(toleranceKeepInEdge=0.0, toleranceMoveToCore=0.5,
      toleranceInterruptSimulation=1e8, maximumEdgeSpecies=200)
options(units='si', generateOutputHTML=False, generatePlots=False,
        saveEdgeSpecies=False, saveSimulationProfiles=True,
        generateChemkin=False, generateRMSYAML=False,
        generateCanteraYAML1=False, generateCanteraYAML2=False)
"""
    return deck.replace('__ELECTRON_KINETICS__', repr(declaration)).replace(
        '__INITIAL_STATE__', repr(initial))


def development_argon_model(core_reactions):
    """Return runtime-only reactions for the development pure-Ar deck.

    The database objects are never edited.  Entry 89 and the radiative
    recombination library are intentionally absent from this development
    deck: its electron reactions are the two table-owned ground-state
    channels plus the explicitly declared entries 88 and 91.
    """
    selected = []
    seen = set()
    for source in core_reactions:
        entry = getattr(getattr(source, 'entry', None), 'index', None)
        if getattr(source, 'library', None) != 'PlasmaArgon':
            continue
        if entry not in DEVELOPMENT_RETAINED_ENTRIES:
            continue
        reaction = copy.copy(source)
        if entry in DEVELOPMENT_BINDINGS:
            reaction.kinetics = EEDFChannel(
                process=DEVELOPMENT_BINDINGS[entry],
                collision_set=DEVELOPMENT_COLLISION_SET,
                side='ine')
        selected.append(reaction)
        seen.add(entry)
    missing = sorted(DEVELOPMENT_RETAINED_ENTRIES - seen)
    if missing:
        raise RuntimeError(
            'development pure-Ar deck is missing PlasmaArgon entries {0!r}'.format(missing))
    return selected
