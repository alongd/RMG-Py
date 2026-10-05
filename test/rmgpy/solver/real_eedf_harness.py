"""Real LoKI-B fixture and pure-Ar EEDF run harness support.

This module owns no scientific inputs.  It reads the spec named by
``EEDF_REAL_SPEC``, verifies that every external file exists, and lets the
existing generator qualify the artifact.  It never edits acceptance,
fingerprints, branch metadata, classifications, or manifests.
"""

import copy
import json
from pathlib import Path

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


def pure_argon_deck(artifact, manifest):
    """Render the dispatch operating point with no machine path in repo data."""
    branch = manifest['branches'][0]
    declaration = {
        'provider': 'loki-table',
        'table': (str(Path(artifact) / 'table.h5'), manifest['artifact_sha256']),
        'branch': branch,
        'initialReducedField': (17.0, 'Td'),
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
    initialMoleFractions={'Ar': 0.999998, 'Arp': 1.0e-6,
                          'e-': 1.0e-6, 'Ars': 0.0},
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
    return deck.replace('__ELECTRON_KINETICS__', repr(declaration))


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
