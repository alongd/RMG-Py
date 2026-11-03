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

"""Synthetic source/thermo inputs are generated in scratch, never physical data."""
import copy
import hashlib
import json
import os
from pathlib import Path
import sys
from types import SimpleNamespace
import zipfile

import pytest

from rmgpy.tools.eedf.n2_vibrational import (
    LadderError, MissingPairsError, N2_GRAPH, assess_truncation, generate,
    main, read_cross_sections, require_complete_pairs,
)


def pin(path):
    return {'path': str(path.resolve()), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'citation': 'Synthetic fixture only; not physical N2 data', 'provenance': 'sourced'}


def sections_text(upper, omit=None):
    text = []
    for i in range(upper + 1):
        for j in range(i, upper + 1):
            if (i, j) == omit:
                continue
            text.append(f'EXCITATION\nN2(X,v={i}) -> N2(X,v={j})\n'
                        f'{0.2 * (j-i):.12f}\nCOMMENT: Synthetic, not physical data\n'
                        '----------\n0.01 0\n10 1e-20\n----------\n')
    return '\n'.join(text)


def partial_sections_text(initial_levels, final_v, omit=None):
    text = []
    for i in initial_levels:
        for j in range(i, final_v + 1):
            if (i, j) == omit:
                continue
            text.append(f'EXCITATION\nN2(X,v={i}) -> N2(X,v={j})\n'
                        f'{0.2 * (j-i):.12f}\nCOMMENT: Synthetic, not physical data\n'
                        '----------\n0.01 0\n10 1e-20\n----------\n')
    return '\n'.join(text)


@pytest.fixture
def source(tmp_path):
    upper = 32
    spectroscopy = {'molecule': '14N2', 'electronic_state': 'X1Sigma_g+',
                    'energy_unit': 'eV', 'energy_zero': 'v=0',
                    'levels': [{'v': v, 'energy_eV': 0.2*v, 'statistical_weight': v+1,
                                'thermo_entry': f'N2v{v}'} for v in range(upper+1)]}
    level_file = tmp_path / 'levels.json'
    level_file.write_text(json.dumps(spectroscopy))
    xs = tmp_path / 'synthetic.txt'
    xs.write_text(sections_text(upper))
    cross_section_pin = pin(xs)
    cross_section_pin['source_id'] = 'synthetic'
    library = tmp_path / 'library.py'
    text = ["name='SyntheticN2'\nshortDesc='Synthetic only, not physical data'\n"]
    for v in range(upper + 1):
        adj = f'vibrationallevel {v}\n' + N2_GRAPH
        text.append(f'entry(index={v}, label="N2v{v}", molecule={adj!r}, '
                    'thermo=NASA(polynomials=[NASAPolynomial('
                    f'coeffs=[3.5,0,0,0,0,{v*1000},3], Tmin=(200,"K"), Tmax=(6000,"K"))], '
                    'Tmin=(200,"K"), Tmax=(6000,"K"), '
                    'Cp0=(29.100619163,"J/(mol*K)"), CpInf=(29.100619163,"J/(mol*K)")))\n')
    library.write_text('\n'.join(text))
    return {'schema': 'n2-vibrational-source-v1', 'level_source': pin(level_file),
            'cross_sections': [cross_section_pin], 'source_domain': {
                'initial_levels': list(range(upper + 1)), 'final_v': upper},
            'source_assignments': {str(v): 'synthetic' for v in range(upper + 1)},
            'thermo_library': pin(library),
            'energy_tolerance_eV': 1e-9,
            'D6': {'set_by': 'owner', 'top_fraction_limit': 0.02,
                   'top_leak_share_limit': 0.03, 'K_step': 5}}


def decoded(source, K=10, production=True):
    return json.loads(generate(source, K, production)['manifest.json'])


def test_default_identities_and_v0_lookup(source, tmp_path, monkeypatch):
    import rmgpy.data.rmg as data_module
    from rmgpy.data.thermo import ThermoDatabase, ThermoLibrary
    from rmgpy.rmg.main import RMG

    monkeypatch.setattr(data_module, 'database', None)
    artifacts = generate(source)
    for name, data in artifacts.items():
        (tmp_path / name).write_bytes(data)
    job = RMG()
    job.load_input(str(tmp_path / 'input.py'))
    assert len(job.initial_species) == 11
    assert job.initial_species[0].molecule[0].vibrational_level == -1
    lib = ThermoLibrary(label='N2Vibrational')
    db = ThermoDatabase()
    lib.load(str(tmp_path / 'thermo.py'), db.local_context, {})
    db.libraries = {lib.label: lib}
    db.library_order = [lib.label]
    monkeypatch.setattr(data_module, 'database', SimpleNamespace(thermo=db, solvation=None))
    ground = job.initial_species[0]
    assert db.get_thermo_data(ground).label == 'N2v0'
    for v, species in enumerate(job.initial_species):
        # Probe fixed-v library lookup directly: model.generate_thermo also
        # requests an energy-transfer model, outside this generator's scope.
        thermo = db.get_thermo_data(species)
        assert thermo.label == f'N2v{v}'
    manifest = json.loads(artifacts['manifest.json'])
    from rmgpy.molecule.molecule import Molecule
    molecules = [Molecule().from_adjacency_list(s['adjacency_list']) for s in manifest['states']]
    assert [m.vibrational_level for m in molecules] == list(range(11))
    assert all(not a.is_isomorphic(b) for i,a in enumerate(molecules) for b in molecules[i+1:])
    assert manifest['states'][-1]['v'] == 10
    assert len(manifest['states']) == 11  # no extra lump species


@pytest.mark.parametrize('placeholder', ['0', '{}'])
def test_non_model_retained_thermo_refuses(source, placeholder):
    path = Path(source['thermo_library']['path'])
    lines = path.read_text().splitlines()
    lines = [line[:line.index('thermo=')] + f'thermo={placeholder})'
             if 'label="N2v1",' in line else line for line in lines]
    path.write_text('\n'.join(lines) + '\n')
    source['thermo_library'] = pin(path)
    with pytest.raises(LadderError, match=r'N2\(X,v=1\): N2v1'):
        generate(source)


def test_retained_thermo_uses_exact_library_lookup(source, monkeypatch):
    from rmgpy.data.thermo import ThermoDatabase

    lookup = ThermoDatabase.get_thermo_data_from_library
    seen = []

    def exact_lookup(database, species, library):
        match = lookup(database, species, library)
        seen.append((species.molecule[0].vibrational_level, match[2].label))
        return match

    def forbid_estimation(*args, **kwargs):
        pytest.fail('ladder validation must use only the exact library lookup')

    monkeypatch.setattr(ThermoDatabase, 'get_thermo_data_from_library', exact_lookup)
    monkeypatch.setattr(ThermoDatabase, 'get_thermo_data', forbid_estimation)
    generate(source)
    assert seen == [(v, f'N2v{v}') for v in range(11)]


@pytest.mark.parametrize('K', [10, 30])
def test_generation_and_regeneration_byte_determinism(source, K):
    first = generate(source, K)
    generate(source, 20)  # an intervening regeneration cannot change source state
    assert generate(copy.deepcopy(source), K) == first
    output = json.loads(first['manifest.json'])
    assert output['K'] == K and len(output['states']) == K+1
    assert output['D6']['top_v'] == K
    assert sum(c['direction'] == 'superelastic' for c in output['channels']) == K*(K+1)//2
    assert all(c['initial_v'] <= K and c['final_v'] <= K for c in output['channels'])
    assert output['production_manifest']
    assert 'qualified EEDF hierarchical N2(X,v) selectors' in output['runtime_requirements']


def test_complete_census_and_every_supplied_transition_audited(source):
    output = decoded(source)
    sections = read_cross_sections(source['cross_sections'][0])
    require_complete_pairs(sections, 32)
    audit = output['transition_audit']
    assert len(audit) == len(sections) == 33*34//2
    assert len({r['pair_id'] for r in audit}) == len(audit)
    assert sum(r['disposition'] == 'retained' for r in audit) == 55
    assert all(r['disposition'] == 'unresolved_initial_above_K'
               for r in audit if r['initial_v'] > 10)
    forward = [c for c in output['channels'] if c['direction'] == 'forward']
    expected = {(r['initial_v'], r['final_v']) for r in sections if r['initial_v'] <= 10}
    assert {(c['initial_v'], c['source_final_v']) for c in forward} == expected
    assert len(forward) == len(expected)
    assert all(c['final_v'] == 10 for c in forward if c['source_final_v'] > 10)
    # Initial top channels remain explicit top->top collisions with energy ledger.
    assert any(c['initial_v'] == c['final_v'] == 10 and c['source_final_v'] == 32 for c in forward)


def test_initial_coverage_is_separate_from_final_domain(tmp_path):
    initials = list(range(30))
    final_v = 58
    path = tmp_path / 'partial-domain.txt'
    path.write_text(partial_sections_text(initials, final_v))
    sections = read_cross_sections(pin(path))
    assert len(sections) == sum(final_v - i + 1 for i in initials)
    require_complete_pairs(sections, initial_levels=initials, final_v=final_v)
    with pytest.raises(MissingPairsError, match=r'\(29,58\)'):
        require_complete_pairs(sections[:-1], initial_levels=initials, final_v=final_v)
    with pytest.raises(LadderError, match='outside declared domain'):
        require_complete_pairs(sections, initial_levels=initials, final_v=57)


def test_archive_member_uses_hash_of_original_snapshot(tmp_path):
    text = sections_text(2)
    path = tmp_path / 'source.zip'
    with zipfile.ZipFile(path, 'w') as archive:
        archive.writestr('Cross section.txt', text)
    source_pin = pin(path)
    source_pin['member'] = 'Cross section.txt'
    assert len(read_cross_sections(source_pin)) == 6
    path.write_bytes(path.read_bytes() + b'changed')
    with pytest.raises(LadderError, match='source hash mismatch'):
        read_cross_sections(source_pin)


def test_generation_refuses_K_requiring_uncovered_initial_level(source):
    source['source_domain']['initial_levels'] = list(range(30))
    source['source_domain']['final_v'] = 32
    source['source_assignments'] = {str(v): 'synthetic' for v in range(30)}
    with pytest.raises(LadderError, match=r'requested K=30 requires uncovered initial v=30'):
        generate(source, 30)


@pytest.mark.parametrize(('field', 'value', 'message'), [
    ('source_domain', {'initial_levels': '0-29', 'final_v': 58}, 'initial_levels'),
    ('source_assignments', {'0': []}, 'keys and values must be strings'),
    ('D6', [], 'D6 policy must be an object'),
    ('cross_sections', {'source_id': 'bad'}, 'nonempty list of source pins'),
])
def test_malformed_source_domain_config_refuses_by_name(source, field, value, message):
    source[field] = value
    with pytest.raises(LadderError, match=message):
        generate(source)


def test_superelastic_only_for_two_retained_source_endpoints(source):
    output = decoded(source)
    reverse = [c for c in output['channels'] if c['direction'] == 'superelastic']
    forward = {c['pair_id']: c for c in output['channels'] if c['direction'] == 'forward'}
    assert len(reverse) == 55
    for partner in reverse:
        original = forward[partner['pair_id']]
        assert original['source_final_v'] <= 10
        assert partner['initial_v'] == original['final_v']
        assert partner['final_v'] == original['initial_v']
        assert partner['rate_binding'] == {'source_pair_index': original['source_pair_index'], 'column': 'k_sup'}
        assert original['rate_binding']['column'] == 'k_ine'
        assert not original['reversible'] and not partner['reversible']
        assert partner['detailed_balance']['source_pair_index'] == original['source_pair_index']
        assert partner['detailed_balance']['weight_ratio_initial_over_final'] == pytest.approx(
            (original['initial_v']+1)/(original['source_final_v']+1))
        assert partner['electron_energy_loss_eV'] == -original['electron_energy_loss_eV']
    assert not any(c['direction'] == 'superelastic' and c['source_final_v'] > 10 for c in output['channels'])


def test_particle_and_declared_energy_bookkeeping(source):
    output = decoded(source)
    for c in output['channels']:
        assert sum(c['reactants'].values()) == sum(c['products'].values()) == 2
        assert c['reactants']['electron'] == c['products']['electron'] == 1
        assert sum(n for s,n in c['reactants'].items() if s.startswith('N2')) == 1
        assert sum(n for s,n in c['products'].items() if s.startswith('N2')) == 1
        assert c['electron_energy_loss_eV'] == pytest.approx(
            c['retained_energy_gain_eV'] + c['top_leak_eV'] + c['source_energy_residual_eV'])
    ledger = {r['pair_id']: r for r in output['top_leak_ledger']}
    pumped = [c for c in output['channels'] if c['direction'] == 'forward' and c['source_final_v'] > 10]
    assert len(ledger) == len(pumped) == 11*22
    for c in pumped:
        entry = ledger[c['pair_id']]
        assert entry['energy_eV_per_event'] == pytest.approx(0.2*(c['source_final_v']-10))
        assert entry['molecular_gain_eV'] == -entry['electron_loss_eV'] == c['top_leak_eV']
        assert entry['gas_heating_eV'] == 0
        assert entry['provenance']['provenance'] == 'sourced'
        assert entry['name'] == 'top-leak:' + c['pair_id']


def test_incomplete_synthetic_refuses_by_named_pair_before_thermo(source):
    xs = Path(source['cross_sections'][0]['path'])
    xs.write_text(sections_text(32, omit=(4,9)))
    source['cross_sections'][0] = pin(xs)
    source['cross_sections'][0]['source_id'] = 'synthetic'
    source['thermo_library']['path'] = '/nonexistent/not-reached'
    with pytest.raises(MissingPairsError, match=r'\(4,9\)') as caught:
        generate(source)
    assert caught.value.pairs == ((4,9),)


def test_real_laporta_incomplete_coverage(source):
    path = os.environ.get('N2_LAPORTA_BATCH1')
    if path is None:
        pytest.skip('set N2_LAPORTA_BATCH1 for the local-only real source refusal witness')
    real_pin = pin(Path(path))
    real_pin['source_id'] = 'laporta_batch1'
    real_pin['citation'] = 'Laporta, Celiberto & Wadehra, PSST 21 055018 (2012); LXCat Laporta'
    sections = read_cross_sections(real_pin)
    assert len(sections) == 1000
    with pytest.raises(MissingPairsError) as caught:
        require_complete_pairs(sections, 58)
    missing = caught.value.pairs
    assert (3,39) in missing and (4,5) in missing and (30,31) in missing
    assert len([p for p in missing if p[1] <= 10]) == 21
    assert len([p for p in missing if p[1] <= 30]) == 141
    assert '(4,5)' in str(caught.value)
    # Synthetic spectroscopy placeholders only: refusal precedes their use.
    level_path = Path(source['level_source']['path'])
    data = json.loads(level_path.read_bytes())
    data['levels'].extend({'v': v, 'energy_eV': 0.2*v, 'statistical_weight': v+1,
                           'thermo_entry': f'N2v{v}'} for v in range(33,59))
    level_path.write_text(json.dumps(data))
    source['level_source'] = pin(level_path)
    source['cross_sections'] = [real_pin]
    source['source_domain'] = {'initial_levels': list(range(59)), 'final_v': 58}
    source['source_assignments'] = {str(v): 'laporta_batch1' for v in range(59)}
    with pytest.raises(MissingPairsError) as refused:
        generate(source)
    assert refused.value.pairs == missing


@pytest.mark.parametrize('what', ['hash', 'estimated', 'citation', 'thermo_identity', 'threshold', 'identity_override'])
def test_source_and_thermo_refusals(source, what):
    if what == 'hash':
        source['level_source']['sha256'] = '0'*64
    elif what == 'estimated':
        source['level_source']['provenance'] = 'estimated'
    elif what == 'citation':
        source['level_source']['citation'] = ''
    elif what == 'thermo_identity':
        p = Path(source['thermo_library']['path'])
        p.write_text(p.read_text().replace('vibrationallevel 0', 'vibrationallevel 99'))
        source['thermo_library'] = pin(p)
    else:
        p = Path(source['level_source']['path'])
        data = json.loads(p.read_bytes())
        if what == 'identity_override':
            data['levels'][1]['adjacency_list'] = 'vibrationallevel 0\n' + N2_GRAPH
        else:
            data['levels'][1]['energy_eV'] += 0.01
        p.write_text(json.dumps(data))
        source['level_source'] = pin(p)
    with pytest.raises(LadderError):
        generate(source)


def test_duplicate_or_outside_source_domain_refuses(source):
    sections = read_cross_sections(source['cross_sections'][0])
    with pytest.raises(LadderError, match='duplicate'):
        require_complete_pairs(sections + [sections[0]], 32)
    with pytest.raises(LadderError, match='outside'):
        require_complete_pairs(sections, 31)


def test_owner_threshold_placeholder_refuses(source):
    source['D6']['top_fraction_limit'] = None
    with pytest.raises(LadderError, match='top_fraction_limit'):
        generate(source)
    draft = decoded(source, production=False)
    assert not draft['production_manifest']
    assert draft['D6']['top_fraction_limit'] is None
    with pytest.raises(LadderError, match='top_fraction_limit'):
        assess_truncation(10, [1]*11, 1, 0, draft['D6'])


def test_d6_accept_raise_and_ceiling(source):
    policy = source['D6']
    assert assess_truncation(10, [1]+[0]*10, 0, 0, policy)['action'] == 'accept'
    assert assess_truncation(10, [0]*10+[1], 1, 0, policy)['K'] == 15
    decision = assess_truncation(28, [1]+[0]*28, 1, 0.04, policy)
    assert decision['K'] == 30 and decision['exceeded'] == ['top_leak_share_limit']
    for K in [10,15,20,25]:
        assert assess_truncation(K, [0]*K+[1], 1, 0, policy)['K'] == K+5
    with pytest.raises(LadderError, match='D6 ceiling K=30 exceeded: top_fraction_limit'):
        assess_truncation(30, [0]*30+[1], 1, 0, policy)
    with pytest.raises(LadderError, match='top_leak_share_limit'):
        assess_truncation(30, [1]+[0]*30, 1, 0.04, policy)
    assert assess_truncation(30, [1]+[0]*30, 1, 0, policy)['action'] == 'accept'
    assert assess_truncation(10, [98]+[0]*9+[2], 100, 3, policy)['action'] == 'accept'


@pytest.mark.parametrize('populations,forward,leak', [([0]*11,1,0), ([-1]+[0]*10,1,0),
    ([float('nan')]+[0]*10,1,0), ([1]*10,1,0), ([1]*11,0,1), ([1]*11,1,-1)])
def test_d6_invalid_runtime_values_refuse(source, populations, forward, leak):
    with pytest.raises(LadderError):
        assess_truncation(10, populations, forward, leak, source['D6'])


@pytest.mark.parametrize('K', [0,31,True,10.5])
def test_invalid_K_refuses(source, K):
    with pytest.raises(LadderError, match='K must be an integer'):
        generate(source, K)


def test_cli_emits_without_payload_and_refuses_overwrite(source, tmp_path, monkeypatch):
    source_path = tmp_path / 'source.json'
    source_path.write_text(json.dumps(source))
    output = tmp_path / 'output'
    monkeypatch.setattr(sys, 'argv', ['n2_vibrational', str(source_path), str(output)])
    main()
    assert {p.name for p in output.iterdir()} == {'manifest.json','input.py','thermo.py'}
    assert (output / 'manifest.json').read_bytes() == generate(source)['manifest.json']
    document = (output / 'manifest.json').read_text()
    assert 'sigma_m2' not in document and 'energy_eV": [' not in document
    assert 'tail_populations' not in document and 'Treanor' not in document
    with pytest.raises(FileExistsError):
        main()
