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

import json

import numpy as np
import pytest

from rmgpy.tools.eedf.report import ReportRefusal, _main, power_table, render_markdown


def channels():
    return [
        {'description': 'e + A(X) -> e + A(X*), Excitation',
         'kind': 'excitation', 'classification': 'B', 'reaction': None,
         'flux_group': 'e + A(X) -> e + A(X*), Excitation'},
        {'description': 'e + B(v0) -> e + B(v1), Vibrational',
         'kind': 'vibrational', 'classification': 'B', 'reaction': None,
         'flux_group': 'e + B(v0) -> e + B(v1), Vibrational'},
        {'description': 'e + C(j0) -> e + C(j1), Rotational',
         'kind': 'rotational', 'classification': 'B', 'reaction': None,
         'flux_group': 'e + C(j0) -> e + C(j1), Rotational'},
        {'description': 'e + D -> e + e + D+, Ionization',
         'kind': 'ionization', 'classification': 'B', 'reaction': None,
         'flux_group': 'e + D -> e + e + D+, Ionization'},
        {'description': 'e + F -> F-, Attachment',
         'kind': 'attachment', 'classification': 'B', 'reaction': None,
         'flux_group': 'e + F -> F-, Attachment'},
    ]


def budget(**overrides):
    declared = channels()
    powers = [20., 5., 3., 0., 0.]
    value = {
        'P_abs': 100., 'Q_inelastic': 28., 'Q_elastic': 10.,
        'Q_wall_electron': 1., 'Q_wall_ion': 1., 'Q_flow': 1., 'dU_dt': 59.,
        'closure': 0., 'Q_inelastic_channels': powers,
        'Q_superelastic_channels': [0., 0., 0., 0., 0.],
        'Q_energy_only': {
            channel['flux_group']: powers[index]
            for index, channel in enumerate(declared)
            if channel['classification'] == 'B'
        },
        'Q_heavy_particle_electron_by_reaction': [0.],
    }
    value.update(overrides)
    return value


def budget_for(declared, **overrides):
    value = budget(**overrides)
    if 'Q_inelastic' not in overrides:
        value['Q_inelastic'] = (
            sum(value['Q_inelastic_channels'])
            - sum(value['Q_superelastic_channels'])
            + sum(value['Q_heavy_particle_electron_by_reaction'])
        )
    value['Q_energy_only'] = {
        channel['flux_group']: (
            value['Q_inelastic_channels'][index]
            - value['Q_superelastic_channels'][index]
        )
        for index, channel in enumerate(declared)
        if channel['classification'] == 'B'
    }
    return value


def test_every_row_and_renderer():
    report = power_table(budget(), channels())
    names = {row['name'] for row in report.rows}
    assert {'absorbed electron power', 'elastic transfer', 'excitation', 'ionisation',
            'attachment/detachment', 'dissociation', 'superelastic return', 'wall loss',
            'heavy-particle electron transfer', 'flow loss', 'stored-energy change',
            'vibrational', 'rotational', 'numerical residual'} <= names
    assert 'a-posteriori qualification: not run' in render_markdown(report)
    assert 'Electron-energy disclosure:' in render_markdown(report)
    assert report.unit == 'W (reactor inventory)'
    assert 'Power — W (reactor inventory)' in render_markdown(report)


@pytest.mark.parametrize('power, expected', [(21.9, False), (22., True)])
def test_excitation_switch(power, expected):
    declared = channels()
    b = budget_for(declared, Q_inelastic_channels=[power, 0., 0., 0., 0.],
                   dU_dt=100. - 10. - 2. - 1. - power, closure=0.)
    rendered = render_markdown(power_table(b, declared))
    assert ('Excitation channels:' in rendered) is expected


def test_closure_mismatch_raises():
    with pytest.raises(ReportRefusal, match='arithmetic'):
        power_table(budget(closure=1.), channels())


def test_relative_arithmetic_tolerance_rejects_small_budget_five_percent_error():
    declared = channels()
    small = budget_for(
        declared,
        P_abs=1e-8,
        Q_inelastic=1.05e-8,
        Q_elastic=0.,
        Q_wall_electron=0.,
        Q_wall_ion=0.,
        Q_flow=0.,
        dU_dt=0.,
        closure=0.,
        Q_inelastic_channels=[1.05e-8, 0., 0., 0., 0.],
    )
    with pytest.raises(ReportRefusal, match='arithmetic'):
        power_table(small, declared)


def test_p1_fail():
    report = power_table(budget(closure=2., dU_dt=57.), channels(), p1_rtol=1e-3)
    rows = {row['name']: row['power_W'] for row in report.rows}
    assert report.p1 == 'FAIL'
    assert rows['numerical residual'] == 2.
    assert report.arithmetic_residual == 0.


def test_p1_uses_declared_absolute_tolerance_at_tiny_power():
    declared = channels()
    tiny = budget_for(
        declared,
        P_abs=1e-16,
        Q_inelastic_channels=[0., 0., 0., 0., 0.],
        Q_elastic=0.,
        Q_wall_electron=0.,
        Q_wall_ion=0.,
        Q_flow=0.,
        dU_dt=-4e-16,
        closure=5e-16,
    )
    assert power_table(tiny, declared, p1_rtol=0.).p1 == 'PASS'


@pytest.mark.parametrize(('absorbed', 'closure', 'stored'), [
    (-100., .1, -141.1),
    (0., 0., -41.),
])
def test_absorbed_power_must_be_strictly_positive(absorbed, closure, stored):
    with pytest.raises(ReportRefusal, match='P_abs must be positive'):
        power_table(budget(P_abs=absorbed, closure=closure, dU_dt=stored), channels())


def test_unclassifiable_refuses():
    bad = channels()
    bad[0] = {'description': 'mystery', 'kind': 'mystery',
              'classification': 'A', 'reaction': None}
    with pytest.raises(ReportRefusal, match='unclassifiable channel: mystery'):
        power_table(budget(), bad)


def test_double_counted_reaction_refuses():
    bad = channels()
    bad[0]['classification'] = 'A'
    bad[0]['reaction'] = {'index': 0, 'core_reaction_index': 0,
                          'repr': ("EEDFChannel(process='e + A(X) -> e + A(X*), Excitation', "
                                   "collision_set='fixture', side='ine')")}
    with pytest.raises(ReportRefusal, match='double-counted core reaction'):
        power_table(budget_for(
            bad, Q_heavy_particle_electron_by_reaction=[1.]), bad)


def test_qualification_missing_is_never_pass():
    assert 'not run' in power_table(budget(), channels()).qualification


def test_qualification_string_is_unverified_and_never_prints_pass():
    rendered = render_markdown(power_table(budget(), channels(), qualification='PASS'))
    assert 'a-posteriori qualification: PASS' not in rendered
    assert 'a-posteriori qualification: unverified record' in rendered


def test_legacy_budget_refuses():
    with pytest.raises(ReportRefusal, match='legacy'):
        power_table({'P_abs': 1.}, channels())


@pytest.mark.parametrize('bad', [
    {20.: None, 5.: None, 3.: None, 0.: None},
    (20., 5., 3., 0., 0.),
    np.asarray([20., 5., 3., 0., 0.]),
    [20, 5., 3., 0., 0.],
    [True, 5., 3., 0., 0.],
    ['20', 5., 3., 0., 0.],
    [None, 5., 3., 0., 0.],
    [float('nan'), 5., 3., 0., 0.],
    [float('inf'), 5., 3., 0., 0.],
])
def test_numeric_vectors_are_lists_of_finite_floats_only(bad):
    with pytest.raises(ReportRefusal, match='list of finite floats'):
        power_table(budget(Q_inelastic_channels=bad), channels())


def test_classification_b_ionization_lands_in_ionization():
    b_channels = channels()
    b_channels[-1] = {'description': 'e + X -> e + e + X+, Ionization', 'kind': 'ionization',
                      'classification': 'B', 'flux_group': 'ionisation B',
                      'reaction': None}
    report = power_table(budget_for(
        b_channels, Q_inelastic_channels=[20., 5., 3., 0., 7.], dU_dt=52.),
        b_channels)
    rows = {row['name']: row for row in report.rows if row.get('included_in_balance', True)}
    assert rows['ionisation']['power_W'] == 7.
    assert {'group': 'ionisation B', 'row': 'ionization', 'power_W': 7.} \
        in report.energy_only_groups


def test_top_level_partition_excludes_energy_only_breakdown():
    b_channels = channels()
    b_channels[0] = {'description': 'e + X -> e + X*, Excitation', 'kind': 'excitation',
                     'classification': 'B', 'flux_group': 'excitation B',
                     'reaction': None}
    report = power_table(budget_for(b_channels), b_channels)
    top_level = sum(row['power_W'] for row in report.rows[1:]
                    if row.get('included_in_balance', True)
                    and row['name'] != 'stored-energy change')
    assert budget()['P_abs'] == pytest.approx(top_level + budget()['dU_dt'] + report.closure)
    assert any(not row.get('included_in_balance', True) and row.get('level') == 1
               for row in report.rows)
    assert '&nbsp;&nbsp;↳ energy-only: excitation B' in render_markdown(report)


def test_energy_only_group_without_declared_row_refuses():
    with pytest.raises(ReportRefusal, match='Q_energy_only inventory'):
        power_table(budget(Q_energy_only={'missing': 1.}), channels())


def test_reaction_identity_routes_dissociation_and_detachment():
    mapped = channels()
    mapped[0]['classification'] = 'A'
    mapped[0]['description'] = 'e + N2(X) -> e + N(4S) + N(4S), Excitation'
    mapped[0]['reaction'] = {'index': 99, 'core_reaction_index': 0,
                             'repr': ("EEDFChannel(process='e + N2(X) -> e + N(4S) + "
                                      "N(4S), Excitation', collision_set='fixture', side='ine')")}
    mapped[4]['description'] = 'e + A-(1) -> A(1) + e + e, Attachment'
    mapped[4]['classification'] = 'A'
    mapped[4]['reaction'] = {'index': 98, 'core_reaction_index': 0,
                             'repr': ("EEDFChannel(process='e + A-(1) -> A(1) + e + e, "
                                      "Attachment', collision_set='fixture', side='ine')")}
    report = power_table(budget_for(
        mapped, Q_inelastic_channels=[20., 5., 3., 0., 7.], dU_dt=52.), mapped)
    rows = {row['name']: row for row in report.rows}
    assert rows['dissociation']['power_W'] == 20.
    assert rows['attachment/detachment']['power_W'] == 7.


def test_reaction_repr_is_not_parsed_when_structured_fields_resolve_channel():
    mapped = channels()
    mapped[0]['classification'] = 'A'
    mapped[0]['reaction'] = {'index': 1, 'repr': 'e + N2(X) -> e + N(4S) + N(4S)'}
    assert power_table(budget_for(mapped), mapped).p1 == 'PASS'


def test_energy_only_dissociation_uses_declared_process_identity():
    declared = channels()
    declared[0] = {
        'description': 'e + N2(X) -> e + N(4S) + N(4S), Excitation',
        'kind': 'excitation', 'classification': 'B',
        'flux_group': 'nitrogen dissociation', 'reaction': None,
    }
    report = power_table(budget_for(declared), declared)
    rows = {row['name']: row['power_W'] for row in report.rows}
    assert rows['excitation'] == 0.
    assert rows['dissociation'] == 20.
    assert {
        'group': 'nitrogen dissociation', 'row': 'dissociation', 'power_W': 20.,
    } in report.energy_only_groups


def test_excitation_breakdown_uses_final_row_assignment():
    declared = channels()
    declared[0]['description'] = 'e + N2(X) -> e + N(4S) + N(4S), Excitation'
    declared[1] = {
        'description': 'e + N2(X) -> e + N2(A), Excitation',
        'kind': 'excitation', 'classification': 'B', 'reaction': None,
        'flux_group': 'e + N2(X) -> e + N2(A), Excitation',
    }
    report = power_table(budget_for(
        declared, Q_inelastic=50., Q_inelastic_channels=[20., 30., 0., 0., 0.],
        dU_dt=37.), declared)
    assert report.excitation_channels == [{
        'channel': 'e + N2(X) -> e + N2(A), Excitation', 'power_W': 30.,
    }]


def test_unresolved_declared_process_refuses_instead_of_using_kind():
    declared = channels()
    declared[0] = {
        'description': 'unresolved collision', 'kind': 'excitation',
        'classification': 'A', 'core_reaction_index': 0,
        'reaction': {
            'index': 1,
            'repr': "EEDFChannel(process='unresolved collision', collision_set='fixture', side='ine')",
        },
    }
    with pytest.raises(ReportRefusal, match='unclassifiable channel: unresolved collision'):
        power_table(budget(), declared)


def test_superelastic_and_unmapped_heavy_powers_are_persisted_rows():
    mapped = channels()
    report = power_table(budget_for(
        mapped, Q_inelastic_channels=[20., 0., 0., 0., 0.],
        Q_superelastic_channels=[0., 5., 0., 0., 0.],
        Q_heavy_particle_electron_by_reaction=[0., 5.], dU_dt=67.), mapped)
    rows = {row['name']: row for row in report.rows}
    assert rows['superelastic return']['power_W'] == -5.
    assert rows['heavy-particle electron transfer']['power_W'] == 5.
    assert report.arithmetic_residual == 0.


def test_energy_only_groups_are_cross_checked_net_of_superelastic_power():
    declared = channels()
    net = budget(
        Q_inelastic=23.,
        Q_superelastic_channels=[5., 0., 0., 0., 0.],
        dU_dt=64.,
    )
    group = declared[0]['flux_group']
    net['Q_energy_only'][group] = 15.
    assert power_table(net, declared).p1 == 'PASS'

    for mixed in (20., 10.):
        invalid = json.loads(json.dumps(net))
        invalid['Q_energy_only'][group] = mixed
        with pytest.raises(ReportRefusal, match='power disagrees'):
            power_table(invalid, declared)


def test_double_count_uses_core_position_not_source_index():
    mapped = channels()
    mapped[0]['classification'] = 'A'
    mapped[0]['reaction'] = {
        'index': 100, 'core_reaction_index': 1,
        'repr': "EEDFChannel(process='e + A(X) -> e + A(X*), Excitation', collision_set='fixture', side='ine')"}
    with pytest.raises(ReportRefusal, match='double-counted core reaction: 1'):
        power_table(budget_for(
            mapped, Q_heavy_particle_electron_by_reaction=[0., 5., 0.]), mapped)


def test_double_count_without_core_position_refuses():
    mapped = channels()
    mapped[0]['classification'] = 'A'
    mapped[0]['reaction'] = {
        'index': 100,
        'repr': "EEDFChannel(process='e + A(X) -> e + A(X*), Excitation', collision_set='fixture', side='ine')"}
    with pytest.raises(ReportRefusal, match='core reaction position unavailable'):
        power_table(budget_for(
            mapped, Q_heavy_particle_electron_by_reaction=[0., 5., 0.]), mapped)


@pytest.mark.parametrize('reaction', [
    {'index': None, 'core_reaction_index': 0,
     'repr': "EEDFChannel(process='e + A(X) -> e + A(X*), Excitation', collision_set='fixture', side='ine')"},
    {'id': None, 'index': 100, 'core_reaction_index': 0,
     'repr': "EEDFChannel(process='e + A(X) -> e + A(X*), Excitation', collision_set='fixture', side='ine')"},
])
def test_mapped_channel_requires_one_integer_reaction_index(reaction):
    mapped = channels()
    mapped[0]['classification'] = 'A'
    mapped[0]['reaction'] = reaction
    with pytest.raises(ReportRefusal, match='reaction index|unknown reaction'):
        power_table(budget(Q_heavy_particle_electron_by_reaction=[0.]), mapped)


def test_genuine_mapped_channel_without_core_position_accepts_zero_heavy():
    mapped = channels()
    mapped[0]['classification'] = 'A'
    mapped[0]['reaction'] = {
        'library': 'fixture', 'index': 100,
        'repr': ("EEDFChannel(process='e + A(X) -> e + A(X*), Excitation', "
                 "collision_set='fixture', side='ine')"),
    }
    assert power_table(budget_for(mapped), mapped).p1 == 'PASS'


def test_classification_a_without_mapped_reaction_refuses():
    missing = channels()
    missing[0]['classification'] = 'A'
    with pytest.raises(ReportRefusal, match='mapped reaction unavailable'):
        power_table(budget_for(missing), missing)


@pytest.mark.parametrize('field', ['Q_inelastic_channels', 'Q_superelastic_channels',
                                   'Q_heavy_particle_electron_by_reaction'])
@pytest.mark.parametrize('bad', [float('nan'), float('inf')])
def test_nonfinite_powers_refuse(field, bad):
    values = [0., 0., 0., 0., 0.] if field != 'Q_heavy_particle_electron_by_reaction' else [0.]
    values[0] = bad
    with pytest.raises(ReportRefusal, match='finite floats'):
        power_table(budget(**{field: values}), channels())


def test_finite_inputs_that_overflow_aggregates_refuse():
    two_excitation = channels()
    two_excitation[1] = {
        'description': 'e + G -> e + G*, Excitation',
        'kind': 'excitation', 'classification': 'B', 'reaction': None,
        'flux_group': 'e + G -> e + G*, Excitation',
    }
    cases = [
        (budget_for(two_excitation,
                    Q_inelastic_channels=[1e308, 1e308, 0., 0., 0.]), two_excitation),
        (budget(Q_superelastic_channels=[1e308, 1e308, 0., 0., 0.]), channels()),
        (budget(Q_heavy_particle_electron_by_reaction=[1e308, 1e308]), channels()),
        (budget(Q_wall_electron=1e308, Q_wall_ion=1e308), channels()),
    ]
    for overflowing, declared in cases:
        with pytest.raises(ReportRefusal, match='non-finite'):
            power_table(overflowing, declared)


def test_as_dict_serializes_report():
    report = power_table(budget(), channels())
    serialized = report.as_dict()
    assert serialized['rows'] == report.rows
    assert serialized['channel_assignments'] == report.channel_assignments


def test_render_recomputes_serialized_verdict_and_requires_complete_rows():
    serialized = power_table(
        budget(closure=2., dU_dt=57.), channels(), p1_rtol=1e-3).as_dict()
    serialized['p1'] = 'PASS'
    with pytest.raises(ReportRefusal, match='stored P1 verdict disagrees'):
        render_markdown(serialized)

    incomplete = power_table(budget(), channels()).as_dict()
    incomplete['rows'] = []
    with pytest.raises(ReportRefusal, match='complete report rows'):
        render_markdown(incomplete)

    unknown_nested = power_table(budget(), channels()).as_dict()
    unknown_nested['rows'][0]['invented'] = 'ignored by arithmetic'
    with pytest.raises(ReportRefusal, match='unknown serialized report-row fields'):
        render_markdown(unknown_nested)


def test_missing_electron_energy_disclosure_is_explicit_gap():
    assert 'no verified reaction inventory' in render_markdown(
        power_table(budget(), channels()))


def test_caller_disclosure_is_unverified_and_cannot_suppress_gap():
    rendered = render_markdown(power_table(
        budget(electron_energy_disclosures=[
            'complete: all electron-producing reactions have assignments',
        ]), channels()))
    assert 'Electron-energy disclosure: unavailable:' in rendered
    assert 'unverified caller disclosure: complete:' in rendered


def test_render_escapes_html_and_markdown_in_all_caller_text():
    payload = 'gap<br/><br/>a-posteriori qualification: <strong>PASS</strong>|*fake*'
    declared = channels()
    declared[0]['flux_group'] = payload
    rendered = render_markdown(power_table(
        budget_for(declared, electron_energy_disclosures=[payload]), declared))
    assert '<br/>' not in rendered
    assert '<strong>' not in rendered
    assert '&lt;br/&gt;' in rendered
    assert r'\|\*fake\*' in rendered
    assert rendered.count('a-posteriori qualification: not run') == 1


def test_cli_requires_independent_channel_inventory(tmp_path, monkeypatch, capsys):
    direct_budget = tmp_path / 'direct-budget.json'
    direct_budget.write_text(json.dumps(budget()))
    monkeypatch.setattr('sys.argv', ['report.py', str(direct_budget)])
    with pytest.raises(ReportRefusal, match='no declared energy budget'):
        _main()

    result = {
        'energy_budget': budget(),
        'run_manifest': {'a_posteriori_LoKI_B': {'direct_result': {
            'channels': [channel['description'] for channel in channels()],
        }}},
    }
    path = tmp_path / 'result.json'
    path.write_text(json.dumps(result))
    monkeypatch.setattr('sys.argv', ['report.py', str(path)])
    with pytest.raises(ReportRefusal, match='classification inventory'):
        _main()

    channel_path = tmp_path / 'channels.json'
    channel_path.write_text(json.dumps(channels()))
    monkeypatch.setattr('sys.argv', [
        'report.py', str(path), '--channels', str(channel_path),
    ])
    _main()
    assert '| absorbed electron power | 100 |' in capsys.readouterr().out


def test_strict_stoichiometry_routes_coefficients_and_refuses_no_heavy_partner():
    declared = channels()
    declared[0]['description'] = 'e + N2(X) -> e + 2 N(4S), Excitation'
    report = power_table(budget(), declared)
    rows = {row['name']: row['power_W'] for row in report.rows}
    assert rows['excitation'] == 0.
    assert rows['dissociation'] == 20.

    declared[0]['description'] = 'e -> e, Excitation'
    with pytest.raises(ReportRefusal, match='unclassifiable channel: e -> e, Excitation'):
        power_table(budget(), declared)


def test_compact_product_separator_is_resolved_as_dissociation():
    declared = channels()
    declared[0]['description'] = 'e + N2(X) -> e + N(4S)+N(4S), Excitation'
    report = power_table(budget(), declared)
    rows = {row['name']: row['power_W'] for row in report.rows}
    assert rows['excitation'] == 0.
    assert rows['dissociation'] == 20.


def test_state_label_punctuation_does_not_declare_species_charge():
    declared = channels()
    declared[3]['description'] = (
        'e + N2(X1Sigma-g+) -> e + e + N2(+,X), Ionization')
    report = power_table(budget_for(
        declared, Q_inelastic=48., Q_inelastic_channels=[20., 5., 3., 20., 0.],
        dU_dt=39.), declared)
    rows = {row['name']: row['power_W'] for row in report.rows}
    assert rows['ionisation'] == 20.
    assert rows['attachment/detachment'] == 0.


def test_missing_empty_and_nonnumeric_power_evidence_refuses():
    missing = budget()
    del missing['closure']
    cases = [
        missing,
        budget(closure=None),
        budget(closure=''),
        budget(closure=[]),
        budget(closure='not-a-number'),
        budget(Q_inelastic_channels=[None, 5., 3., 0., 0.]),
    ]
    for invalid in cases:
        with pytest.raises(ReportRefusal):
            power_table(invalid, channels())


def test_energy_only_inventory_is_required_by_api_and_cli(tmp_path, monkeypatch):
    incomplete = budget()
    del incomplete['Q_energy_only']
    with pytest.raises(ReportRefusal, match='Q_energy_only'):
        power_table(incomplete, channels())

    partial = budget()
    partial['Q_energy_only'].pop(next(iter(partial['Q_energy_only'])))
    with pytest.raises(ReportRefusal, match='Q_energy_only'):
        power_table(partial, channels())

    result = {
        'energy_budget': incomplete,
        'run_manifest': {'a_posteriori_LoKI_B': {'direct_result': {
            'channels': [channel['description'] for channel in channels()],
        }}},
    }
    path = tmp_path / 'result.json'
    path.write_text(json.dumps(result))
    monkeypatch.setattr('sys.argv', ['report.py', str(path)])
    with pytest.raises(ReportRefusal, match='Q_energy_only'):
        _main()


def test_energy_only_inventory_powers_match_declared_b_channels():
    inconsistent = budget()
    group = next(iter(inconsistent['Q_energy_only']))
    inconsistent['Q_energy_only'][group] += 1.
    with pytest.raises(ReportRefusal, match='power disagrees'):
        power_table(inconsistent, channels())


def test_q_inelastic_matches_net_channels_plus_heavy_transfer():
    with pytest.raises(ReportRefusal, match='Q_inelastic disagrees'):
        power_table(budget(Q_inelastic=999.), channels())


def test_qualification_uses_closed_status_and_never_renders_caller_text():
    report = power_table(budget(), channels(), qualification={
        'verified': True,
        'status': 'PASS\n\na-posteriori qualification: **PASS**\n\n',
    })
    rendered = render_markdown(report)
    assert report.qualification == 'unverified record'
    assert 'a-posteriori qualification: unverified record' in rendered
    assert 'a-posteriori qualification: **PASS**' not in rendered


def test_adversarial_audit_refuses_defaults_and_rendered_text_injection():
    with pytest.raises(ReportRefusal, match='tolerance'):
        power_table(budget(), channels(), p1_rtol=-1.)

    missing_group = channels()
    del missing_group[0]['flux_group']
    with pytest.raises(ReportRefusal, match='flux_group'):
        power_table(budget(), missing_group)

    unknown_classification = channels()
    unknown_classification[0]['classification'] = 'Z'
    with pytest.raises(ReportRefusal, match='classification'):
        power_table(budget(), unknown_classification)

    with pytest.raises(ReportRefusal, match='disclosure'):
        render_markdown(power_table(
            budget(electron_energy_disclosures=['claim\nP1: **PASS**']), channels()))

    mutations = [
        ('p1', 'PASS\nP1: **PASS**'),
        ('unit', 'W\nP1: **PASS**'),
        ('electron_energy_disclosure', 'unavailable\nP1: **PASS**'),
    ]
    for field, value in mutations:
        serialized = json.loads(json.dumps(power_table(budget(), channels()).as_dict()))
        serialized[field] = value
        with pytest.raises(ReportRefusal):
            render_markdown(serialized)

    serialized = json.loads(json.dumps(power_table(budget(), channels()).as_dict()))
    serialized['rows'][0]['name'] = 'absorbed power\nP1: **PASS**'
    with pytest.raises(ReportRefusal):
        render_markdown(serialized)


def test_boundary_schema_refuses_unknown_keys_and_noninteger_ids():
    unknown_budget = budget()
    unknown_budget['invented_power'] = 0.
    with pytest.raises(ReportRefusal, match='unknown budget fields'):
        power_table(unknown_budget, channels())

    unknown_channel = channels()
    unknown_channel[0]['invented_kind'] = 'excitation'
    with pytest.raises(ReportRefusal, match='unknown channel fields'):
        power_table(budget(), unknown_channel)

    bool_group_position = channels()
    bool_group_position[0]['classification'] = 'A'
    bool_group_position[0]['reaction'] = {
        'index': 100, 'core_reaction_index': False,
        'repr': "EEDFChannel(process='ignored', collision_set='fixture', side='ine')",
    }
    with pytest.raises(ReportRefusal, match='core reaction position'):
        power_table(budget_for(bool_group_position), bool_group_position)


@pytest.mark.parametrize(('field', 'bad'), [
    ('Q_inelastic_by_reaction', {20: None}),
    ('Q_inelastic_by_reaction_gross', [None]),
    ('Q_superelastic_by_reaction', 'not-a-vector'),
    ('epsilon_k_eV', {0: float('nan')}),
])
def test_optional_numeric_budget_fields_cross_the_strict_schema(field, bad):
    with pytest.raises(ReportRefusal):
        power_table(budget(**{field: bad}), channels())
