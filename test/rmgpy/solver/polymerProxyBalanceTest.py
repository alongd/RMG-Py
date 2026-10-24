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

"""
The solver refuses, by name, a pool-coupled row whose proxies do not balance:
- the proxies change the row's net charge (``unbalanced:charge``);
- the proxies read another reaction channel than the one stamped from the
  row's raw molecules at generation (``unbalanced:radicals``), unless the row
  folds into the pool it reacts from (a chip, a chain-end fold), which drops
  the chain's radical by design.
Refused rows are listed in ``proxy_balance_census``. The check runs before the
kernel supersession, so a defective row is named for its defect.

Hand-built rows on small test proxies (test data only).
"""

import logging

import numpy as np
import pytest

import rmgpy.solver.polymer as sp_mod
from rmgpy.reaction import Reaction
from rmgpy.solver.polymer import HybridPolymerSystem

from polymerChannelIdentityTest import KIN, _pool, _spc, _stamp

red = pytest.mark.xfail(strict=True, raises=AssertionError,
                        reason='no proxy-balance guard: an unbalanced pool-coupled row runs live')


def _system(t_proxy, x='[CH2]CC', channel=None, fold=False, kernel=None):
    """Pool A (closed proxy) decomposes: A -> T + X, or as a fold A -> A + X."""
    core = [_spc('N#N', 'N2'), _spc(x, 'X')]
    pools = []
    for n, (label, proxy) in enumerate((('A', 'CCCCC'), ('T', t_proxy))):
        stem = 'C' * (n + 1)
        base = len(core)
        core += [_spc(proxy, label), _spc(stem + 'O', f'{label}_mu0'),
                 _spc(stem + 'N', f'{label}_mu1'), _spc(stem + 'S', f'{label}_mu2')]
        kw = {'k_scission': 0.5} if (label == 'A' and kernel == 'k_scission') else {}
        pools.append(_pool(label, base + 1, **kw))
    sp = {s.label: s for s in core}
    if fold:
        rxn = Reaction(reactants=[sp['A']], products=[sp['A'], sp['X']], kinetics=KIN, reversible=False)
        rxn.polymer_chip_units = 1
        rxn.polymer_flux_archetype = sp_mod.FLUX_DISCRETE_CHIP
    else:
        rxn = Reaction(reactants=[sp['A']], products=[sp['T'], sp['X']], kinetics=KIN, reversible=False)
        rxn.polymer_eject_units = 1.0
        rxn.polymer_flux_archetype = sp_mod.FLUX_VOLATILE_EJECTION
    if channel is not None:
        _stamp(rxn, channel)
    mask = np.array([s.label in ('N2', 'X') for s in core], dtype=bool)
    rs = HybridPolymerSystem(T=700.0, P=1.0e5, initial_mole_fractions={core[0]: 1.0}, V_poly=1.0,
                             polymer_pools=pools, mass_transfer=[], gas_species_mask=mask,
                             constant_gas_volume=False, initial_polymer_moments={'A': (1.0, 5.0, 30.0)},
                             termination=[])
    rs.initialize_model(core, [rxn], [], [])
    return rs, rxn


def _census(rs):
    return getattr(rs, 'proxy_balance_census', None)


@red
def test_a_row_whose_proxy_carries_a_charge_is_refused_by_name(caplog):
    from rmgpy import polymer
    getattr(polymer, '_proxy_balance_warned', set()).clear()
    with caplog.at_level(logging.WARNING):
        rs, rxn = _system('CCCC[CH+]C', channel='bond_fission')
    assert (rs.reaction_refused[0], rs.reaction_flux_archetype[0]) == (1, sp_mod.FLUX_NONE)
    assert _census(rs) == [{'reaction': str(rxn), 'reason': 'unbalanced:charge',
                            'charge': [0, 1], 'channel': ['bond_fission', 'unidentified']}]
    assert sum('PROXY BALANCE CENSUS' in r.getMessage() for r in caplog.records) == 1


@red
def test_a_row_whose_proxies_lose_the_stamped_channel_is_refused_by_name():
    """A fission stamped from the raw molecules whose pool product proxy is
    closed-shell: one radical appears from none."""
    rs, rxn = _system('CCCCCC', channel='bond_fission')
    assert (rs.reaction_refused[0], rs.reaction_flux_archetype[0]) == (1, sp_mod.FLUX_NONE)
    assert _census(rs) == [{'reaction': str(rxn), 'reason': 'unbalanced:radicals',
                            'charge': [0, 0], 'channel': ['bond_fission', 'unidentified']}]


@red
def test_the_defect_is_named_before_a_kernel_supersedes_the_row():
    rs, rxn = _system('CCCC[CH+]C', channel='bond_fission', kernel='k_scission')
    assert rs.reaction_refused[0] == 1
    assert [c['reason'] for c in _census(rs) or []] == ['unbalanced:charge']
    assert rs.kernel_supersession_census == []


def test_a_balanced_fission_stays_live():
    rs, _ = _system('[CH2]CCCCC', channel='bond_fission')
    assert (rs.reaction_refused[0], rs.reaction_flux_archetype[0]) == (0, sp_mod.FLUX_VOLATILE_EJECTION)
    assert not _census(rs)


def test_a_fold_into_its_own_pool_stays_live():
    """A chip keeps the chain in the pool it reacts from: the chain end's
    radical is dropped by design, not by a defect."""
    rs, _ = _system('[CH2]CCCCC', channel='bond_fission', fold=True)
    assert (rs.reaction_refused[0], rs.reaction_flux_archetype[0]) == (0, sp_mod.FLUX_DISCRETE_CHIP)
    assert not _census(rs)


def test_an_unstamped_row_is_checked_for_charge_only():
    rs, _ = _system('CCCCCC', channel=None)
    assert rs.reaction_refused[0] == 0
    assert not _census(rs)
