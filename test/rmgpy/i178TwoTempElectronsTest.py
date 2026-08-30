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

"""
I-178 — :class:`rmgpy.kinetics.arrhenius.TwoTemperaturePlasma` carries a signed-net
``electrons`` field, so power-law plasma channels that change the free-electron count
(dissociative recombination ``Ar2+ + e- -> Ar + Ar``, three-body recombination) are
representable rather than refused.

The convention is established in the class docstring and pinned here: the field is the
same net quantity as :attr:`rmgpy.reaction.Reaction.electrons`, ``BadnellRRArrhenius``
and ``VoronovEIArrhenius``. It is NOT the "electrons transferred" quantity the
charge-transfer laws carry, which is why those stay refused as a second placement
source while this one is validated as a third source that must agree.

These tests cover the three verifier gates of that ticket:

* the loader (``KineticsLibrary.load``) refuses an electron-changing entry when the
  rate law carries no net count, and accepts it once the field declares one;
* the placement resolver validates a declared ``TwoTemperaturePlasma`` net against the
  family declaration, accepting agreement and refusing disagreement;
* a charge-transfer rate law's ``electrons`` stays fatal — the list was widened, not
  punctured.
"""

import os
import pickle

import pytest

from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.electron_placement import (
    _NET_ELECTRON_KINETICS_CLASSES,
    _declares_net_electron_count,
    resolve_electron_placement,
)
from rmgpy.exceptions import DatabaseError, ElectronPlacementError
from rmgpy.kinetics.arrhenius import ArrheniusChargeTransfer, TwoTemperaturePlasma
from rmgpy.species import Species

# ``Plasma_Electron_Attachment`` ships in FAMILY_ELECTRON_PLACEMENT as ``(1, 0)`` —
# one electron on the reactant side, net -1 — so a reaction attributed to it reaches
# the resolver without injecting any declaration.
ATTACHMENT_FAMILY = "Plasma_Electron_Attachment"


def _electron():
    return Species(label="e").from_adjacency_list("1 e u1 p0 c-1")


def _o2():
    return Species(label="O2").from_smiles("[O][O]")


def _o2_anion():
    return Species(label="O2-").from_smiles("[O][O-]")


def _two_temp(electrons=0, units="m^3/(mol*s)"):
    """A TwoTemperaturePlasma whose A-factor units fix its reaction order: order 2
    for ``m^3/(mol*s)``, which matches attachment's one-heavy-plus-one-electron view."""
    return TwoTemperaturePlasma(A=(1.0e12, units), n=0.0, Ea_g=(0.0, "J/mol"),
                                Ea_e=(0.0, "J/mol"), electrons=electrons)


# --------------------------------------------------------------------------- #
# The field itself                                                            #
# --------------------------------------------------------------------------- #
class TestTwoTemperaturePlasmaElectronField:

    def test_default_is_zero_generic_form_asserts_no_electron_change(self):
        k = TwoTemperaturePlasma(A=(1.0e12, "m^3/(mol*s)"), n=0.0)
        assert k.electrons.value_si == 0.0

    def test_declared_net_is_signed_and_stored(self):
        assert _two_temp(electrons=-1).electrons.value_si == -1.0
        assert _two_temp(electrons=2).electrons.value_si == 2.0

    def test_field_survives_pickle_roundtrip(self):
        k = _two_temp(electrons=-1)
        assert pickle.loads(pickle.dumps(k, -1)).electrons.value_si == -1.0

    def test_repr_emits_nonzero_electrons_only(self):
        # A non-zero count decides whether a persisted entry balances, so it must
        # round-trip through repr; the generic default of 0 is left off.
        assert "electrons" not in repr(_two_temp(electrons=0))
        assert "electrons=-1" in repr(_two_temp(electrons=-1))

    def test_electrons_distinguishes_identity(self):
        assert not _two_temp(electrons=-1).is_identical_to(_two_temp(electrons=0))
        assert _two_temp(electrons=-1).is_identical_to(_two_temp(electrons=-1))

    def test_registered_as_a_net_electron_kinetics_class(self):
        assert "TwoTemperaturePlasma" in _NET_ELECTRON_KINETICS_CLASSES
        assert _declares_net_electron_count(_two_temp(electrons=-1))


# --------------------------------------------------------------------------- #
# The loader gate: KineticsLibrary.load                                       #
# --------------------------------------------------------------------------- #
def _write_ar_dr_library(directory, declare_electrons):
    """Write a one-entry dissociative-recombination library ``Ar2+ + e- -> Ar + Ar``.

    ``Ar2+`` builds as ``[Ar][Ar+]`` (net charge +1, a doublet); the free electron is
    carried as ``Reaction.electrons``, exactly as radiative recombination carries its
    own. When ``declare_electrons`` is None the rate law leaves ``electrons`` at its
    generic default of 0, which is the state that cannot balance a +1 -> 0 charge
    change and is refused."""
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, "dictionary.txt"), "w") as f:
        f.write("[Ar2p]\nmultiplicity 2\n1 Ar u0 p3 c+1 {2,S}\n2 Ar u1 p3 c0 {1,S}\n\n"
                "[Ar]\n1 Ar u0 p4 c0\n")
    electrons_kw = "" if declare_electrons is None else f", electrons={declare_electrons}"
    with open(os.path.join(directory, "reactions.py"), "w") as f:
        f.write(
            'name = "ArDissociativeRecombination"\n'
            'shortDesc = u""\n'
            'longDesc = u""""""\n'
            'entry(\n'
            '    index = 0,\n'
            '    label = "[Ar2p] => [Ar] + [Ar]",\n'
            '    degeneracy = 1,\n'
            '    reversible = False,\n'
            f'    kinetics = TwoTemperaturePlasma(A=(1.0e19, "cm^3/(mol*s)"), n=-0.5,\n'
            f'        Ea_g=(0.0, "kJ/mol"), Ea_e=(0.0, "kJ/mol"){electrons_kw},\n'
            f'        Tmin=(300.0, "K"), Tmax=(30000.0, "K")),\n'
            '    shortDesc = u"placeholder rate, not sourced",\n'
            '    longDesc = u"""scratch probe only""",\n'
            ')\n'
        )


def _load(directory):
    lib = KineticsLibrary()
    lib.load(os.path.join(directory, "reactions.py"),
             local_context={"TwoTemperaturePlasma": TwoTemperaturePlasma})
    return list(lib.entries.values())[0].item


class TestLoaderGate:

    def test_undeclared_electron_change_is_refused_as_unbalanced(self, tmp_path):
        # The verifier's "refusal reproduces" gate: without a net count on the rate
        # law the electron cannot be carried, so the loader rejects the charge change.
        d = str(tmp_path / "undeclared")
        _write_ar_dr_library(d, declare_electrons=None)
        with pytest.raises(DatabaseError, match="was not balanced"):
            _load(d)

    def test_declared_electron_change_loads_with_correct_stoichiometry(self, tmp_path):
        # The verifier's "accepted with the correct stoichiometry" gate.
        d = str(tmp_path / "declared")
        _write_ar_dr_library(d, declare_electrons=-1)
        rxn = _load(d)
        assert [s.label for s in rxn.reactants] == ["[Ar2p]"]
        assert [s.label for s in rxn.products] == ["[Ar]", "[Ar]"]
        assert rxn.electrons == -1
        assert rxn.is_balanced()


# --------------------------------------------------------------------------- #
# The resolver gate                                                           #
# --------------------------------------------------------------------------- #
class TestResolverGate:

    def _attachment(self, kinetics):
        o2, o2m = _o2(), _o2_anion()
        reaction = TemplateReaction(
            reactants=[o2], products=[o2m], family=ATTACHMENT_FAMILY,
            electrons=-1, reversible=False, is_forward=True, kinetics=kinetics)
        return reaction, [_electron(), o2, o2m]

    def test_matching_two_temp_net_is_accepted_as_a_third_source(self):
        # electrons=-1 agrees with the family net of -1: the resolver builds the view
        # with the electron placed on the reactant side, unchanged canonical reaction.
        reaction, species = self._attachment(_two_temp(electrons=-1))
        view = resolve_electron_placement(reaction, species)
        assert [str(s) for s in view.reactants] == ["O2", "e"]
        assert [str(s) for s in view.products] == ["O2-"]
        assert view.electrons == 0  # the electron is now explicit, not metadata
        assert reaction.electrons == -1  # canonical reaction untouched

    def test_disagreeing_two_temp_net_is_refused(self):
        # electrons=+1 contradicts the family net of -1: the third source must AGREE.
        reaction, species = self._attachment(_two_temp(electrons=1))
        with pytest.raises(ElectronPlacementError, match="disagrees with the family"):
            resolve_electron_placement(reaction, species)

    def test_charge_transfer_electron_field_stays_fatal(self):
        # The widen-not-puncture gate: a charge-transfer law's ``electrons`` is a
        # different quantity (electrons transferred, a live Butler-Volmer parameter)
        # and is refused as a second placement source, exactly as before I-178.
        ct = ArrheniusChargeTransfer(A=(1.0e12, "m^3/(mol*s)"), n=0.0,
                                     Ea=(0.0, "kJ/mol"), electrons=-1)
        reaction, species = self._attachment(ct)
        with pytest.raises(ElectronPlacementError,
                           match="second placement source|double-represent"):
            resolve_electron_placement(reaction, species)
