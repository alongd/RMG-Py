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
I-126 -- the export path resolves the electron from the same declaration the
reactor does.

The defect these tests pin: RMG carried TWO converters from the canonical
database representation (``Reaction.electrons`` as a scalar) to the explicit
one a solver needs.

* the reactor's -- :func:`rmgpy.electron_placement.resolve_electron_placement`,
  driven by the two-sided ``FAMILY_ELECTRON_PLACEMENT`` declaration;
* the export path's -- :func:`rmgpy.electron_balance.expand_electrons`, derived
  from the net scalar alone.

The net scalar under-determines placement. ``electrons = +1`` is equally
consistent with ``Li + e- => Li+ + 2 e-`` (order 2) and ``Li => Li+ + e-``
(order 1), and the net-derived converter always produces the second. For
electron-impact ionisation that is the wrong file -- the rate is off by a factor
of the electron density while the equation still balances in ``E`` -- and
``check_electron_reactant_order`` refused to write it. The guard was right; the
converter was wrong.

The fix makes the export converter consult the SAME declaration, so the two
halves of the codebase cannot disagree about which representation a reaction is
in. The tests below pin, in order: the ionisation shape now exports with its
incident electron; every one-sided declared owner exports EXACTLY as before (the
declaration and the net rule agree there, so nothing may move); an undeclared
charged reaction still gets the net-derived rule; a declaration that contradicts
the reaction's own net count is a named refusal rather than a laundered guess;
and a written Chemkin file survives a round trip with its rate intact.
"""

import math

import pytest

from rmgpy.chemkin import write_reaction_string, write_kinetics_entry
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.electron_balance import (expand_electrons, check_electron_balance,
                                    check_electron_reactant_order,
                                    get_placement_declaration)
from rmgpy.electron_placement import FAMILY_ELECTRON_PLACEMENT, resolve_electron_placement
from rmgpy.exceptions import MechanismWriterError
from rmgpy.kinetics import Arrhenius, BadnellRRArrhenius, VoronovEIArrhenius
from rmgpy.reaction import Reaction
from rmgpy.species import Species

IONISATION = 'PlasmaElectronImpactIonization'
RECOMBINATION = 'PlasmaRadiativeRecombination'

ADJACENCY_LISTS = {
    'e': '1 e u1 p0 c-1',
    'Li': '1 Li u1 p0 c0',
    'Liplus': '1 Li u0 p0 c+1',
}


def _species():
    """A fresh label -> Species map, indexed so ``get_species_identifier``
    produces the same ``label(index)`` spellings a real run writes."""
    out = {}
    for index, (label, adjlist) in enumerate(ADJACENCY_LISTS.items(), start=1):
        spc = Species(label=label).from_adjacency_list(adjlist)
        spc.index = index
        out[label] = spc
    return out


def _ionisation(spc, electrons=1):
    """``[Li] => [Lip]``, ``electrons = +1``: the canonical database form of
    electron-impact ionisation, exactly as ``KineticsLibrary.load`` builds it
    from ``PlasmaElectronImpactIonization/reactions.py``. The declaration for
    this owner is ``(1, 2)``."""
    return LibraryReaction(reactants=[spc['Li']], products=[spc['Liplus']],
                           library=IONISATION, reversible=False, electrons=electrons,
                           kinetics=VoronovEIArrhenius(Z=3, N=3))


def _recombination(spc, electrons=-1):
    """``[Lip] => [Li]``, ``electrons = -1``: the canonical form of radiative
    recombination. Its declaration is ``(1, 0)`` -- one-sided, so the net rule
    and the declaration agree and its export must not move."""
    return LibraryReaction(reactants=[spc['Liplus']], products=[spc['Li']],
                           library=RECOMBINATION, reversible=False, electrons=electrons,
                           kinetics=BadnellRRArrhenius(Z=3, N=2))


class TestExportConsultsTheDeclaration:
    """The headline: the export converter places from the family declaration."""

    def test_ionisation_expands_to_the_two_sided_placement(self):
        spc = _species()
        reactants, products = expand_electrons(_ionisation(spc), list(spc.values()))
        assert sum(1 for s in reactants if s.is_electron()) == 1, (
            'the incident electron must be on the reactant side -- it is what makes '
            'the exported equation second order')
        assert sum(1 for s in products if s.is_electron()) == 2
        assert len(reactants) == 2

    def test_export_view_is_the_reactor_view(self):
        """Not "equivalent to": the same objects, in the same order. This is the
        claim the ticket exists to make true -- the writer and the reactor read
        one declaration, so they cannot drift."""
        spc = _species()
        reaction = _ionisation(spc)
        species_list = list(spc.values())
        reactants, products = expand_electrons(reaction, species_list)
        view = resolve_electron_placement(reaction, species_list)
        assert [id(s) for s in reactants] == [id(s) for s in view.reactants]
        assert [id(s) for s in products] == [id(s) for s in view.products]

    def test_chemkin_writes_the_ionisation(self):
        spc = _species()
        equation = write_reaction_string(_ionisation(spc), species_list=list(spc.values()))
        assert equation == 'Li(2)+e(1)=>Liplus(3)+e(1)+e(1)', equation

    def test_cantera_writes_the_ionisation(self):
        from rmgpy.yaml_cantera2 import get_reaction_equation
        spc = _species()
        equation = get_reaction_equation(_ionisation(spc), list(spc.values()))
        assert equation == 'Li(2) + e(1) => Liplus(3) + e(1) + e(1)', equation

    def test_written_equation_passes_both_export_guards(self):
        """The guards are the file's correctness criteria, not obstacles: E must
        balance AND the reactant side must carry the electron the rate law is
        proportional to."""
        spc = _species()
        reaction = _ionisation(spc)
        reactants, products = expand_electrons(reaction, list(spc.values()))
        equation = write_reaction_string(reaction, species_list=list(spc.values()))
        check_electron_balance(reaction, reactants, products, equation)
        check_electron_reactant_order(reaction, reactants, equation)

    def test_a_factor_dimensionality_matches_the_exported_reactant_count(self):
        """The written rate and the written equation must agree about order, or
        the file is wrong by a factor of the electron density in exactly the way
        the guard exists to prevent.

        The Voronov coefficient is ``cm^3/(mol*s)``, second order, and the
        Chemkin line carries its modified-Arrhenius reduction converted with
        ``get_conversion_factor_from_si_to_cm_mol_s()``. That factor is
        ``1e6**(order-1)``, so a second-order number is the SI value times 1e6 --
        and the equation beside it must therefore have two reactants, which is
        also what ``write_kinetics_entry`` computes as ``num_reactants`` from
        ``expand_electrons``."""
        from rmgpy.chemkin import _plasma_arrhenius_for_chemkin
        spc = _species()
        reaction = _ionisation(spc)
        assert len(expand_electrons(reaction, list(spc.values()))[0]) == 2
        arrhenius, _note = _plasma_arrhenius_for_chemkin(reaction.kinetics)
        assert arrhenius.A.get_conversion_factor_from_si_to_cm_mol_s() == pytest.approx(1.0e6)
        entry = write_kinetics_entry(reaction, species_list=list(spc.values()), verbose=False)
        written = float(entry.split('\n')[0].split()[-3])
        assert written == pytest.approx(arrhenius.A.value_si * 1.0e6, rel=1e-3)


class TestOneSidedOwnersDoNotMove:
    """Negative control at the converter. For every declared owner whose
    declaration is one-sided, the declaration and the net rule agree by
    construction, so no exported equation may change."""

    @pytest.mark.parametrize('owner', sorted(o for o, d in FAMILY_ELECTRON_PLACEMENT.items()
                                             if 0 in d))
    def test_every_one_sided_declaration_agrees_with_the_net_rule(self, owner):
        reactant_count, product_count = FAMILY_ELECTRON_PLACEMENT[owner]
        spc = _species()
        net = product_count - reactant_count
        undeclared = expand_electrons(
            Reaction(reactants=[spc['Li']], products=[spc['Liplus']], electrons=net),
            list(spc.values()))
        declared = expand_electrons(
            TemplateReaction(reactants=[spc['Li']], products=[spc['Liplus']],
                             electrons=net, family=owner),
            list(spc.values()))
        assert [str(s) for s in declared[0]] == [str(s) for s in undeclared[0]]
        assert [str(s) for s in declared[1]] == [str(s) for s in undeclared[1]]

    def test_radiative_recombination_equation_is_unchanged(self):
        spc = _species()
        equation = write_reaction_string(_recombination(spc), species_list=list(spc.values()))
        assert equation == 'Liplus(3)+e(1)=>Li(2)', equation


class TestUndeclaredReactionsKeepTheNetRule:
    """An owner absent from the registry -- every charge-transfer reaction, and
    every reaction with no family attribution at all -- keeps today's
    net-derived placement exactly. This is what bounds the blast radius to the
    four declared owners."""

    def test_no_family_attribution(self):
        """A plain :class:`Reaction` carries no ``family`` attribute at all --
        which is also why the reactor's placement VIEW (a plain ``Reaction``)
        can never be re-expanded: it has no owner to look a declaration up
        under, and ``electrons = 0`` besides."""
        spc = _species()
        reaction = Reaction(reactants=[spc['Li']], products=[spc['Liplus']], electrons=-1,
                            kinetics=Arrhenius(A=(1e12, 'cm^3/(mol*s)'), n=0, Ea=(0, 'kJ/mol')))
        assert get_placement_declaration(reaction) is None
        reactants, products = expand_electrons(reaction, list(spc.values()))
        assert [str(s) for s in reactants] == ['Li(2)', 'e(1)']
        assert [str(s) for s in products] == ['Liplus(3)']

    def test_undeclared_family(self):
        spc = _species()
        reaction = TemplateReaction(reactants=[spc['Li']], products=[spc['Liplus']],
                                    electrons=-1,
                                    family='Some_Family_With_No_Placement_Declaration')
        assert get_placement_declaration(reaction) is None
        reactants, _ = expand_electrons(reaction, list(spc.values()))
        assert [str(s) for s in reactants] == ['Li(2)', 'e(1)']

    def test_a_neutral_reaction_is_untouched(self):
        spc = _species()
        reaction = _ionisation(spc, electrons=0)
        reactants, products = expand_electrons(reaction, list(spc.values()))
        assert [str(s) for s in reactants] == ['Li(2)']
        assert [str(s) for s in products] == ['Liplus(3)']


class TestContradictionIsNamed:
    """A declared owner whose reaction carries a net count the declaration does
    not predict is corrupt. The export path refuses it by name rather than
    falling back to the net rule -- a silent fallback is how a wrong file gets
    written while the export reports success."""

    def test_net_count_contradicting_the_declaration_raises(self):
        spc = _species()
        reaction = _ionisation(spc, electrons=2)
        with pytest.raises(MechanismWriterError) as exc:
            expand_electrons(reaction, list(spc.values()))
        assert IONISATION in str(exc.value)
        assert '(1, 2)' in str(exc.value)

    def test_the_message_names_both_numbers(self):
        spc = _species()
        with pytest.raises(MechanismWriterError) as exc:
            expand_electrons(_ionisation(spc, electrons=-1), list(spc.values()))
        message = str(exc.value)
        assert 'electrons=-1' in message
        assert 'net +1' in message


class TestChemkinRoundTrip:
    """Write the mechanism, read it back, and show the electron placement and
    the rate both survive. A green writer is not a correct file.

    The trip stops one step short of ``load_chemkin_file``, and the reason is a
    defect this ticket did not introduce and does not own: RMG's Chemkin READER
    has no case for the ``TDEP/<electron>/`` auxiliary line RMG's own Chemkin
    WRITER emits for every plasma rate law. ``test_tdep_line_is_what_blocks_the
    _full_trip`` pins that separately, on a reaction whose exported equation is
    byte-identical before and after this ticket's change -- which is how it is
    established as pre-existing rather than caused here. Everything the format
    itself supports is round-tripped below.
    """

    def _write(self, tmp_path, reaction, species_list):
        from rmgpy.chemkin import save_chemkin_file, save_species_dictionary
        for s in species_list:
            s.thermo = _flat_nasa()
        chem = str(tmp_path / 'chem.inp')
        dictionary = str(tmp_path / 'species_dictionary.txt')
        save_chemkin_file(chem, species_list, [reaction], verbose=False,
                          check_for_duplicates=False)
        save_species_dictionary(dictionary, species_list)
        return chem, dictionary

    def test_round_trip(self, tmp_path):
        from rmgpy.chemkin import load_chemkin_file

        spc = _species()
        reaction = _ionisation(spc)
        species_list = [spc['e'], spc['Li'], spc['Liplus']]
        chem, dictionary = self._write(tmp_path, reaction, species_list)

        text = open(chem).read()
        assert 'Li(2)+e(1)=>Liplus(3)+e(1)+e(1)' in text, text
        assert 'TDEP/e(1)/' in text, text

        # Drop the TDEP auxiliary line -- what any Chemkin parser without plasma
        # support sees, and what RMG's own reader would see if it had a case for
        # it. Everything this ticket is about is on the two lines that remain.
        with open(chem, 'w') as f:
            f.write('\n'.join(line for line in text.split('\n')
                              if not line.strip().startswith('TDEP/')))

        _read_species, read_reactions = load_chemkin_file(chem, dictionary)
        assert len(read_reactions) == 1
        read = read_reactions[0]
        assert sum(1 for s in read.reactants if s.is_electron()) == 1, (
            'the incident electron did not survive the round trip')
        assert sum(1 for s in read.products if s.is_electron()) == 2
        assert len(read.reactants) == 2, 'the reaction came back at the wrong order'

        # The rate survives the trip. What the file carries is NOT the Voronov
        # rate law -- Chemkin has no such form, and the writer says so in the
        # TDEP note -- but its modified-Arrhenius reduction along T = Te. So the
        # round trip is measured against that reduction, to the four significant
        # figures the Chemkin line has room for. Comparing against the Voronov
        # coefficient instead would measure the REDUCTION's fit error (up to ~14%
        # over its 11 604 K - 2.3e8 K window), which the round trip neither
        # causes nor can fix.
        from rmgpy.chemkin import _plasma_arrhenius_for_chemkin
        written, _note = _plasma_arrhenius_for_chemkin(reaction.kinetics)
        for Te in (12000.0, 20000.0, 100000.0):
            assert math.isclose(read.kinetics.get_rate_coefficient(Te),
                                written.get_rate_coefficient(Te),
                                rel_tol=0.02), (
                Te, read.kinetics.get_rate_coefficient(Te), written.get_rate_coefficient(Te))

    def test_tdep_line_is_what_blocks_the_full_trip(self, tmp_path):
        """Pre-existing, and not this ticket's: the same refusal appears on
        radiative recombination, whose declaration is one-sided and whose
        exported equation is therefore unchanged by this ticket."""
        from rmgpy.chemkin import load_chemkin_file
        from rmgpy.exceptions import ChemkinError

        spc = _species()
        chem, dictionary = self._write(tmp_path, _recombination(spc),
                                       [spc['e'], spc['Li'], spc['Liplus']])
        assert 'Liplus(3)+e(1)=>Li(2)' in open(chem).read()
        with pytest.raises(ChemkinError) as exc:
            load_chemkin_file(chem, dictionary)
        assert 'TDEP' in str(exc.value)


def _flat_nasa():
    from rmgpy.thermo import NASA, NASAPolynomial
    coeffs = [2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    return NASA(polynomials=[NASAPolynomial(coeffs=coeffs, Tmin=(100, 'K'), Tmax=(1000, 'K')),
                             NASAPolynomial(coeffs=coeffs, Tmin=(1000, 'K'), Tmax=(5000, 'K'))],
                Tmin=(100, 'K'), Tmax=(5000, 'K'))
