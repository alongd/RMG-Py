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
I-135 -- the mechanism RMG writes can be read back, with the electron
temperature still attached.

RMG's Chemkin writer marks every plasma rate law with a ``TDEP/<electron>/``
auxiliary line. That is a real Chemkin construct, not an invention: CHEMKIN-III,
section IV, "Species Temperature Dependence" --

    When solving multi-fluid problems that involve multiple temperatures (for
    example, electron temperature and neutral gas temperature), the auxiliary
    information data may follow the reaction to specify the species on whose
    temperature the reaction depends. Here the species name follows the
    auxiliary keyword TDEP. This option causes the reaction rate constant to be
    evaluated using the specified species temperature and the rate parameters
    given in the reaction data.

-- and the manual's own worked argon/chlorine plasma mechanism writes exactly
``TDEP/E/`` on its own line beneath the reaction. What RMG emits conforms.

RMG's *reader* did not implement it. Its auxiliary-keyword dispatch had no case
for the name, so the line fell through to the collider-efficiency default and
the electron's species label was parsed as a number:

    ChemkinError: 'e-(1)' doesn't look like a collision efficiency for species
    TDEP in line 'TDEP/E-(1)/'

The important thing about the fix is what it must NOT do. Deleting the line, or
consuming it and ignoring it, makes ``load_chemkin_file`` succeed -- and returns
a plain :class:`~rmgpy.kinetics.Arrhenius` evaluated at the GAS temperature. At
the deck's own conditions (1000 K gas, 11 604.5 K electrons) that is not a
lossier reading of the same reaction, it is a different reaction, some ten
orders of magnitude slower. So these tests pin the rate the reloaded object
returns at two distinct temperatures, not merely that the file loads.

What the round trip still loses is pinned too, deliberately, in
:class:`TestWhatTheRoundTripStillLoses`: the electron count and the original
functional form do not come back. Both losses are the writer's, documented on
the ``TDEP`` line itself; neither is introduced here. A test that pins a loss is
how the loss stops being a surprise.
"""

import math

import pytest

from rmgpy.chemkin import (_plasma_arrhenius_for_chemkin, load_chemkin_file,
                           save_chemkin_file, save_species_dictionary)
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.exceptions import ChemkinError
from rmgpy.kinetics import (Arrhenius, BadnellRRArrhenius, ElectronCollisionPlasma,
                            TwoTemperaturePlasma, VoronovEIArrhenius)
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial

#: The deck's own conditions. The gap between them is the whole point: a rate
#: law read back on the gas temperature is wrong by orders of magnitude here.
T_GAS = 1000.0
T_ELECTRON = 11604.5

ADJACENCY_LISTS = {
    'e': '1 e u1 p0 c-1',
    'Li': '1 Li u1 p0 c0',
    'Liplus': '1 Li u0 p0 c+1',
}


def _flat_nasa():
    coeffs = [2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    return NASA(polynomials=[NASAPolynomial(coeffs=coeffs, Tmin=(100, 'K'), Tmax=(1000, 'K')),
                             NASAPolynomial(coeffs=coeffs, Tmin=(1000, 'K'), Tmax=(5000, 'K'))],
                Tmin=(100, 'K'), Tmax=(5000, 'K'))


def _species():
    """A fresh label -> Species map, indexed so the written labels are the
    ``label(index)`` spellings a real run produces."""
    out = {}
    for index, (label, adjlist) in enumerate(ADJACENCY_LISTS.items(), start=1):
        spc = Species(label=label).from_adjacency_list(adjlist)
        spc.index = index
        spc.thermo = _flat_nasa()
        out[label] = spc
    return out


def _ionisation(spc, kinetics=None):
    """``[Li] => [Lip]``, ``electrons = +1`` -- the canonical database form of
    electron-impact ionisation, which the writer expands to
    ``Li + e- => Li+ + 2 e-``."""
    return LibraryReaction(reactants=[spc['Li']], products=[spc['Liplus']],
                           library='PlasmaElectronImpactIonization', reversible=False,
                           electrons=1,
                           kinetics=kinetics if kinetics is not None else VoronovEIArrhenius(Z=3, N=3))


def _write(tmp_path, reaction, spc, name='chem.inp'):
    species_list = [spc['e'], spc['Li'], spc['Liplus']]
    chem = str(tmp_path / name)
    dictionary = str(tmp_path / 'species_dictionary.txt')
    save_chemkin_file(chem, species_list, [reaction], verbose=False,
                      check_for_duplicates=False)
    save_species_dictionary(dictionary, species_list)
    return chem, dictionary


def _rewrite_tdep_line(chem, replacement):
    """Replace the file's ``TDEP/.../`` line with `replacement`, in place."""
    with open(chem) as f:
        lines = f.read().splitlines(True)
    out = []
    replaced = False
    for line in lines:
        if line.strip().upper().startswith('TDEP/'):
            out.append(replacement + '\n')
            replaced = True
        else:
            out.append(line)
    assert replaced, 'no TDEP line in the written file'
    with open(chem, 'w') as f:
        f.writelines(out)
    return chem


class TestTheFileReadsBack:
    """The headline: the artifact RMG writes is one RMG can read."""

    def test_the_written_file_loads(self, tmp_path):
        spc = _species()
        chem, dictionary = _write(tmp_path, _ionisation(spc), spc)
        assert 'TDEP/e(1)/' in open(chem).read()

        _read_species, read_reactions = load_chemkin_file(chem, dictionary)
        assert len(read_reactions) == 1

    def test_the_rate_comes_back_on_the_electron_temperature(self, tmp_path):
        """Not merely "it loads". The reloaded rate law must be a function of
        Te, and must give the electron-temperature answer when the reactor asks
        it for one at Te != T."""
        spc = _species()
        reaction = _ionisation(spc)
        chem, dictionary = _write(tmp_path, reaction, spc)
        written, _note = _plasma_arrhenius_for_chemkin(reaction.kinetics)

        _read_species, read_reactions = load_chemkin_file(chem, dictionary)
        kinetics = read_reactions[0].kinetics

        assert isinstance(kinetics, TwoTemperaturePlasma), type(kinetics).__name__
        assert kinetics.uses_electron_temperature is True

        # The reactor's question. `rmgpy/solver/plasma.pyx` routes every
        # `uses_electron_temperature` rate law through this call.
        # 2%, not machine precision: the Chemkin reaction line carries A to
        # four significant figures and n and Ea to three decimals, so the trip
        # cannot be tighter than the format is. Same tolerance as
        # test/rmgpy/i126ChemkinElectronOrderTest.py's round trip.
        k_two_temp = kinetics.get_rate_coefficient_two_temp(T_GAS, T_ELECTRON)
        assert math.isclose(k_two_temp, written.get_rate_coefficient(T_ELECTRON), rel_tol=0.02), (
            k_two_temp, written.get_rate_coefficient(T_ELECTRON))

    def test_a_gas_temperature_reading_would_be_a_different_reaction(self, tmp_path):
        """The failure mode this ticket exists to avoid, made quantitative.

        Reading the same three numbers at the gas temperature is not a lossy
        answer, it is the wrong one, and the gap is enormous at plasma
        conditions. Pinned as a floor so that a future "just ignore the TDEP
        line" cannot slip past the test above."""
        spc = _species()
        reaction = _ionisation(spc)
        chem, dictionary = _write(tmp_path, reaction, spc)
        written, _note = _plasma_arrhenius_for_chemkin(reaction.kinetics)

        _read_species, read_reactions = load_chemkin_file(chem, dictionary)
        kinetics = read_reactions[0].kinetics

        k_correct = kinetics.get_rate_coefficient_two_temp(T_GAS, T_ELECTRON)
        k_gas_temperature = written.get_rate_coefficient(T_GAS)
        assert k_correct / k_gas_temperature > 1e6, (k_correct, k_gas_temperature)

    @pytest.mark.parametrize('kinetics', [
        VoronovEIArrhenius(Z=3, N=3),
        BadnellRRArrhenius(Z=3, N=2),
        ElectronCollisionPlasma(energies=([0.0, 1.0, 2.0, 5.0, 10.0], 'eV/molecule'),
                                sigma=([0.0, 1.0e-21, 3.0e-21, 5.0e-21, 4.0e-21], 'm^2')),
        TwoTemperaturePlasma(A=(1.0e10, 'cm^3/(mol*s)'), n=0.5,
                             Ea_g=(50.0, 'kJ/mol'), Ea_e=(50.0, 'kJ/mol')),
    ], ids=['voronov', 'badnell', 'electron-collision', 'two-temperature'])
    def test_every_plasma_rate_law_survives(self, tmp_path, kinetics):
        """The writer marks all four with TDEP, so the reader must handle all
        four -- not just the ionisation shape the deck happens to exercise."""
        spc = _species()
        chem, dictionary = _write(tmp_path, _ionisation(spc, kinetics=kinetics), spc)
        written, _note = _plasma_arrhenius_for_chemkin(kinetics)

        _read_species, read_reactions = load_chemkin_file(chem, dictionary)
        read_kinetics = read_reactions[0].kinetics

        assert read_kinetics.uses_electron_temperature is True
        assert math.isclose(read_kinetics.get_rate_coefficient_two_temp(T_GAS, T_ELECTRON),
                            written.get_rate_coefficient(T_ELECTRON), rel_tol=0.02)


class TestWhatTheRoundTripStillLoses:
    """Pinned, not hidden. Both losses are the Chemkin form's, and both were
    already there before this ticket; a test is how they stay visible."""

    def test_the_electron_count_does_not_come_back(self, tmp_path):
        """``Reaction.electrons`` is RMG's canonical scalar representation of
        the electron. The Chemkin equation carries the electrons explicitly
        instead, and nothing collapses them back on the way in, so the reloaded
        reaction has ``electrons == 0`` with electron species in its lists.

        Not fixed here: reaction identity and the electron representation are a
        separate ticket's, and collapsing them would change what the model
        builder sees."""
        spc = _species()
        chem, dictionary = _write(tmp_path, _ionisation(spc), spc)

        _read_species, read_reactions = load_chemkin_file(chem, dictionary)
        read = read_reactions[0]

        assert read.electrons == 0
        assert sum(1 for s in read.reactants if s.is_electron()) == 1
        assert sum(1 for s in read.products if s.is_electron()) == 2

    def test_the_original_functional_form_does_not_come_back(self, tmp_path):
        """Chemkin has no Voronov form. What the file carries is the
        modified-Arrhenius reduction of ``k(Te)``, and that is what comes back
        -- as the writer's own TDEP note says on the line beside the numbers."""
        spc = _species()
        reaction = _ionisation(spc)
        chem, dictionary = _write(tmp_path, reaction, spc)

        _read_species, read_reactions = load_chemkin_file(chem, dictionary)

        assert not isinstance(read_reactions[0].kinetics, VoronovEIArrhenius)
        assert 'not\n' not in open(chem).read()  # the note is one wrapped comment
        assert 'not representable in Chemkin' in open(chem).read()


class TestRefusalsRatherThanGuesses:
    """Everything the reader cannot rebuild faithfully raises. The alternative
    -- consume the line, keep the numbers, drop the temperature -- is the exact
    failure this ticket is about, so it must not be reachable by any input."""

    def test_another_auxiliary_keyword_on_the_line_is_refused(self, tmp_path):
        """``TDEP/E/ MOME`` is legal Chemkin (the manual's own example), and
        MOME changes what the rate coefficient means. RMG has no rate law for
        it, so the line must keep raising rather than be reduced to its TDEP
        part."""
        spc = _species()
        chem, dictionary = _write(tmp_path, _ionisation(spc), spc)
        _rewrite_tdep_line(chem, '    TDEP/e(1)/ MOME')

        with pytest.raises(ChemkinError) as exc:
            load_chemkin_file(chem, dictionary)
        # The message has to be RMG's own. Before the TDEP branch existed this
        # line also raised -- from the collider-efficiency default, complaining
        # that ' MOME' is not a number -- so asserting only that something
        # raised would pass against the defect this file exists to close.
        assert 'RMG understands TDEP/<species>/ on its own' in str(exc.value)

    def test_tdep_on_a_species_that_is_not_the_electron_is_refused(self, tmp_path):
        """Chemkin's TDEP is general over species temperatures. RMG's only
        non-gas temperature is the electron's, so anything else has no rate law
        behind it and must not be read back as one."""
        spc = _species()
        chem, dictionary = _write(tmp_path, _ionisation(spc), spc)
        _rewrite_tdep_line(chem, '    TDEP/Li(2)/')

        with pytest.raises(ChemkinError) as exc:
            load_chemkin_file(chem, dictionary)
        assert 'not the electron' in str(exc.value)

    def test_tdep_on_an_unknown_species_is_refused(self, tmp_path):
        spc = _species()
        chem, dictionary = _write(tmp_path, _ionisation(spc), spc)
        _rewrite_tdep_line(chem, '    TDEP/notaspecies(9)/')

        with pytest.raises(ChemkinError) as exc:
            load_chemkin_file(chem, dictionary)
        assert 'species dictionary' in str(exc.value)

    def test_a_malformed_tdep_line_is_refused(self, tmp_path):
        spc = _species()
        chem, dictionary = _write(tmp_path, _ionisation(spc), spc)
        _rewrite_tdep_line(chem, '    TDEP/')

        with pytest.raises(ChemkinError) as exc:
            load_chemkin_file(chem, dictionary)
        assert 'TDEP/<species>/' in str(exc.value)

    def test_tdep_combined_with_falloff_is_refused(self, tmp_path):
        """RMG has no pressure-dependent rate law evaluated at the electron
        temperature. Silently keeping the falloff and dropping the temperature
        is precisely the laundering this ticket removes."""
        spc = _species()
        chem, dictionary = _write(tmp_path, _ionisation(spc), spc)
        _rewrite_tdep_line(chem, '    TDEP/e(1)/\n    LOW/ 1.0e12 0.0 0.0 /')

        with pytest.raises(ChemkinError) as exc:
            load_chemkin_file(chem, dictionary)
        assert 'Lindemann' in str(exc.value)
