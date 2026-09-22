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
Estimation of a species whose SATURATED FORM no atom type owns.

Estimators that work by hydrogen bond increments -- thermo, transport, solute data -- build a
saturated copy of a radical as an intermediate. That step assumes every radical electron is an
unfilled chemical valence. For an electronically excited atom it is not: metastable argon,
``Ar u2 p3 c0``, saturates to ArH2, which is electron-count-consistent and owned by no atom
type, so ``update_atomtypes`` raises ``AtomTypeError`` many frames below whatever the caller
asked for.

These tests pin the two things that were wrong about that:

1. transport estimation used to die on such a species, taking down ordinary ``tran.dat``
   rendering and pressure-dependent collision calculations with it. It now falls back to the
   Lennard-Jones parameters it already falls back to for every other failure.
2. neither failure said what had happened. Both now name the species, its saturated form, the
   estimate that could not be made, and what to supply instead.

Ground-state argon is the control that keeps this honest: it is the same element, it survives
today, and it must keep surviving -- ``u0`` saturates to itself, so it never enters this path.
An ordinary radical is the second control.
"""

import os

import pytest

from rmgpy import settings
from rmgpy.data.base import saturate_for_estimation
from rmgpy.data.solvation import SolvationDatabase
from rmgpy.data.transport import TransportDatabase
from rmgpy.exceptions import AtomTypeError, SaturatedStructureError
from rmgpy.molecule import Molecule
from rmgpy.species import Species

#: metastable argon, Ar(4s 3P2): two unpaired electrons that are not valences
AR_METASTABLE = "1 Ar u2 p3 c0"
#: ground state argon, which saturates to itself and must be unaffected by any of this
AR_GROUND = "1 Ar u0 p4 c0"
METHYL = """
1 C u1 p0 c0 {2,S} {3,S} {4,S}
2 H u0 p0 c0 {1,S}
3 H u0 p0 c0 {1,S}
4 H u0 p0 c0 {1,S}
"""


def molecule(adjlist):
    return Molecule().from_adjacency_list(adjlist)


def species(adjlist):
    return Species(molecule=[molecule(adjlist)])


class TestSaturateForEstimation:
    """The shared helper, with no database loaded."""

    def test_ordinary_radical_saturates(self):
        """A radical whose saturated form is a molecule is unaffected by the guard."""
        saturated, added = saturate_for_estimation(molecule(METHYL), "thermodynamic data")
        assert saturated.is_isomorphic(Molecule().from_smiles("C"))
        assert len(added) == 1

    def test_input_is_not_mutated(self):
        """The helper saturates a copy, as the code it replaced did."""
        methyl = molecule(METHYL)
        saturate_for_estimation(methyl, "thermodynamic data")
        assert methyl.get_radical_count() == 1

    def test_unownable_saturated_form_raises_named_error(self):
        """The whole point: a diagnosis, not an AtomTypeError from three frames down."""
        with pytest.raises(SaturatedStructureError) as exc_info:
            saturate_for_estimation(molecule(AR_METASTABLE), "transport data")
        message = str(exc_info.value)
        # what failed, for which species, and what to do instead
        assert "transport data" in message
        assert "no RMG atom type describes" in message
        assert "library entry" in message
        assert "Ar u2 p3 c0" in message, "the message must show the species"
        assert "Ar u0 p3 c0" in message, "the message must show the saturated form"

    def test_named_error_is_not_an_atom_type_error(self):
        """
        A modeller filtering on AtomTypeError would still be told nothing. The new exception
        deliberately does not inherit from it.
        """
        assert not issubclass(SaturatedStructureError, AtomTypeError)


class TestTransportOfUnsaturableSpecies:
    @classmethod
    def setup_class(cls):
        cls.database = TransportDatabase()
        cls.database.load(os.path.join(settings["database.directory"], "transport"), libraries=None)

    def test_estimation_falls_back_instead_of_crashing(self):
        """
        This is the defect: standard tran.dat rendering and pdep collision calculations both
        reach this call, and a family-level forbidden group cannot intercept either.
        """
        transport = self.database.get_transport_properties(species(AR_METASTABLE))[0]
        assert transport is not None
        assert transport.shapeIndex == 0, "a monatomic species is a sphere"
        assert transport.sigma.value_si > 0
        assert transport.epsilon.value_si > 0

    def test_the_fallback_says_it_is_a_fallback(self):
        """An unexplained number is the failure mode this replaces, not an improvement on it."""
        transport = self.database.get_transport_properties(species(AR_METASTABLE))[0]
        assert "saturated form" in transport.comment
        assert "guess" in transport.comment

    def test_group_estimate_still_raises_for_whoever_asks_directly(self):
        """
        The fallback lives in get_transport_properties. The group estimator itself must keep
        reporting that it could not do the job, or a caller of the inner method gets silence.
        """
        with pytest.raises(SaturatedStructureError):
            self.database.get_transport_properties_via_group_estimates(species(AR_METASTABLE))

    def test_ground_state_argon_is_unaffected(self):
        """Control: same element, no unpaired electrons, never saturated."""
        transport = self.database.get_transport_properties(species(AR_GROUND))[0]
        assert "saturated form" not in transport.comment

    def test_ordinary_radical_is_unaffected(self):
        """
        Control: the fallback must not have become the general answer for radicals. 1-butyl
        saturates to butane and gets a real Joback estimate off it.

        It is 1-butyl rather than CH3 on purpose: the transport group tree has no data at or
        above ``C_centered``, so a methyl group estimate raises KeyError and reaches the same
        Lennard-Jones fallback for an unrelated reason -- which would make this control unable
        to tell the two fallbacks apart.
        """
        transport = self.database.get_transport_properties_via_group_estimates(
            Species().from_smiles("[CH2]CCC"))[0]
        assert "Joback" in transport.comment
        assert "guess" not in transport.comment


class TestSoluteDataOfUnsaturableSpecies:
    @classmethod
    def setup_class(cls):
        cls.database = SolvationDatabase()
        cls.database.load(os.path.join(settings["database.directory"], "solvation"))

    def test_estimation_names_the_cause(self):
        """
        There is no defensible solute estimate for such a species, so this one reports rather
        than invents -- but it reports the species and the remedy, which AtomTypeError did not.
        """
        with pytest.raises(SaturatedStructureError) as exc_info:
            self.database.get_solute_data(species(AR_METASTABLE))
        assert "solute data" in str(exc_info.value)

    def test_ground_state_argon_is_unaffected(self):
        """Control."""
        assert self.database.get_solute_data(species(AR_GROUND)) is not None

    def test_ordinary_radical_is_unaffected(self):
        """Control: HBI estimation for real radicals is untouched."""
        assert self.database.get_solute_data(species(METHYL)) is not None
