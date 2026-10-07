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
Depository membership is provenance, not a quality ranking.

``training`` does not mean "more accurate" and ``NIST`` does not mean "preferred"; neither
do filesystem order, dictionary order, load order or branch history.  When two
independently sourced depository entries supply kinetics for the same oriented reaction,
the engine must refuse to choose, and say what it is refusing to choose between.

The behaviour these tests pin, and the file that was here before them, are opposites.  The
previous contract was "``training`` wins", enforced by a load order; the tests below assert
that nothing wins and that the run stops with both candidates named.  There is therefore no
version of the engine that passes both, which is the point: the RED log for this file was
taken against the tip that still ranked ``training`` first
(``docs/depository-order/logs/30-conflict-fixture-RED-pre-repair.*``), because after the
repair the silent choice it documents does not exist to be documented.

**Assert on kinetics and provenance, never on a label.**  A test that only checks that the
returned depository is called ``training`` stays green if somebody swaps the two
directories' contents, which makes it a check that cannot fail.  Where a test needs to tie
a label to a rate, it pins the number from the database on disk, named in
``DATABASE_PROVENANCE`` below, and fails loudly if that number moves.
"""

import os
import re

import numpy as np
import pytest

from rmgpy import settings
from rmgpy.data.base import Entry
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.depository import KineticsDepository
from rmgpy.data.kinetics.family import KineticsFamily, TemplateReaction
from rmgpy.electron_balance import get_electron_placement_counts
from rmgpy.exceptions import DatabaseError
from rmgpy.kinetics import Arrhenius, ArrheniusEP, KineticsData, MultiArrhenius, PDepArrhenius, StickingCoefficient
from rmgpy.kinetics.uncertainties import RateUncertainty
from rmgpy.quantity import ScalarQuantity
from rmgpy.reaction import Reaction
from rmgpy.species import Species

# NOT imported at module level on purpose.  ``KineticsDepositoryConflictError`` does not
# exist on the tip that still ranked ``training`` first, and importing it up here turns the
# RED run into a single collection error -- one line that says nothing about which of the
# twelve behaviours below is missing, and that hides the eight tests which fail for their
# own reasons rather than for the exception's absence.


def _conflict_error():
    from rmgpy.exceptions import KineticsDepositoryConflictError
    return KineticsDepositoryConflictError

################################################################################
# Synthetic fixtures.  Two depositories, built in memory, that disagree about one
# reaction -- no database, no filesystem, no load order.
################################################################################

#: The reaction every synthetic fixture is about, as SMILES.  Ordinary thermal chemistry
#: on purpose: the conflict has nothing to do with plasma, surfaces or charge.
REACTANTS = ('[H]', 'C')  # H + CH4
PRODUCTS = ('[H][H]', '[CH3]')  # H2 + CH3

#: Two rates that differ by a factor of 100, so "which one came back" is a question about
#: the number and not about the label attached to it.
A_ALPHA = 1.0e6
A_BETA = 1.0e8


def _oriented_reaction(reversible=True, electrons=0, flip=False):
    reactants = [Species(smiles=s) for s in (PRODUCTS if flip else REACTANTS)]
    products = [Species(smiles=s) for s in (REACTANTS if flip else PRODUCTS)]
    return Reaction(reactants=reactants, products=products, reversible=reversible, electrons=electrons)


def _entry(index, a_factor, reversible=True, electrons=0, flip=False, rank=5,
           reference='Somebody et al., J. Made Up Data 1 (1999) 1', tmax=2000.0):
    return Entry(
        index=index,
        label='H + CH4 <=> H2 + CH3',
        item=_oriented_reaction(reversible=reversible, electrons=electrons, flip=flip),
        data=Arrhenius(A=(a_factor, 'cm^3/(mol*s)'), n=0.0, Ea=(40.0, 'kJ/mol'), T0=(1, 'K'),
                       Tmin=(300.0, 'K'), Tmax=(tmax, 'K')),
        rank=rank,
        reference=reference,
        reference_type='theory',
        short_desc='synthetic entry {0}'.format(index),
    )


def _depository(label, entries):
    depository = KineticsDepository(label=label, electrons=0)
    depository.entries = {entry.index: entry for entry in entries}
    return depository


def _family(depositories):
    """
    The smallest object ``get_kinetics`` will run against: a family that owns nothing but
    its depositories.  ``retrieve_template`` is replaced because the template is only used
    to decorate the kinetics comment, and building a real group tree here would add a
    dependency on the database that every one of these tests is designed not to have.
    """
    family = KineticsFamily(label='Synthetic')
    family.electrons = 0
    family.depositories = list(depositories)
    family.retrieve_template = lambda template_labels: []
    return family


def _ask(family, reaction=None):
    """Ask the family for kinetics exactly the way the model enlargement path does."""
    return family.get_kinetics(
        reaction if reaction is not None else _oriented_reaction(),
        template_labels=[],
        degeneracy=1,
        estimator='',
        return_all_kinetics=False,
    )


def _loaded_array_depositories(tmp_path):
    """Build two independent records through the public depository file loader."""
    depositories = []
    for name in ('training', 'NIST'):
        label = 'Synthetic/{0}'.format(name)
        entry = _entry(1, A_ALPHA)
        for species, species_label in zip(entry.item.reactants + entry.item.products,
                                          ('H', 'CH4', 'H2', 'CH3')):
            species.label = species_label
        entry.reference = None
        entry.data = KineticsData(
            Tdata=([300.0, 1000.0], 'K'),
            kdata=([1.0, 2.0], 's^-1', '*|/', [1.1, 1.2]),
        )
        source = _depository(label, [entry])
        root = tmp_path / name
        root.mkdir()
        path = root / 'reactions.py'
        source.save(str(path), reindex=False)

        loaded = KineticsDepository(label=label, electrons=0)
        loaded.load(str(path), KineticsDatabase().local_context, {})
        depositories.append(loaded)
    return depositories


def _unsupported_array_pair(case):
    if case == 'structured-object':
        class OpaqueValue:
            def __init__(self, payload):
                self.payload = payload

            def __eq__(self, other):
                return True

        dtype = np.dtype([('rate', 'f8'), ('payload', 'O')])
        return (np.array([(1.0, OpaqueValue('left-1')),
                          (2.0, OpaqueValue('left-2'))], dtype=dtype),
                np.array([(1.0, OpaqueValue('right-1')),
                          (2.0, OpaqueValue('right-2'))], dtype=dtype))
    if case == 'field-metadata':
        dtypes = [
            np.dtype({'names': ['rate'],
                      'formats': [np.dtype('float64', metadata={'source': source})]})
            for source in ('left', 'right')
        ]
        return tuple(np.array([(1.0,), (2.0,)], dtype=dtype) for dtype in dtypes)
    dtype = {'non-native-float64': '>f8',
             'float16': 'f2',
             'complex128': 'c16'}[case]
    return (np.array([1.0, 2.0], dtype=dtype),
            np.array([1.0, 2.0], dtype=dtype))


class TestStoredArrayDtypeAllowlist:
    @pytest.mark.parametrize('reverse_load_order', [False, True])
    @pytest.mark.parametrize(('case', 'storage_field'), [
        ('structured-object', 'uncertainty_si'),
        ('field-metadata', 'uncertainty_si'),
        ('non-native-float64', 'value_si'),
        ('non-native-float64', 'uncertainty_si'),
        ('float16', 'value_si'),
        ('float16', 'uncertainty_si'),
        ('complex128', 'value_si'),
        ('complex128', 'uncertainty_si'),
    ])
    def test_unreviewed_array_storage_is_not_assumed_identical(
            self, tmp_path, case, storage_field, reverse_load_order):
        depositories = _loaded_array_depositories(tmp_path)
        for depository, array in zip(depositories, _unsupported_array_pair(case)):
            setattr(depository.entries[1].data.kdata, storage_field, array)
        if reverse_load_order:
            depositories.reverse()

        with pytest.raises(_conflict_error()):
            _ask(_family(depositories))

    @pytest.mark.parametrize('storage_field', ['value_si', 'uncertainty_si'])
    def test_independently_loaded_plain_float64_arrays_are_identical(
            self, tmp_path, storage_field):
        depositories = _loaded_array_depositories(tmp_path)
        for depository in depositories:
            array = np.array([1.0, 2.0], dtype=np.float64)
            setattr(depository.entries[1].data.kdata, storage_field, array)

        kinetics, _, _, _ = _ask(_family(depositories))
        assert type(getattr(kinetics.kdata, storage_field)) is np.ndarray
        assert getattr(kinetics.kdata, storage_field).dtype == np.dtype('float64')


class TestTheEngineRefusesToRankDepositories:
    """§6.1-§6.3, §6.7-§6.11 -- the policy itself, on fixtures that own no filesystem."""

    def test_one_matching_source_loads_normally(self):
        """§6.1.  One source is not a conflict, and the rate that comes back is its rate."""
        family = _family([_depository('Synthetic/alpha', [_entry(1, A_ALPHA)])])
        kinetics, depository, entry, is_forward = _ask(family)
        assert is_forward is True
        assert entry.index == 1
        assert kinetics.A.value_si == pytest.approx(Arrhenius(A=(A_ALPHA, 'cm^3/(mol*s)')).A.value_si)

    def test_two_depositories_for_one_oriented_reaction_raise_and_name_both(self):
        """§6.2.  The headline: two sources, one reaction, no winner."""
        family = _family([
            _depository('Synthetic/training', [_entry(1, A_ALPHA)]),
            _depository('Synthetic/NIST', [_entry(7, A_BETA)]),
        ])
        with pytest.raises(_conflict_error()) as caught:
            _ask(family)
        message = str(caught.value)
        assert 'Synthetic/training' in message
        assert 'Synthetic/NIST' in message
        # Both *rates* are named, not merely both labels -- the user has to be able to see
        # what they are choosing between without going back to the database.
        assert '1e+06' in message and '1e+08' in message
        assert 'index 1' in message and 'index 7' in message

    @pytest.mark.parametrize('reverse_load_order', [False, True])
    def test_the_conflict_does_not_depend_on_load_order(self, reverse_load_order):
        """
        §6.3.  Reversing the order the depositories were loaded in changes nothing: the
        same error, naming the same two candidates.  Under the old policy this is exactly
        the knob that decided the rate.
        """
        depositories = [
            _depository('Synthetic/training', [_entry(1, A_ALPHA)]),
            _depository('Synthetic/NIST', [_entry(7, A_BETA)]),
        ]
        if reverse_load_order:
            depositories.reverse()
        with pytest.raises(_conflict_error()) as caught:
            _ask(_family(depositories))
        named = set(re.findall(r'Synthetic/\w+', str(caught.value)))
        assert named == {'Synthetic/training', 'Synthetic/NIST'}

    def test_two_explicit_opposite_irreversible_directions_stay_distinct(self):
        """
        §6.7.  One source declares ``A => B`` one-way, the other declares ``B => A``
        one-way.  Those are two processes, not two accounts of one, so this is not a
        conflict -- and the rate returned for a query oriented ``A => B`` is the one whose
        own declared direction is ``A => B``.  That is identity matching, not a ranking:
        swap which depository holds which direction and the answer swaps with it.
        """
        forward_in_training = _family([
            _depository('Synthetic/training', [_entry(1, A_ALPHA, reversible=False)]),
            _depository('Synthetic/NIST', [_entry(7, A_BETA, reversible=False, flip=True)]),
        ])
        kinetics, depository, entry, is_forward = _ask(forward_in_training)
        assert is_forward is True
        assert entry.index == 1
        assert depository.label == 'Synthetic/training'
        assert kinetics.A.value_si == pytest.approx(
            Arrhenius(A=(A_ALPHA, 'cm^3/(mol*s)')).A.value_si)
        assert 'Matched reaction 1' in kinetics.comment
        assert 'Synthetic/training' in kinetics.comment

        forward_in_nist = _family([
            _depository('Synthetic/training', [_entry(1, A_ALPHA, reversible=False, flip=True)]),
            _depository('Synthetic/NIST', [_entry(7, A_BETA, reversible=False)]),
        ])
        kinetics, depository, entry, is_forward = _ask(forward_in_nist)
        assert is_forward is True
        assert entry.index == 7
        assert depository.label == 'Synthetic/NIST'
        assert kinetics.A.value_si == pytest.approx(
            Arrhenius(A=(A_BETA, 'cm^3/(mol*s)')).A.value_si)
        assert 'Matched reaction 7' in kinetics.comment
        assert 'Synthetic/NIST' in kinetics.comment

    def test_different_electron_metadata_does_not_collapse_into_one_candidate(self):
        """
        §6.8.  A reaction that consumes an electron and one that does not are different
        reactions, however alike their heavy species look.  The entry declaring an electron
        must not be offered as a competing account of the electron-free query.
        """
        family = _family([
            _depository('Synthetic/training', [_entry(1, A_ALPHA, electrons=0)]),
            _depository('Synthetic/NIST', [_entry(7, A_BETA, electrons=-1)]),
        ])
        kinetics, depository, entry, is_forward = _ask(family, _oriented_reaction(electrons=0))
        assert entry.index == 1
        assert depository.label == 'Synthetic/training'
        assert is_forward is True
        assert kinetics.A.value_si == pytest.approx(
            Arrhenius(A=(A_ALPHA, 'cm^3/(mol*s)')).A.value_si)
        assert 'Matched reaction 1' in kinetics.comment
        assert 'Synthetic/training' in kinetics.comment

    def test_the_same_record_reached_twice_is_not_a_conflict(self):
        """
        §6.9/§3, first clause: literally one record, arrived at by two routes.  The same
        depository object listed twice is the degenerate case, and it must not be reported
        as two sources disagreeing with themselves.
        """
        depository = _depository('Synthetic/training', [_entry(1, A_ALPHA)])
        family = _family([depository, depository])
        kinetics, chosen, entry, is_forward = _ask(family)
        assert entry is depository.entries[1]

    def test_a_proven_identical_duplicate_is_not_a_conflict(self):
        """
        §6.9/§3, second clause: two records, proven identical in direction, kinetics class
        and parameters, validity range, rank and citation.  Deduplicating that is
        permitted, because there is provably one source and not two.
        """
        family = _family([
            _depository('Synthetic/training', [_entry(1, A_ALPHA)]),
            _depository('Synthetic/copy-of-training', [_entry(1, A_ALPHA)]),
        ])
        kinetics, depository, entry, is_forward = _ask(family)
        assert kinetics.A.value_si == pytest.approx(Arrhenius(A=(A_ALPHA, 'cm^3/(mol*s)')).A.value_si)

    def test_nearly_equal_stored_coefficients_are_not_exact_duplicates(self):
        """A repr-rounded coefficient is still a distinct source record."""
        family = _family([
            _depository('Synthetic/training', [_entry(1, A_ALPHA)]),
            # Keep every other compared field identical so only the stored coefficient
            # can make this a conflict.
            _depository('Synthetic/NIST', [_entry(1, A_ALPHA * (1 + 1e-8))]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    def test_adjacent_stored_si_coefficients_are_not_exact_duplicates(self):
        """Identity compares stored SI floats, not a unit-converted reconstruction."""
        first = _entry(1, A_ALPHA)
        second = _entry(1, A_ALPHA)
        first.data.A.value_si = 0.7865786249459419
        second.data.A.value_si = 0.786578624945942
        family = _family([
            _depository('Synthetic/training', [first]),
            _depository('Synthetic/NIST', [second]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    def test_equivalent_si_quantities_with_different_units_are_identical(self):
        first = _entry(1, A_ALPHA)
        second = _entry(1, A_ALPHA)
        second.data.A = (1.0, 'm^3/(mol*s)')
        second.data.A.value_si = first.data.A.value_si
        assert first.data.A.value_si == second.data.A.value_si
        family = _family([
            _depository('Synthetic/training', [first]),
            _depository('Synthetic/copy-of-training', [second]),
        ])
        kinetics, _, _, _ = _ask(family)
        assert kinetics.A.value_si == first.data.A.value_si

    @pytest.mark.parametrize('kinetics_type', [
        'Arrhenius',
        'KineticsData',
        'PDepArrhenius',
        'MultiArrhenius',
    ])
    def test_equal_si_values_with_incompatible_dimensions_are_not_identical(self, kinetics_type):
        """An equal float cannot erase whether a rate is first- or second-order."""
        def make_kinetics(units):
            if kinetics_type == 'Arrhenius':
                return Arrhenius(A=(1.0, units), n=0.0, Ea=(10.0, 'kJ/mol'))
            if kinetics_type == 'KineticsData':
                return KineticsData(Tdata=([300.0, 1000.0], 'K'),
                                    kdata=([1.0, 2.0], units))
            nested = Arrhenius(A=(1.0, units), n=0.0, Ea=(10.0, 'kJ/mol'))
            if kinetics_type == 'PDepArrhenius':
                return PDepArrhenius(pressures=([1.0], 'bar'), arrhenius=[nested])
            return MultiArrhenius(arrhenius=[nested])

        first = _entry(1, A_ALPHA)
        second = _entry(1, A_ALPHA)
        first.data = make_kinetics('s^-1')
        second.data = make_kinetics('m^3/(mol*s)')
        first_quantity = (first.data.kdata if kinetics_type == 'KineticsData'
                          else first.data.A if kinetics_type == 'Arrhenius'
                          else first.data.arrhenius[0].A)
        second_quantity = (second.data.kdata if kinetics_type == 'KineticsData'
                           else second.data.A if kinetics_type == 'Arrhenius'
                           else second.data.arrhenius[0].A)
        second_quantity.value_si = first_quantity.value_si
        family = _family([
            _depository('Synthetic/training', [first]),
            _depository('Synthetic/NIST', [second]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    def test_rate_uncertainty_subclass_is_not_assumed_identical(self):
        class TaggedRateUncertainty(RateUncertainty):
            pass

        entries = []
        for tag in ('left', 'right'):
            entry = _entry(1, A_ALPHA)
            uncertainty = TaggedRateUncertainty(
                mu=0.0, var=1.0, Tref=1000.0, N=1, correlation='same')
            uncertainty.tag = tag
            entry.data.uncertainty = uncertainty
            entries.append(entry)
        family = _family([
            _depository('Synthetic/training', [entries[0]]),
            _depository('Synthetic/NIST', [entries[1]]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    def test_quantity_subclass_is_not_assumed_identical(self):
        class TaggedScalarQuantity(ScalarQuantity):
            pass

        entries = []
        for tag in ('left', 'right'):
            entry = _entry(1, A_ALPHA)
            quantity = TaggedScalarQuantity(A_ALPHA, 'cm^3/(mol*s)')
            quantity.tag = tag
            entry.data._A = quantity
            entries.append(entry)
        family = _family([
            _depository('Synthetic/training', [entries[0]]),
            _depository('Synthetic/NIST', [entries[1]]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    @pytest.mark.parametrize('storage_field', ['value_si', 'uncertainty_si'])
    def test_array_subclass_cannot_override_exact_identity(self, storage_field):
        class AlwaysEqualArray(np.ndarray):
            def __array_function__(self, function, types, args, kwargs):
                if function is np.array_equal:
                    return True
                return super().__array_function__(function, types, args, kwargs)

        arrays = [
            np.array([1.0, 2.0]).view(AlwaysEqualArray),
            np.array([100.0, 200.0]).view(AlwaysEqualArray),
        ]
        entries = []
        for array in arrays:
            entry = _entry(1, A_ALPHA)
            entry.data = KineticsData(Tdata=([300.0, 1000.0], 'K'),
                                      kdata=([1.0, 2.0], 's^-1'))
            setattr(entry.data.kdata, storage_field, array)
            entries.append(entry)
        family = _family([
            _depository('Synthetic/training', [entries[0]]),
            _depository('Synthetic/NIST', [entries[1]]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    @pytest.mark.parametrize('storage_field', ['value_si', 'uncertainty_si'])
    def test_array_dtype_metadata_is_not_assumed_identical(self, storage_field):
        arrays = [
            np.array([1.0, 2.0], dtype=np.dtype('float64', metadata={'source': 'left'})),
            np.array([1.0, 2.0], dtype=np.dtype('float64', metadata={'source': 'right'})),
        ]
        entries = []
        for array in arrays:
            entry = _entry(1, A_ALPHA)
            entry.data = KineticsData(Tdata=([300.0, 1000.0], 'K'),
                                      kdata=([1.0, 2.0], 's^-1'))
            setattr(entry.data.kdata, storage_field, array)
            entries.append(entry)
        family = _family([
            _depository('Synthetic/training', [entries[0]]),
            _depository('Synthetic/NIST', [entries[1]]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    @pytest.mark.parametrize('storage_field', ['value_si', 'uncertainty_si'])
    def test_opaque_object_array_is_not_assumed_identical(self, storage_field):
        class OpaqueValue:
            def __init__(self, payload):
                self.payload = payload

            def __eq__(self, other):
                return True

            def __mul__(self, other):
                return self

            __rmul__ = __mul__

            def __truediv__(self, other):
                return self

        arrays = [
            np.array([OpaqueValue('left-1'), OpaqueValue('left-2')], dtype=object),
            np.array([OpaqueValue('right-1'), OpaqueValue('right-2')], dtype=object),
        ]
        entries = []
        for array in arrays:
            entry = _entry(1, A_ALPHA)
            entry.data = KineticsData(Tdata=([300.0, 1000.0], 'K'),
                                      kdata=([1.0, 2.0], 's^-1'))
            setattr(entry.data.kdata, storage_field, array)
            entries.append(entry)
        family = _family([
            _depository('Synthetic/training', [entries[0]]),
            _depository('Synthetic/NIST', [entries[1]]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    @pytest.mark.parametrize('storage_field', ['value_si', 'uncertainty_si'])
    def test_independent_plain_float64_arrays_are_identical(self, storage_field):
        entries = []
        for _ in range(2):
            entry = _entry(1, A_ALPHA)
            entry.data = KineticsData(Tdata=([300.0, 1000.0], 'K'),
                                      kdata=([1.0, 2.0], 's^-1'))
            array = np.array([1.0, 2.0], dtype=np.float64)
            setattr(entry.data.kdata, storage_field, array)
            entries.append(entry)
        assert type(getattr(entries[0].data.kdata, storage_field)) is np.ndarray
        family = _family([
            _depository('Synthetic/training', [entries[0]]),
            _depository('Synthetic/copy-of-training', [entries[1]]),
        ])
        kinetics, _, _, _ = _ask(family)
        assert type(getattr(kinetics.kdata, storage_field)) is np.ndarray

    @pytest.mark.parametrize(('field', 'left', 'right'), [
        ('Tmin', (300.0, 'K'), (301.0, 'K')),
        ('Tmax', (2000.0, 'K'), (2001.0, 'K')),
        ('Pmin', (0.1, 'bar'), (0.2, 'bar')),
        ('Pmax', (10.0, 'bar'), (11.0, 'bar')),
    ])
    def test_validity_bounds_are_part_of_exact_identity(self, field, left, right):
        first = _entry(1, A_ALPHA)
        second = _entry(1, A_ALPHA)
        setattr(first.data, field, left)
        setattr(second.data, field, right)
        family = _family([
            _depository('Synthetic/training', [first]),
            _depository('Synthetic/NIST', [second]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    def test_different_rate_uncertainty_is_not_an_exact_duplicate(self):
        """Uncertainty is part of a source record even when the rate itself is identical."""
        first = _entry(1, A_ALPHA)
        second = _entry(1, A_ALPHA)
        first.data.uncertainty = RateUncertainty(mu=0.0, var=1.0, Tref=1000.0, N=1,
                                                 correlation='same')
        second.data.uncertainty = RateUncertainty(mu=0.0, var=100.0, Tref=1000.0, N=1,
                                                  correlation='same')
        family = _family([
            _depository('Synthetic/training', [first]),
            _depository('Synthetic/NIST', [second]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    def test_identical_kinetics_data_records_are_proven_identical(self):
        """A known tabulated kinetics shape can be proven equal without object identity."""
        entries = []
        for _ in range(2):
            entry = _entry(1, A_ALPHA)
            entry.data = KineticsData(
                Tdata=([300.0, 1000.0], 'K'),
                kdata=([1.0e6, 2.0e6], 'cm^3/(mol*s)'),
                Tmin=(300.0, 'K'),
                Tmax=(2000.0, 'K'),
                comment='same tabulation',
            )
            entries.append(entry)
        family = _family([
            _depository('Synthetic/training', [entries[0]]),
            _depository('Synthetic/copy-of-training', [entries[1]]),
        ])
        kinetics, _, _, _ = _ask(family)
        assert kinetics.kdata.value_si.tolist() == entries[0].data.kdata.value_si.tolist()

    def test_kinetics_data_uncertainty_is_part_of_exact_identity(self):
        entries = []
        for variance in (1.0, 100.0):
            entry = _entry(1, A_ALPHA)
            entry.data = KineticsData(
                Tdata=([300.0, 1000.0], 'K'),
                kdata=([1.0e6, 2.0e6], 'cm^3/(mol*s)'),
            )
            entry.data.uncertainty = RateUncertainty(
                mu=0.0, var=variance, Tref=1000.0, N=1, correlation='same')
            entries.append(entry)
        family = _family([
            _depository('Synthetic/training', [entries[0]]),
            _depository('Synthetic/NIST', [entries[1]]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    def test_identical_pdep_arrhenius_records_are_proven_identical(self):
        entries = []
        for _ in range(2):
            entry = _entry(1, A_ALPHA)
            entry.data = PDepArrhenius(
                pressures=([0.1, 1.0], 'bar'),
                arrhenius=[
                    Arrhenius(A=(1.0e6, 'cm^3/(mol*s)'), n=0.0, Ea=(10.0, 'kJ/mol')),
                    Arrhenius(A=(2.0e6, 'cm^3/(mol*s)'), n=0.0, Ea=(20.0, 'kJ/mol')),
                ],
                Tmin=(300.0, 'K'),
                Tmax=(2000.0, 'K'),
                comment='same pressure dependence',
            )
            entries.append(entry)
        family = _family([
            _depository('Synthetic/training', [entries[0]]),
            _depository('Synthetic/copy-of-training', [entries[1]]),
        ])
        kinetics, _, _, _ = _ask(family)
        assert kinetics.pressures.value_si.tolist() == entries[0].data.pressures.value_si.tolist()

    def test_identical_sticking_coefficient_records_are_proven_identical(self):
        entries = []
        for _ in range(2):
            entry = _entry(1, A_ALPHA)
            entry.data = StickingCoefficient(
                A=0.25,
                n=0.5,
                Ea=(5.0, 'kJ/mol'),
                T0=(1.0, 'K'),
                Tmin=(300.0, 'K'),
                Tmax=(2000.0, 'K'),
                comment='same sticking law',
            )
            entries.append(entry)
        family = _family([
            _depository('Synthetic/training', [entries[0]]),
            _depository('Synthetic/copy-of-training', [entries[1]]),
        ])
        kinetics, _, _, _ = _ask(family)
        assert kinetics.A.value_si == entries[0].data.A.value_si

    def test_identical_multi_arrhenius_records_are_proven_identical(self):
        entries = []
        for _ in range(2):
            entry = _entry(1, A_ALPHA)
            entry.data = MultiArrhenius(
                arrhenius=[
                    Arrhenius(A=(1.0e6, 'cm^3/(mol*s)'), n=0.0, Ea=(10.0, 'kJ/mol')),
                    Arrhenius(A=(2.0e6, 'cm^3/(mol*s)'), n=0.0, Ea=(20.0, 'kJ/mol')),
                ],
                Tmin=(300.0, 'K'),
                Tmax=(2000.0, 'K'),
                comment='same summed law',
            )
            entries.append(entry)
        family = _family([
            _depository('Synthetic/training', [entries[0]]),
            _depository('Synthetic/copy-of-training', [entries[1]]),
        ])
        kinetics, _, _, _ = _ask(family)
        assert len(kinetics.arrhenius) == 2

    @pytest.mark.parametrize('kinetics_type', [
        'PDepArrhenius',
        'StickingCoefficient',
        'MultiArrhenius',
    ])
    def test_inherited_uncertainty_is_part_of_exact_identity(self, kinetics_type):
        def make_kinetics():
            arrhenius = [Arrhenius(A=(1.0e6, 'cm^3/(mol*s)'), n=0.0,
                                   Ea=(10.0, 'kJ/mol'))]
            if kinetics_type == 'PDepArrhenius':
                return PDepArrhenius(pressures=([1.0], 'bar'), arrhenius=arrhenius)
            if kinetics_type == 'StickingCoefficient':
                return StickingCoefficient(A=0.25, n=0.5, Ea=(5.0, 'kJ/mol'))
            return MultiArrhenius(arrhenius=arrhenius)

        entries = []
        for variance in (1.0, 100.0):
            entry = _entry(1, A_ALPHA)
            entry.data = make_kinetics()
            entry.data.uncertainty = RateUncertainty(
                mu=0.0, var=variance, Tref=1000.0, N=1, correlation='same')
            entries.append(entry)
        family = _family([
            _depository('Synthetic/training', [entries[0]]),
            _depository('Synthetic/NIST', [entries[1]]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    @pytest.mark.parametrize('kinetics_type', ['PDepArrhenius', 'MultiArrhenius'])
    def test_nested_subkinetics_uncertainty_is_part_of_exact_identity(self, kinetics_type):
        entries = []
        for variance in (1.0, 100.0):
            nested = Arrhenius(A=(1.0e6, 'cm^3/(mol*s)'), n=0.0, Ea=(10.0, 'kJ/mol'))
            nested.uncertainty = RateUncertainty(
                mu=0.0, var=variance, Tref=1000.0, N=1, correlation='same')
            entry = _entry(1, A_ALPHA)
            entry.data = (PDepArrhenius(pressures=([1.0], 'bar'), arrhenius=[nested])
                          if kinetics_type == 'PDepArrhenius'
                          else MultiArrhenius(arrhenius=[nested]))
            entries.append(entry)
        family = _family([
            _depository('Synthetic/training', [entries[0]]),
            _depository('Synthetic/NIST', [entries[1]]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    def test_reaction_degeneracy_is_part_of_exact_identity(self):
        first = _entry(1, A_ALPHA)
        second = _entry(1, A_ALPHA)
        first.item.degeneracy = 1
        second.item.degeneracy = 2
        family = _family([
            _depository('Synthetic/training', [first]),
            _depository('Synthetic/NIST', [second]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    def test_duplicate_flag_is_part_of_exact_identity(self):
        first = _entry(1, A_ALPHA)
        second = _entry(1, A_ALPHA)
        first.item.duplicate = False
        second.item.duplicate = True
        family = _family([
            _depository('Synthetic/training', [first]),
            _depository('Synthetic/NIST', [second]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    def test_unreviewed_kinetics_class_fails_closed(self):
        entries = []
        for _ in range(2):
            entry = _entry(1, A_ALPHA)
            entry.data = ArrheniusEP(A=(1.0e6, 'cm^3/(mol*s)'), n=0.0,
                                     alpha=0.5, E0=(10.0, 'kJ/mol'))
            entries.append(entry)
        family = _family([
            _depository('Synthetic/training', [entries[0]]),
            _depository('Synthetic/NIST', [entries[1]]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    def test_opposite_one_way_entry_cannot_hide_forward_collision_by_rank(self):
        """Normalize direction before selecting an entry within a depository."""
        family = _family([
            _depository('Synthetic/training', [
                _entry(1, A_ALPHA, reversible=False, rank=5),
                _entry(2, A_BETA, reversible=False, flip=True, rank=1),
            ]),
            _depository('Synthetic/NIST', [_entry(7, A_BETA, reversible=False, rank=5)]),
        ])
        with pytest.raises(_conflict_error()) as caught:
            _ask(family)
        assert 'index 1' in str(caught.value) and 'index 7' in str(caught.value)

    def test_numerically_equal_but_independently_sourced_entries_are_still_two_sources(self):
        """
        §3, second sentence.  Two entries that happen to agree at one temperature remain
        two sources.  Same numbers, different citation and different validity range: still
        a conflict.
        """
        family = _family([
            _depository('Synthetic/training', [_entry(1, A_ALPHA, reference='Alpha, 1990')]),
            _depository('Synthetic/NIST', [_entry(7, A_ALPHA, reference='Beta, 2004', tmax=1200.0)]),
        ])
        with pytest.raises(_conflict_error()):
            _ask(family)

    def test_the_report_survives_a_citation_it_cannot_print(self):
        """
        An exception's message is the last thing between a user and a run they cannot
        explain, so it must not be able to fail.

        Not hypothetical: 59 entries in the shipped database -- all of
        ``2+2_cycloaddition/NIST`` and part of ``H_Abstraction/training`` -- carry byte
        strings in their reference author lists, and ``str()`` on one of those raises
        ``TypeError`` inside ``rmgpy/data/reference.py``. Before this guard, a depository
        conflict on any of those reactions surfaced as a string-formatting error rather
        than as the conflict it was. The data defect is real and is reported in
        ``docs/depository-order/census.md``; this test pins that the reporter survives it
        instead of hiding behind it.
        """
        class Unprintable:
            def __str__(self):
                raise TypeError('sequence item 0: expected str instance, bytes found')

            def __repr__(self):
                raise TypeError('and repr is no better')

        hostile = _entry(7, A_BETA)
        hostile.reference = Unprintable()
        family = _family([
            _depository('Synthetic/training', [_entry(1, A_ALPHA)]),
            _depository('Synthetic/NIST', [hostile]),
        ])

        with pytest.raises(_conflict_error()) as caught:
            _ask(family)
        message = str(caught.value)
        assert 'unprintable' in message
        # ... and the rest of the report is still there, so one bad field costs one field.
        assert 'Synthetic/training' in message and 'Synthetic/NIST' in message
        assert 'index 1' in message and 'index 7' in message

    def test_two_unprintable_citations_are_not_assumed_identical(self):
        """
        The other half of the same data defect.  The duplicate test also reads the
        citation, and a citation nobody can render has proven nothing -- so two of them
        must not be taken as proof that two entries are one source.  Erring towards
        reporting a conflict is the safe direction; erring the other way deletes a source
        silently.
        """
        class Unprintable:
            def __str__(self):
                raise TypeError('nope')

            def __repr__(self):
                raise TypeError('also nope')

        # Identical in every field the duplicate test looks at EXCEPT the citation, so the
        # citation is the only thing that can decide the verdict. Without this the entries
        # differ in short_desc, the conflict is raised for that reason instead, and the
        # test passes whatever the citation fallback does -- a check that cannot fail.
        first, second = _entry(1, A_ALPHA), _entry(7, A_ALPHA)
        second.short_desc = first.short_desc
        first.reference = Unprintable()
        second.reference = Unprintable()
        family = _family([
            _depository('Synthetic/training', [first]),
            _depository('Synthetic/NIST', [second]),
        ])

        with pytest.raises(_conflict_error()):
            _ask(family)

    def test_non_overlapping_depositories_are_unchanged(self):
        """§6.10.  A second depository that does not carry this reaction changes nothing."""
        other = _entry(4, A_BETA)
        other.item = Reaction(reactants=[Species(smiles='[OH]'), Species(smiles='C')],
                              products=[Species(smiles='O'), Species(smiles='[CH3]')])
        family = _family([
            _depository('Synthetic/training', [_entry(1, A_ALPHA)]),
            _depository('Synthetic/NIST', [other]),
        ])
        kinetics, depository, entry, is_forward = _ask(family)
        assert entry.index == 1
        assert depository.label == 'Synthetic/training'
        assert is_forward is True
        assert kinetics.A.value_si == pytest.approx(
            Arrhenius(A=(A_ALPHA, 'cm^3/(mol*s)')).A.value_si)
        assert 'Matched reaction 1' in kinetics.comment
        assert 'Synthetic/training' in kinetics.comment

    def test_the_error_is_deterministic_and_carries_enough_provenance_to_decide(self):
        """
        §6.11 and §1's field list.  Every field the ruling requires is in the message, and
        two runs of the same conflict produce the same message character for character.
        """
        def build():
            return _family([
                _depository('Synthetic/training', [_entry(1, A_ALPHA)]),
                _depository('Synthetic/NIST', [_entry(7, A_BETA, rank=None, reference='NIST, 1994')]),
            ])

        messages = []
        for _ in range(2):
            with pytest.raises(_conflict_error()) as caught:
                _ask(build())
            messages.append(str(caught.value))
        assert messages[0] == messages[1], 'the conflict report must be deterministic'

        message = messages[0]
        # the canonical oriented reaction
        assert 'H2' in message and 'CH4' in message
        # electron metadata and directionality
        assert 'electron' in message.lower()
        assert 'direction' in message.lower()
        # depository name, entry index and label
        assert 'Synthetic/training' in message and 'index 1' in message
        assert 'Synthetic/NIST' in message and 'index 7' in message
        # kinetics class
        assert 'Arrhenius' in message
        # source/provenance citation
        assert 'NIST, 1994' in message
        assert 'Somebody et al.' in message
        # validity ranges
        assert '300' in message and '2000' in message
        # uncertainty / rank metadata
        assert 'rank' in message.lower()
        # and it tells the user what to do about it
        assert 'kineticsDepositories' in message


################################################################################
# The same policy against the real database on disk.
################################################################################

#: Pinned from the RMG-database checkout ``rmgrc`` points at.  These numbers are the
#: content the labels ``training`` and ``NIST`` are supposed to be pointing at; asserting
#: them is what stops this file becoming a test that passes when the two directories'
#: contents are swapped.  If the database moves and these fail, read the new entries and
#: update them deliberately -- do not relax the assertion.
DATABASE_PROVENANCE = {
    'family': 'HO2_Elimination_from_PeroxyRadical',
    #: C3H7O2 <=> C3H6 + HO2, carried by both depositories in the same direction.  Written
    #: out as SMILES rather than read back out of one of the depositories, so that the
    #: query does not come from the thing under test.
    'reactants': ['CCCO[O]'],
    'products': ['C=CC', '[O]O'],
    'training_index': 3,
    'training_k1000': 8.632e5,
    'nist_index': 4,
    'nist_k1000': 6.644e4,
    'nist_author': 'DeSain',
}


def _load_family(depositories):
    database = KineticsDatabase()
    database.load_families(
        os.path.join(settings['database.directory'], 'kinetics', 'families'),
        families=[DATABASE_PROVENANCE['family']],
        depositories=depositories,
    )
    return database.families[DATABASE_PROVENANCE['family']]


def _contested_query():
    """The oriented reaction both shipped depositories claim, stated rather than derived."""
    return Reaction(
        reactants=[Species(smiles=s) for s in DATABASE_PROVENANCE['reactants']],
        products=[Species(smiles=s) for s in DATABASE_PROVENANCE['products']],
    )


def _ask_real(family, reaction):
    return family.get_kinetics(
        reaction,
        template_labels=[g.label for g in family.forward_template.reactants],
        degeneracy=1,
        estimator='',
        return_all_kinetics=False,
    )


@pytest.mark.database
class TestThePolicyAgainstTheShippedDatabase:
    """§6.4-§6.6 and §6.12 -- the same rules, against real depositories."""

    def test_all_does_not_select_a_winner(self):
        """
        §6.4.  ``kineticsDepositories = 'all'`` is a statement about which depositories are
        included, not permission to pick one of them.
        """
        family = _load_family('all')
        assert len(family.depositories) == 2, 'fixture family must really carry two depositories'
        with pytest.raises(_conflict_error()) as caught:
            _ask_real(family, _contested_query())
        message = str(caught.value)
        assert '/training' in message and '/NIST' in message

    def test_selecting_only_training_uses_the_training_entry(self):
        """§6.5.  Named unambiguously, the run proceeds -- with training's number."""
        family = _load_family(['training'])
        assert [d.label for d in family.depositories] == \
            ['{0}/training'.format(DATABASE_PROVENANCE['family'])]
        kinetics, depository, entry, is_forward = _ask_real(family, _contested_query())
        assert entry.index == DATABASE_PROVENANCE['training_index']
        assert kinetics.get_rate_coefficient(1000.0) == \
            pytest.approx(DATABASE_PROVENANCE['training_k1000'], rel=1e-3)

    def test_selecting_only_nist_uses_the_nist_entry(self):
        """
        §6.6.  The mirror image, and the place the input semantics bite: a bare
        ``['NIST']`` is silently turned into ``['NIST', 'training']`` by
        ``KineticsFamily.load``, so selecting NIST alone requires the documented
        ``'!training'`` flag.  The rate that comes back is NIST's, and it is a different
        number from training's -- which is what makes this an assertion about content and
        not about a label.
        """
        family = _load_family(['NIST', '!training'])
        assert [d.label for d in family.depositories] == \
            ['{0}/NIST'.format(DATABASE_PROVENANCE['family'])]
        kinetics, depository, entry, is_forward = _ask_real(family, _contested_query())
        assert entry.index == DATABASE_PROVENANCE['nist_index']
        assert kinetics.get_rate_coefficient(1000.0) == \
            pytest.approx(DATABASE_PROVENANCE['nist_k1000'], rel=1e-3)
        assert DATABASE_PROVENANCE['nist_author'] in str(entry.reference)
        assert DATABASE_PROVENANCE['nist_k1000'] != DATABASE_PROVENANCE['training_k1000']
        assert depository.label == '{0}/NIST'.format(DATABASE_PROVENANCE['family'])
        assert is_forward is True
        assert 'Matched reaction {0}'.format(DATABASE_PROVENANCE['nist_index']) in kinetics.comment
        assert depository.label in kinetics.comment

    def test_two_sided_ionization_entry_matches_the_generated_representation(self, tmp_path):
        """A loaded family depository keeps the owner needed for electron identity."""
        family_label = 'Plasma_Electron_Impact_Ionization'
        root = tmp_path / family_label / 'training'
        root.mkdir(parents=True)
        (root / 'dictionary.txt').write_text("""Ar
1 Ar u0 p4 c0

Arp
multiplicity 2
1 Ar u1 p3 c+1
""")
        (root / 'reactions.py').write_text("""entry(
    index = 1,
    label = 'Ar => Arp',
    reversible = False,
    kinetics = Arrhenius(A=(2.5e12, 's^-1'), n=0, Ea=(0, 'kJ/mol'), T0=(1, 'K')),
    rank = 5,
    shortDesc = 'ionization fixture',
)
""")
        depository = KineticsDepository(label=family_label + '/training', electrons=1)
        depository.load(str(root / 'reactions.py'), KineticsDatabase().local_context, {})
        loaded = depository.entries[1]
        generated = TemplateReaction(
            reactants=list(loaded.item.reactants),
            products=list(loaded.item.products),
            family=family_label,
            electrons=1,
            reversible=False,
        )

        assert get_electron_placement_counts(loaded.item) == (1, 2)
        assert get_electron_placement_counts(generated) == (1, 2)
        assert loaded.item.is_same_reaction(generated, either_direction=False)

        family = _family([depository])
        family.label = family_label
        family.electrons = 1
        kinetics, source, entry, is_forward = _ask(family, generated)
        assert source is depository
        assert entry is loaded
        assert is_forward is True
        assert kinetics.A.value_si == pytest.approx(2.5e12)

    #: label -> the depositories the family must load under ``'all'``, and how many entries
    #: each must carry.  The electrochemical fixture is deliberately non-empty so loading
    #: only a directory shell cannot satisfy §6.12; naming every count also catches a
    #: depository silently dropped.
    LOADING_PATHS = {
        'H_Abstraction': {'training': 3117, 'NIST': 1312},          # ordinary thermal
        'Plasma_Electron_Attachment': {'training': 3},              # plasma, electron metadata
        'Surface_Abstraction': {'training': 4},                     # surface
        'Surface_Proton_Electron_Reduction_Alpha': {'training': 6},  # electrochemical
    }

    @pytest.mark.parametrize('family_label', sorted(LOADING_PATHS))
    def test_ordinary_loading_paths_remain_green(self, family_label):
        """
        §6.12.  Thermal, plasma, surface and electrochemical families still load under
        ``'all'``, with exactly the depositories and entry counts they had before.
        """
        expected = self.LOADING_PATHS[family_label]
        database = KineticsDatabase()
        database.load_families(
            os.path.join(settings['database.directory'], 'kinetics', 'families'),
            families=[family_label],
            depositories='all',
        )
        family = database.families[family_label]
        loaded = {d.label.split('/')[-1]: len(d.entries) for d in family.depositories}
        assert loaded == expected


@pytest.mark.database
class TestProvenanceExtractionIsNotTrainingOnly:
    """
    A rate taken from a depository other than ``training`` must still be traceable.

    ``get_kinetics_from_depository`` writes "Matched reaction <index> <label> in
    <family>/<depository>" into the kinetics comment for every depository.
    ``extract_source_from_comments`` used to look for that line with a pattern ending in
    ``.*training``, so a rate from any other depository missed it, fell through to the
    rate-rule parse, and died on "Could not find rate rule in comments" -- a reaction whose
    provenance the comment states in full, reported as having none.

    This is also the retraction of a claim made while this branch was being written: that
    a model carries no record of which depository a rate came from. It does. The comment
    above is written into ``chem_annotated.inp`` verbatim; what was broken was the reader,
    not the record.
    """

    @staticmethod
    def _matched(depositories, index):
        database = KineticsDatabase()
        database.load_families(
            os.path.join(settings['database.directory'], 'kinetics', 'families'),
            families=[DATABASE_PROVENANCE['family']],
            depositories=depositories,
        )
        family = database.families[DATABASE_PROVENANCE['family']]
        depository = family.depositories[0]
        entry = depository.entries[index]
        reaction = TemplateReaction(reactants=list(entry.item.reactants),
                                    products=list(entry.item.products),
                                    family=DATABASE_PROVENANCE['family'])
        kinetics, _, _, _ = family.get_kinetics(
            reaction,
            template_labels=[g.label for g in family.forward_template.reactants],
            degeneracy=1, estimator='', return_all_kinetics=False)
        reaction.kinetics = kinetics
        return family, depository, reaction

    @pytest.mark.parametrize('depositories, index, expected_k1000', [
        (['training'], DATABASE_PROVENANCE['training_index'], DATABASE_PROVENANCE['training_k1000']),
        (['NIST', '!training'], DATABASE_PROVENANCE['nist_index'], DATABASE_PROVENANCE['nist_k1000']),
    ], ids=['training', 'NIST'])
    def test_the_source_is_recovered_from_the_comment(self, depositories, index, expected_k1000):
        family, depository, reaction = self._matched(depositories, index)

        assert 'Matched reaction {0}'.format(index) in reaction.kinetics.comment
        assert depository.label in reaction.kinetics.comment

        from_depository, source = family.extract_source_from_comments(reaction)

        assert from_depository is True
        assert source[0] == DATABASE_PROVENANCE['family']
        assert source[1] is depository.entries[index]
        # The recovered entry is the one that supplied the number, not merely one with the
        # right index in some depository.
        assert source[1].data.get_rate_coefficient(1000.0) == pytest.approx(expected_k1000, rel=1e-3)
        assert source[2] is False  # matched in the forward direction

    def test_a_depository_named_in_the_comment_but_not_loaded_says_so(self):
        """
        Reading a mechanism back with a different set of depositories loaded is a real
        mismatch, and silently attributing the rate to whatever happens to be loaded would
        put one depository's number under another's name.
        """
        family, _, reaction = self._matched(['NIST', '!training'],
                                            DATABASE_PROVENANCE['nist_index'])
        family.depositories = []

        with pytest.raises(DatabaseError) as caught:
            family.extract_source_from_comments(reaction)
        assert 'NIST' in str(caught.value)
        assert 'not loaded' in str(caught.value)
