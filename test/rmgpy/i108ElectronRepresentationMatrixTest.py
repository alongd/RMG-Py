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
The electron-representation matrix: which lithium plasma reaction shapes survive
which of the three gates.

Two prior investigations disagreed about whether ``Reaction.is_balanced``
treating the free electron as a conserved pseudo-element is a defect or a
load-bearing invariant, and each measured a different base. This module is the
re-measurement, pinned as tests so the answer stops being a transcription.

Eight reaction shapes are put through three INDEPENDENT gates:

1. ``Reaction.is_balanced()`` -- a pure function of the participants plus the
   signed scalar ``Reaction.electrons``.
2. The kinetics-library loader, ``KineticsLibrary.load`` -- which runs
   ``is_balanced()`` itself, but only AFTER copying ``entry.data.electrons``
   onto the reaction when the rate law carries that attribute
   (``rmgpy/data/kinetics/library.py``). The escape is therefore borne by the
   RATE LAW, not by the ``entry()``/``load_entry`` signature.
3. :class:`rmgpy.solver.plasma.PlasmaReactor`, via ``initialize_model`` --
   which refuses any nonzero metadata electron count outright, and routes it
   first through :func:`rmgpy.electron_placement.resolve_electron_placement`,
   whose ``FAMILY_ELECTRON_PLACEMENT`` declares only two attachment-shaped
   families.

The kinetics used for the ionisation and recombination rows are the real
lithium fits shipped in the plasma database -- ``VoronovEIArrhenius(Z=3, N=3)``
reads ``input/kinetics/voronov.yaml`` and ``BadnellRRArrhenius(Z=3, N=2)``
reads ``input/kinetics/badnell.yaml``, both resolved through
``settings['database.directory']`` -- so gate 2 is answered against genuine
electron-impact-ionisation and radiative-recombination data rather than a
placeholder.

Run ``python test/rmgpy/i108ElectronRepresentationMatrixTest.py`` to print the
measured table; run it under pytest to assert it.
"""

import os
import shutil
import tempfile

import pytest

from rmgpy import settings
from rmgpy.data.kinetics.library import KineticsLibrary
from rmgpy.exceptions import DatabaseError, ElectronPlacementError
from rmgpy.kinetics.arrhenius import (Arrhenius, BadnellRRArrhenius, TwoTemperaturePlasma,
                                      VoronovEIArrhenius)
from rmgpy.reaction import Reaction
from rmgpy.solver.plasma import PlasmaReactor
from rmgpy.species import Species

################################################################################

#: Adjacency lists for every participant. ``Li`` is the neutral atom (one
#: unpaired valence electron), ``Liplus`` the closed-shell cation; both spellings
#: match the plasma database's own families (e.g. ``Li_Abstraction/groups.py``
#: and ``Cation_Li_Abstraction/groups.py``).
ADJACENCY_LISTS = {
    'e': '1 e u1 p0 c-1',
    'Li': '1 Li u1 p0 c0',
    'Liplus': '1 Li u0 p0 c+1',
    'Ar': '1 Ar u0 p4 c0',
    'Arplus': '1 Ar u1 p3 c+1',
}

T_GAS = 1000.0
T_E = 11000.0
P0 = 1.0e5
Y_E0 = 1.0e-6


def _species():
    """A fresh label -> Species map. Fresh per row: the reactor indexes species
    by identity, and a shape that appears twice must not share objects."""
    return {label: Species(label=label).from_adjacency_list(adjlist)
            for label, adjlist in ADJACENCY_LISTS.items()}


def _voronov_li(electrons):
    """Real electron-impact ionisation kinetics for Li I -> Li II, read from the
    plasma database's ``voronov.yaml`` (Z=3, N=3)."""
    return VoronovEIArrhenius(Z=3, N=3, electrons=electrons)


def _badnell_li(electrons):
    """Real radiative recombination kinetics for Li II -> Li I, read from the
    plasma database's ``badnell.yaml`` (Z=3, N=2, i.e. the two electrons Li+
    carries before capture)."""
    return BadnellRRArrhenius(Z=3, N=2, electrons=electrons)


def _three_body():
    """A well-formed third-order plasma rate law. Carries no ``electrons``
    attribute of its own, which is exactly what makes it a control for the
    rate-law-borne metadata escape."""
    return TwoTemperaturePlasma(A=(1.0e4, 'm^6/(mol^2*s)'), n=-1.0,
                                Ea_g=(0.0, 'J/mol'), Ea_e=(0.0, 'J/mol'))


def _thermal_bimolecular():
    """An ordinary Arrhenius rate law: no electron temperature, no electron
    density, no ``electrons`` attribute."""
    return Arrhenius(A=(1.0e12, 'cm^3/(mol*s)'), n=0.0, Ea=(0.0, 'kJ/mol'))


################################################################################

#: The eight rows. ``electrons`` is the signed scalar the shape is *written*
#: with; ``kinetics`` is a zero-argument factory for the chemically appropriate
#: rate law, whose own ``electrons`` field (where it has one) is set to the same
#: value, since that is the only channel by which the library route can deliver
#: the scalar.
MATRIX = [
    {
        'family': 'EI', 'form': 'explicit',
        'equation': 'Li + e => Li+ + e + e',
        'reactants': ['Li', 'e'], 'products': ['Liplus', 'e', 'e'],
        'electrons': 0,
        'kinetics': lambda: _voronov_li(0),
    },
    {
        'family': 'EI', 'form': 'half',
        'equation': 'Li => Li+ + e',
        'reactants': ['Li'], 'products': ['Liplus', 'e'],
        'electrons': 0,
        'kinetics': lambda: _voronov_li(0),
    },
    {
        'family': 'EI', 'form': 'metadata',
        'equation': 'Li => Li+  (electrons=+1)',
        'reactants': ['Li'], 'products': ['Liplus'],
        'electrons': 1,
        'kinetics': lambda: _voronov_li(1),
    },
    {
        'family': 'RR', 'form': 'explicit',
        'equation': 'Li+ + e => Li',
        'reactants': ['Liplus', 'e'], 'products': ['Li'],
        'electrons': 0,
        'kinetics': lambda: _badnell_li(0),
    },
    {
        'family': 'RR', 'form': 'metadata',
        'equation': 'Li+ => Li  (electrons=-1)',
        'reactants': ['Liplus'], 'products': ['Li'],
        'electrons': -1,
        'kinetics': lambda: _badnell_li(-1),
    },
    {
        'family': '3b', 'form': 'e-e-ion',
        'equation': 'Li+ + e + e => Li + e',
        'reactants': ['Liplus', 'e', 'e'], 'products': ['Li', 'e'],
        'electrons': 0,
        'kinetics': _three_body,
    },
    {
        'family': '3b', 'form': 'neutral M',
        'equation': 'Li+ + e + Ar => Li + Ar',
        'reactants': ['Liplus', 'e', 'Ar'], 'products': ['Li', 'Ar'],
        'electrons': 0,
        'kinetics': _three_body,
    },
    {
        'family': 'CT', 'form': 'ion-atom',
        'equation': 'Li+ + Ar => Li + Ar',
        'reactants': ['Liplus', 'Ar'], 'products': ['Li', 'Ar'],
        'electrons': 0,
        'kinetics': _thermal_bimolecular,
    },
]


#: Rows the eight-row matrix does not contain, kept separate so the matrix stays
#: exactly the shape it was asked for.
#:
#: ``CT/charge-balanced`` exists because the matrix's own ``CT ion-atom`` row,
#: ``Li+ + Ar => Li + Ar``, does not conserve charge (+1 on the left, 0 on the
#: right) and has no electron anywhere in it. It is refused for a reason that has
#: nothing to do with the electron representation, and the shape that DOES
#: describe ion-atom charge transfer is written here so the matrix's negative
#: result is not mistaken for one.
SUPPLEMENTARY = [
    {
        'family': 'CT', 'form': 'charge-balanced',
        'equation': 'Li+ + Ar => Li + Ar+',
        'reactants': ['Liplus', 'Ar'], 'products': ['Li', 'Arplus'],
        'electrons': 0,
        'kinetics': _thermal_bimolecular,
    },
    {
        # The form ``check_electron_reactant_order``'s own docstring calls the
        # correct RMG representation of electron-impact ionisation -- "put the
        # consumed electron in reaction.reactants and count only surplus produced
        # electrons in reaction.electrons" -- and the form the shipped
        # plasma-local-context fixture writes (entry 4, ``e + H => proton + e``).
        'family': 'EI', 'form': 'half-explicit',
        'equation': 'Li + e => Li+ + e  (electrons=+1)',
        'reactants': ['Li', 'e'], 'products': ['Liplus', 'e'],
        'electrons': 1,
        'kinetics': lambda: _voronov_li(1),
    },
    {
        'family': '3b', 'form': 'e-e-ion+meta',
        'equation': 'Li+ + e + e => Li + e  (electrons=-1)',
        'reactants': ['Liplus', 'e', 'e'], 'products': ['Li', 'e'],
        'electrons': -1,
        'kinetics': lambda: _badnell_li(-1),
    },
]


################################################################################
# Gate 1 -- Reaction.is_balanced()
################################################################################

def _element_counts(species_list):
    """Per-element atom counts over a participant list, counted the way
    ``Reaction.is_balanced`` counts them -- which includes the free-electron
    pseudo-species, since ``element_list`` begins with element ``e``."""
    counts = {}
    charge = 0
    for spc in species_list:
        for atom in spc.molecule[0].atoms:
            counts[atom.element.symbol] = counts.get(atom.element.symbol, 0) + 1
            charge += atom.charge
    return counts, charge


def _balance_diagnosis(row):
    """Say WHICH of ``is_balanced``'s two comparisons a row fails.

    The distinction is the entire subject of the disagreement this module
    re-measures: an ``E``-element mismatch means the shape is structurally
    unrepresentable with an explicit electron and no scalar can rescue it,
    whereas a charge mismatch means the chemistry itself does not conserve
    charge.
    """
    spc = _species()
    r_counts, r_charge = _element_counts([spc[label] for label in row['reactants']])
    p_counts, p_charge = _element_counts([spc[label] for label in row['products']])
    electrons = row['electrons']
    if electrons < 0:
        r_charge += electrons
    elif electrons > 0:
        p_charge -= electrons

    element_mismatch = sorted(
        sym for sym in set(r_counts) | set(p_counts)
        if r_counts.get(sym, 0) != p_counts.get(sym, 0))
    parts = []
    if element_mismatch:
        parts.append('element mismatch on {0} (L {1} / R {2})'.format(
            ','.join(element_mismatch),
            {s: r_counts.get(s, 0) for s in element_mismatch},
            {s: p_counts.get(s, 0) for s in element_mismatch}))
    if r_charge != p_charge:
        parts.append('charge mismatch (L {0:+d} / R {1:+d} after folding '
                     'electrons={2:+d})'.format(r_charge, p_charge, electrons))
    return '; '.join(parts) if parts else 'balanced'


def gate_is_balanced(row):
    """Answer ``(accepted, detail)`` for ``Reaction.is_balanced()``."""
    spc = _species()
    rxn = Reaction(
        reactants=[spc[label] for label in row['reactants']],
        products=[spc[label] for label in row['products']],
        electrons=row['electrons'],
        reversible=False,
        kinetics=row['kinetics'](),
    )
    return bool(rxn.is_balanced()), 'electrons={0:+d}; {1}'.format(
        row['electrons'], _balance_diagnosis(row))


################################################################################
# Gate 2 -- the kinetics-library loader
################################################################################

#: Rendered into the temporary library's ``reactions.py``. The loader ``exec``s
#: that file against ``KineticsDatabase.local_context``, so the rate law is
#: reconstructed from its ``repr`` -- which is also what proves the ``electrons``
#: field round-trips through the persisted form.
LIBRARY_TEMPLATE = '''#!/usr/bin/env python
# encoding: utf-8

name = "{name}"
shortDesc = u"I-108 electron-representation matrix probe"
longDesc = u"""Generated by test/rmgpy/i108ElectronRepresentationMatrixTest.py."""

entry(
    index = 1,
    label = "{label}",
    reversible = False,
    kinetics = {kinetics!r},
    shortDesc = u"""probe""",
)
'''


def _library_label(row):
    """The ``a + b => c`` string the loader parses. Species names must match the
    dictionary keys, so the scalar electron count cannot appear here -- which is
    the whole point of the gate."""
    return ' + '.join(row['reactants']) + ' => ' + ' + '.join(row['products'])


def gate_library_loader(row, kinetics=None, name='i108-probe'):
    """Write a one-entry kinetics library for ``row`` into a temporary
    directory, load it, and answer ``(accepted, detail)``.

    ``detail`` carries the ``electrons`` value the loader actually stored on the
    reaction when it succeeded, and the loader's own refusal message when it did
    not.
    """
    kinetics = kinetics if kinetics is not None else row['kinetics']()
    tempdir = tempfile.mkdtemp(prefix='i108-')
    try:
        libdir = os.path.join(tempdir, name)
        os.makedirs(libdir)
        with open(os.path.join(libdir, 'dictionary.txt'), 'w') as f:
            f.write('\n\n'.join('{0}\n{1}'.format(label, adjlist)
                                for label, adjlist in ADJACENCY_LISTS.items()) + '\n')
        with open(os.path.join(libdir, 'reactions.py'), 'w') as f:
            f.write(LIBRARY_TEMPLATE.format(
                name=name, label=_library_label(row), kinetics=kinetics))

        library = KineticsLibrary(label=name)
        try:
            library.load(os.path.join(libdir, 'reactions.py'),
                         local_context=_local_context(), global_context={})
        except DatabaseError as exc:
            return False, str(exc).strip()
        rxn = list(library.entries.values())[0].item
        return True, 'stored electrons={0:+d} on {1}'.format(
            rxn.electrons, rxn.__class__.__name__)
    finally:
        shutil.rmtree(tempdir, ignore_errors=True)


def _local_context():
    """The subset of ``KineticsDatabase.local_context`` these entries name."""
    return {
        'Arrhenius': Arrhenius,
        'BadnellRRArrhenius': BadnellRRArrhenius,
        'VoronovEIArrhenius': VoronovEIArrhenius,
        'TwoTemperaturePlasma': TwoTemperaturePlasma,
    }


################################################################################
# Gate 3 -- PlasmaReactor.initialize_model
################################################################################

def gate_reactor(row):
    """Hand ``row`` to a :class:`PlasmaReactor` as a core reaction and answer
    ``(accepted, detail)``.

    The reaction is a bare :class:`Reaction` with no family attribution, which
    is what a hand-written or library-sourced reaction is at this boundary: a
    ``LibraryReaction`` sets ``family`` to the LIBRARY's label, and no library
    label appears in ``FAMILY_ELECTRON_PLACEMENT`` either, so both land on a
    named refusal for metadata-carrying rows.
    """
    spc = _species()
    rxn = Reaction(
        reactants=[spc[label] for label in row['reactants']],
        products=[spc[label] for label in row['products']],
        electrons=row['electrons'],
        reversible=False,
        kinetics=row['kinetics'](),
    )
    core_species = [spc['e']] + [spc[label] for label in
                                 dict.fromkeys(row['reactants'] + row['products'])
                                 if label != 'e']
    imf = {spc['e']: Y_E0}
    others = [s for s in core_species if not s.is_electron()]
    for s in others:
        imf[s] = (1.0 - Y_E0) / len(others)
    reactor = PlasmaReactor(T_GAS, P0, imf, (T_E, 'K'), n_sims=1, termination=[])
    try:
        reactor.initialize_model(core_species, [rxn], [], [])
    except Exception as exc:
        return False, '{0}: {1}'.format(type(exc).__name__, str(exc).strip())
    return True, 'initialized'


################################################################################
# The measurement
################################################################################

def measure(rows=None):
    """Run all three gates over ``rows`` (the eight-row matrix by default).
    Returns a list of dicts."""
    results = []
    for row in (MATRIX if rows is None else rows):
        balanced, balanced_detail = gate_is_balanced(row)
        loaded, load_detail = gate_library_loader(row)
        reactor_ok, reactor_detail = gate_reactor(row)
        results.append({
            'row': row,
            'is_balanced': balanced, 'is_balanced_detail': balanced_detail,
            'library': loaded, 'library_detail': load_detail,
            'reactor': reactor_ok, 'reactor_detail': reactor_detail,
        })
    return results


def format_table(results):
    lines = []
    header = '{0:<4} {1:<16} {2:<40} {3:>11} {4:>9} {5:>9}'.format(
        '', 'form', 'equation', 'is_balanced', 'library', 'reactor')
    lines.append(header)
    lines.append('-' * len(header))
    for r in results:
        lines.append('{0:<4} {1:<16} {2:<40} {3:>11} {4:>9} {5:>9}'.format(
            r['row']['family'], r['row']['form'], r['row']['equation'],
            'ACCEPT' if r['is_balanced'] else 'refuse',
            'ACCEPT' if r['library'] else 'refuse',
            'ACCEPT' if r['reactor'] else 'refuse'))
    return '\n'.join(lines)


################################################################################
# Pinned expectations
################################################################################

#: ``(is_balanced, library, reactor)`` per row, keyed ``family/form``. These are
#: measured values, pinned so a change to any of the three gates has to be a
#: deliberate edit here rather than a silent drift.
EXPECTED = {
    # Every explicit-electron shape whose free-electron COUNT changes fails the
    # per-element comparison on element E, before charge is ever looked at.
    'EI/explicit': (False, False, True),
    'EI/half': (False, False, False),
    'RR/explicit': (False, False, True),
    '3b/e-e-ion': (False, False, True),
    '3b/neutral M': (False, False, True),
    # The scalar forms balance and load.
    'EI/metadata': (True, True, False),
    'RR/metadata': (True, True, False),
    # The matrix's CT row carries no electron at all and fails on CHARGE.
    'CT/ion-atom': (False, False, True),
    # Supplementary.
    'CT/charge-balanced': (True, True, True),
    'EI/half-explicit': (True, True, False),
    '3b/e-e-ion+meta': (False, False, False),
}


@pytest.mark.parametrize('row', MATRIX + SUPPLEMENTARY,
                         ids=lambda r: '{0}/{1}'.format(r['family'], r['form']))
def test_matrix_row(row):
    """Each shape gets the same verdict from each gate that it did when I-108
    measured it."""
    key = '{0}/{1}'.format(row['family'], row['form'])
    expected = EXPECTED[key]
    measured = (gate_is_balanced(row)[0],
                gate_library_loader(row)[0],
                gate_reactor(row)[0])
    assert measured == expected, (
        '{0}: measured (is_balanced, library, reactor) = {1}, expected {2}'.format(
            key, measured, expected))


def test_library_metadata_escape_is_rate_law_borne():
    """The library route stores an electron-count-changing reaction ONLY because
    the rate law carries an ``electrons`` attribute for ``KineticsLibrary.load``
    to copy. The same equation with a plain Arrhenius loses the count and is
    refused by the very ``is_balanced()`` call the loader makes."""
    ei = next(r for r in MATRIX if r['family'] == 'EI' and r['form'] == 'metadata')

    accepted, detail = gate_library_loader(ei)
    assert accepted, detail
    assert 'electrons=+1' in detail

    accepted, detail = gate_library_loader(ei, kinetics=_thermal_bimolecular())
    assert not accepted
    assert 'not balanced' in detail


def _row(key):
    family, form = key.split('/')
    return next(r for r in MATRIX + SUPPLEMENTARY
                if r['family'] == family and r['form'] == form)


def test_explicit_electron_is_a_conserved_element():
    """``is_balanced()`` refuses every explicit-electron form whose free
    electron count changes -- ionisation, radiative recombination and both
    three-body recombination shapes alike -- and refuses each of them on the
    ``E`` per-element comparison, which runs before charge is compared at all."""
    for key in ('EI/explicit', 'EI/half', 'RR/explicit', '3b/e-e-ion', '3b/neutral M'):
        row = _row(key)
        accepted, detail = gate_is_balanced(row)
        assert not accepted, key
        assert 'element mismatch on e' in detail, (key, detail)


def test_the_scalar_cannot_rescue_an_explicit_electron_shape():
    """The signed scalar is folded into the CHARGE comparison only. A shape that
    already fails the per-element comparison on ``E`` is therefore beyond its
    reach: the two representations do not compose, they exclude each other."""
    row = _row('3b/e-e-ion+meta')
    accepted, detail = gate_is_balanced(row)
    assert not accepted
    assert 'element mismatch on e' in detail, detail


def test_matrix_charge_transfer_row_fails_on_charge_not_on_electrons():
    """``Li+ + Ar => Li + Ar`` -- the matrix's own charge-transfer row -- has no
    electron in it anywhere and is refused for not conserving charge. The
    charge-conserving spelling of the same chemistry is accepted by all three
    gates, so the negative result is about that row's stoichiometry and says
    nothing about electron representation."""
    accepted, detail = gate_is_balanced(_row('CT/ion-atom'))
    assert not accepted
    assert 'charge mismatch' in detail, detail
    assert 'element mismatch' not in detail, detail

    assert gate_is_balanced(_row('CT/charge-balanced'))[0]
    assert gate_library_loader(_row('CT/charge-balanced'))[0]
    assert gate_reactor(_row('CT/charge-balanced'))[0]


def test_reactor_gate_is_orthogonal_to_the_balance_gate():
    """The plasma reactor never calls ``is_balanced()``. It accepts explicit
    electron shapes the balance check refuses, and refuses scalar shapes the
    balance check accepts. Neither gate subsumes the other, which is why the
    matrix has to ask all three questions separately."""
    accepted, detail = gate_reactor(_row('RR/explicit'))
    assert accepted, detail
    assert not gate_is_balanced(_row('RR/explicit'))[0]

    assert gate_is_balanced(_row('RR/metadata'))[0]
    accepted, detail = gate_reactor(_row('RR/metadata'))
    assert not accepted
    assert 'ElectronPlacementError' in detail, detail


################################################################################
# The two incompatible definitions of "the E element"
################################################################################

def e_counts_by_both_rules(row):
    """Count element ``E`` on both sides of ``row`` under each of the two rules
    the codebase actually uses, and return them.

    ``Reaction.is_balanced`` counts ATOMS of element ``e``: only the electron
    pseudo-species has any, so a free electron is a conserved element and a
    reaction that creates or destroys one cannot balance.

    ``rmgpy.electron_balance.get_species_electron_count`` -- used by the Chemkin
    and Cantera writers via ``check_electron_balance``, and by
    ``resolve_electron_placement`` step 10 -- counts MINUS THE NET CHARGE for a
    charged species. Under that rule the electron a cation gains is already
    accounted for in the cation's charge, so an ionisation written with explicit
    electrons balances.

    The two rules disagree on exactly the shapes this ticket is about, which is
    the whole substance of the disagreement being re-measured.
    """
    from rmgpy.electron_balance import get_species_electron_count
    spc = _species()
    reactants = [spc[label] for label in row['reactants']]
    products = [spc[label] for label in row['products']]

    atom_rule = (
        sum(1 for s in reactants for a in s.molecule[0].atoms if a.element.symbol == 'e'),
        sum(1 for s in products for a in s.molecule[0].atoms if a.element.symbol == 'e'),
    )
    charge_rule = (
        sum(get_species_electron_count(s) for s in reactants),
        sum(get_species_electron_count(s) for s in products),
    )
    return atom_rule, charge_rule


def test_the_two_E_rules_disagree_on_explicit_ionisation():
    """``Li + e => Li+ + e + e`` is unbalanced under ``is_balanced``'s atom rule
    (1 vs 2) and balanced under the writers' charge rule (1 vs 1).

    So "does the electron balance?" has two different answers in this codebase
    depending on which boundary is asking, and the shape that the balance check
    refuses is the shape the export boundary and the placement resolver's own
    step-10 verification would both accept.
    """
    atom_rule, charge_rule = e_counts_by_both_rules(_row('EI/explicit'))
    assert atom_rule == (1, 2)
    assert charge_rule == (1, 1)

    # Not a one-off. EVERY explicit-electron shape in the matrix is refused by
    # the atom rule and balanced by the charge rule -- the divergence is total,
    # not incidental to ionisation.
    for key in ('EI/explicit', 'EI/half', 'RR/explicit', '3b/e-e-ion', '3b/neutral M'):
        atom_rule, charge_rule = e_counts_by_both_rules(_row(key))
        assert atom_rule[0] != atom_rule[1], (key, atom_rule)
        assert charge_rule[0] == charge_rule[1], (key, charge_rule)
        assert not gate_is_balanced(_row(key))[0], key


################################################################################
# The three reactor-side refusals
################################################################################

#: A family label that is deliberately NOT in ``FAMILY_ELECTRON_PLACEMENT``.
UNDECLARED_FAMILY = 'Li_Electron_Impact_Ionization'
#: The consumption-declaring family the resolver does know about.
DECLARED_CONSUMING_FAMILY = 'Plasma_Electron_Attachment'


def _template_reaction(family, electrons, reactants, products, kinetics):
    from rmgpy.data.kinetics.family import TemplateReaction
    spc = _species()
    return TemplateReaction(
        reactants=[spc[label] for label in reactants],
        products=[spc[label] for label in products],
        family=family, electrons=electrons, reversible=False,
        is_forward=True, kinetics=kinetics)


def refusal_undeclared_family():
    """Refusal 1: an ionisation-shaped Li reaction whose family carries no
    placement declaration at all."""
    rxn = _template_reaction(UNDECLARED_FAMILY, 1, ['Li'], ['Liplus'], _voronov_li(1))
    try:
        from rmgpy.electron_placement import resolve_electron_placement
        resolve_electron_placement(rxn, [_species()['e']])
    except ElectronPlacementError as exc:
        return False, str(exc).strip()
    return True, 'resolved'


def refusal_ionisation_against_consuming_family():
    """Refusal 2: the same ionisation-shaped reaction attributed to the family
    that DOES declare placement -- but declares consumption, net -1."""
    rxn = _template_reaction(DECLARED_CONSUMING_FAMILY, 1, ['Li'], ['Liplus'],
                             _voronov_li(1))
    try:
        from rmgpy.electron_placement import resolve_electron_placement
        resolve_electron_placement(rxn, [_species()['e']])
    except ElectronPlacementError as exc:
        return False, str(exc).strip()
    return True, 'resolved'


def refusal_reactor_metadata_only():
    """Refusal 3: the reactor's own metadata-count guard.

    Reached by calling ``_validate_reactions`` directly, because
    ``initialize_model`` routes every nonzero-``electrons`` reaction through the
    placement resolver first -- so on the whole-pipeline path this guard is the
    backstop behind the resolver, not the first thing a metadata reaction meets.
    """
    spc = _species()
    rxn = Reaction(reactants=[spc['Liplus']], products=[spc['Li']],
                   electrons=-1, reversible=False, kinetics=_badnell_li(-1))
    core_species = [spc['e'], spc['Liplus'], spc['Li']]
    imf = {spc['e']: Y_E0, spc['Liplus']: 0.5, spc['Li']: 0.5 - Y_E0}
    reactor = PlasmaReactor(T_GAS, P0, imf, (T_E, 'K'), n_sims=1, termination=[])
    try:
        reactor._validate_reactions(core_species, [], [rxn], [])
    except Exception as exc:
        return False, '{0}: {1}'.format(type(exc).__name__, str(exc).strip())
    return True, 'validated'


def refusal_double_representation():
    """The half-explicit canonical form -- incident electron as a species, only
    the surplus as the scalar -- put to the resolver with a DECLARED family, so
    that step 1's family check cannot mask what step 2 does with it.

    This is the form ``check_electron_reactant_order`` prescribes and the shipped
    plasma-local-context fixture writes. The resolver classifies it as double
    representation and refuses it outright.
    """
    rxn = _template_reaction(DECLARED_CONSUMING_FAMILY, 1, ['Li', 'e'],
                            ['Liplus', 'e'], _voronov_li(1))
    try:
        from rmgpy.electron_placement import resolve_electron_placement
        resolve_electron_placement(rxn, [_species()['e']])
    except ElectronPlacementError as exc:
        return False, str(exc).strip()
    return True, 'resolved'


REFUSALS = [
    ('undeclared family', refusal_undeclared_family, 'no electron-placement declaration'),
    ('ionisation vs consuming family', refusal_ionisation_against_consuming_family,
     'ionization-shaped'),
    ('reactor metadata-only count', refusal_reactor_metadata_only,
     'cannot distinguish incident-electron order from net electron production'),
]

#: Not one of the three the ticket named -- a fourth, found by measuring.
EXTRA_REFUSALS = [
    ('half-explicit double representation', refusal_double_representation,
     'represents the electron twice'),
]


@pytest.mark.parametrize('name,probe,expected_phrase', EXTRA_REFUSALS,
                         ids=[r[0] for r in EXTRA_REFUSALS])
def test_extra_refusal_fires(name, probe, expected_phrase):
    """The representation the export boundary prescribes is the representation
    the placement boundary refuses. Pinned so the collision is visible."""
    accepted, detail = probe()
    assert not accepted, '{0}: expected a refusal, got {1}'.format(name, detail)
    assert expected_phrase in detail, '{0}: {1}'.format(name, detail)


def test_family_electron_placement_still_holds_exactly_two_consuming_entries():
    """The declaration table still holds exactly two families, both
    single-electron on the reactant side. Nothing ionisation-shaped has been
    added.

    The VALUES were respelled by I-113, which widened the declaration from
    ``(side, count)`` to ``(reactant_count, product_count)`` so that incident
    order could be declared separately from net change; ``('reactants', 1)`` and
    ``(1, 0)`` are the same chemistry in the two spellings. This is the only row
    of this measurement that I-113 deliberately moved. The CHEMISTRY the row
    measures -- two families, one consumed electron each, no ionisation family --
    is unchanged, which is what the row is here to pin.
    """
    from rmgpy.electron_placement import FAMILY_ELECTRON_PLACEMENT
    assert FAMILY_ELECTRON_PLACEMENT == {
        'Plasma_Electron_Attachment': (1, 0),
        'Cation_R_Recombination': (1, 0),
    }


@pytest.mark.parametrize('name,probe,expected_phrase', REFUSALS,
                         ids=[r[0] for r in REFUSALS])
def test_reactor_refusal_still_fires(name, probe, expected_phrase):
    """Each of the three refusals still fires, and still says what it said."""
    accepted, detail = probe()
    assert not accepted, '{0}: expected a refusal, got {1}'.format(name, detail)
    assert expected_phrase in detail, '{0}: {1}'.format(name, detail)


if __name__ == '__main__':
    results = measure()
    supplementary = measure(SUPPLEMENTARY)
    print('THE EIGHT-ROW MATRIX')
    print(format_table(results))
    print()
    print('SUPPLEMENTARY')
    print(format_table(supplementary))
    print()
    for r in results + supplementary:
        key = '{0}/{1}'.format(r['row']['family'], r['row']['form'])
        print('== {0}: {1}'.format(key, r['row']['equation']))
        print('   is_balanced : {0}  ({1})'.format(
            'ACCEPT' if r['is_balanced'] else 'refuse', r['is_balanced_detail']))
        print('   library     : {0}  {1}'.format(
            'ACCEPT' if r['library'] else 'refuse', r['library_detail']))
        print('   reactor     : {0}  {1}'.format(
            'ACCEPT' if r['reactor'] else 'refuse', r['reactor_detail']))
        print()
    print('THE REFUSALS')
    for name, probe, _phrase in REFUSALS + EXTRA_REFUSALS:
        accepted, detail = probe()
        print('== {0}: {1}'.format(name, 'ACCEPTED (!)' if accepted else 'REFUSED'))
        print('   {0}'.format(detail))
        print()
    print('database.directory = {0}'.format(settings['database.directory']))
