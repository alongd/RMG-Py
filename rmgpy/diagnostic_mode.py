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

"""Canonical run-mode qualifications shared by artifact writers.

Mode is derived from each reactor, never stored as process-global run state.
Non-plasma callers retain their existing unlabelled output.
"""

#: The name a prescribed-electron-temperature run must call itself by, verbatim.
PRESCRIBED_TE_MODE = 'LEGACY: prescribed-Te transport diagnostic'

#: The disclaimer that must appear verbatim in every artifact such a run writes.
#: The subscript is U+2091 (LATIN SUBSCRIPT SMALL LETTER E), so the sentence reads
#: "... prescribed Te and nₑ; ...".
PRESCRIBED_TE_DISCLAIMER = (
    'Species evolution conditional on prescribed Te and nₑ; '
    'absolute electron density and ionisation degree are not predicted.'
)

#: The same sentence for artifact formats that cannot carry non-ASCII.  Chemkin
#: input decks are the case that forces this: they are consumed by Fortran-era
#: parsers with no encoding declaration.  Only the subscript changes; the rest of
#: the sentence is identical, byte for byte.
PRESCRIBED_TE_DISCLAIMER_ASCII = (
    'Species evolution conditional on prescribed Te and n_e; '
    'absolute electron density and ionisation degree are not predicted.'
)

SOLVED_TE_MODE = 'LEGACY: solved-Te energy-balance engineering intermediate'
SOLVED_TE_DISCLAIMER = (
    'Electron temperature and species evolution are solved conditional on specified '
    'absorbed power and the declared Maxwellian, wall and sheath closures; coupling '
    'of an external power source to the plasma is not predicted.'
)
EEDF_MODE = 'EEDF energy-balance engineering intermediate'
EEDF_DISCLAIMER = (
    'Electron kinetics use the declared EEDF table; reduced field and species '
    'evolution are solved conditional on specified absorbed power, empirical laws '
    'and wall/sheath closures. This label does not establish terminal qualification '
    'or predict external power-source coupling.'
)
MAXWELLIAN_MODE = 'MODEL-FORM COMPARATOR: Maxwellian-table energy-balance diagnostic'
MAXWELLIAN_DISCLAIMER = (
    'Electron kinetics use a prescribed Maxwellian on the source-table grid. '
    'This is a model-form comparator; it is not a production EEDF result and '
    'does not establish EEDF self-consistency or a physical reduced electric field.'
)
_MODE_DISCLAIMERS = {
    PRESCRIBED_TE_MODE: PRESCRIBED_TE_DISCLAIMER,
    SOLVED_TE_MODE: SOLVED_TE_DISCLAIMER,
    EEDF_MODE: EEDF_DISCLAIMER,
    MAXWELLIAN_MODE: MAXWELLIAN_DISCLAIMER,
}


#: Prefix used for the mode-declaration line, so the mode name can be found in an
#: artifact by an exact-string search without the reader having to know the
#: artifact's comment syntax.
MODE_LINE_PREFIX = 'RMG diagnostic mode: '


def diagnostic_mode_of(reaction_system):
    """
    Return the diagnostic-mode name a reaction system declares for itself, or
    ``None`` if it declares none.

    Duck-typed on purpose.  Reactors that are unconditional predictions -- simple,
    liquid, surface, mb-sampled -- declare nothing and get ``None``, which is what
    keeps their artifacts byte-identical.
    """
    mode = getattr(reaction_system, 'diagnostic_mode', None)
    if mode is not None and (not isinstance(mode, str) or mode not in _MODE_DISCLAIMERS):
        raise ValueError('Unknown diagnostic mode: {0!r}'.format(mode))
    return mode


def _reaction_systems_of(obj):
    """
    Normalise ``obj`` to a list of reaction systems.

    Accepts an :class:`~rmgpy.rmg.main.RMG` object (anything carrying a
    ``reaction_systems`` attribute), a single reaction system, a list/tuple of
    them, or ``None``.
    """
    if obj is None:
        return []
    systems = getattr(obj, 'reaction_systems', None)
    if systems is not None:
        return list(systems)
    if isinstance(obj, (list, tuple, set)):
        return list(obj)
    return [obj]


def run_diagnostic_modes(obj):
    """
    Return the ordered, de-duplicated list of diagnostic-mode names declared by
    the reaction systems of ``obj``.  Empty for a run with no such reactor.
    """
    modes = []
    for system in _reaction_systems_of(obj):
        mode = diagnostic_mode_of(system)
        if mode and mode not in modes:
            modes.append(mode)
    return modes


def run_provenance_lines(obj, ascii_only=False):
    """
    Return the provenance lines that every artifact of this run must carry, as a
    list of bare strings with no comment syntax applied -- the caller wraps them
    in whatever its format uses.

    Returns ``[]`` when no reaction system declares a diagnostic mode, which is
    the case for every non-plasma run.  Callers therefore need no conditional of
    their own, and a run that is not in a diagnostic mode writes exactly the bytes
    it wrote before.

    Set ``ascii_only=True`` for formats that cannot carry the U+2091 subscript.
    """
    modes = run_diagnostic_modes(obj)
    if not modes:
        return []
    lines = []
    for mode in modes:
        disclaimer = _MODE_DISCLAIMERS[mode]
        if ascii_only:
            disclaimer = disclaimer.replace('nₑ', 'n_e')
        lines.extend([MODE_LINE_PREFIX + mode, disclaimer])
    return lines


def comment_block(obj, prefix, ascii_only=False, trailing_blank=True):
    """
    Return the provenance lines of ``obj`` rendered as a comment block using
    ``prefix`` (e.g. ``'! '`` for Chemkin, ``'# '`` for CSV/Python/YAML), as a
    single string ending in a newline -- or ``''`` when there is nothing to say.

    ``''`` is the important return value: concatenating it in front of a file's
    body is a no-op, so a non-diagnostic run's artifact is unchanged.
    """
    lines = run_provenance_lines(obj, ascii_only=ascii_only)
    if not lines:
        return ''
    block = ''.join(prefix.rstrip() + ' ' + line + '\n' if line else prefix.rstrip() + '\n'
                    for line in lines)
    if trailing_blank:
        block += prefix.rstrip() + '\n'
    return block
