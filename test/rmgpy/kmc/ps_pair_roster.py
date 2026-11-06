"""The PS reference pair roster (moved here from ``rmgpy.kmc.windows``; nothing
in the library called it, only the reference windows tests and acceptance
scripts do).  It builds the polystyrene reference chains by name."""

from itertools import combinations_with_replacement

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import windows as win


def ps_pair_roster(certifier, family, probe_dps=(2, 3)) -> tuple:
    """The participant-kind combinations a bimolecular direction needs, from
    the family's own reactant templates: RMG applies them (``member_cases``,
    both slot orders) to the representative PS chains of every kind pair at
    the probe DPs; a pair is in the roster of that direction iff some member
    has every species in the declared domain.  ((direction, (kind, kind)),
    ...), kinds sorted."""
    constants = win.family_constants(family)
    reference = b.PROVIDERS["polystyrene"]
    out = []
    for direction, molecularity in constants.directions:
        if molecularity != 2:
            continue
        for kinds in combinations_with_replacement(win.PARTICIPANT_KINDS, 2):
            found = False
            for dp in probe_dps:
                for first in win.kind_specs(kinds[0], dp)[:2]:
                    for second in win.kind_specs(kinds[1], dp)[:2]:
                        builds = [reference.build("k0", first), reference.build("k1", second)]
                        if not all(win.dom.check(x.molecule).accepted for x in builds):
                            continue
                        cases = win.member_cases(certifier, family, builds, direction)
                        if any(not win.species_domain_refusal(c.witness) for c in cases):
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                out.append((direction, tuple(kinds)))
    return tuple(out)
