"""
Species identity by RMG's own isomorphism, shared by the polymer-mom
serializer (``rmgpy.polymer``: the rebuild signature) and the kMC producer
(``rmgpy.kmc.mom_producer``: the recorded gas mask's structures).  This
module is byte-identical in both checkouts.
"""

from rmgpy.molecule import resonance


# These identity verdicts deliberately are not truthy/falsey stand-ins.  A
# caller must refuse rather than silently treating an unavailable comparison
# as a different structure.
UNDECIDABLE = object()
MALFORMED = object()


def _forms(molecule):
    """Copies of ``molecule`` and of all its resonance forms, unfiltered
    (aromatic and Kekule forms both kept), each with the original's
    multiplicity: RMG's resonance generation resets it to the radical count
    plus one, which turns a singlet diradical's forms into triplets."""
    original = molecule.copy(deep=True)
    multiplicity = original.multiplicity
    forms = [original]
    candidate = molecule.copy(deep=True)
    if candidate.get_net_charge() != 0:
        raise _Undecidable
    candidate.reactive = True
    for form in resonance.generate_resonance_structures(candidate, keep_isomorphic=False,
                                                        filter_structures=False):
        form.multiplicity = multiplicity
        forms.append(form)
    return forms


class _Undecidable(Exception):
    """Resonance expansion is unavailable, rather than a mismatch."""


def same_structure(first, second):
    """Whether two lists of molecules (a species' resonance forms, or one
    drawing) are one species by RMG's own test, the one it merges species
    by: the same multiplicity (of the given molecules, never of generated
    forms), and some resonance form of the one isomorphic (strict:
    electrons, charges, multiplicity) to some resonance form of the other,
    the forms of both sides generated unfiltered on copies.  Symmetric.
    Tautomers are distinct graphs, so never one species.  An empty list or
    generation unavailable: ``UNDECIDABLE``.  A mixed-spin list is
    ``MALFORMED``.  The inputs are left untouched."""
    try:
        if not first or not second:
            return False
        first = [m.copy(deep=True) for m in first]
        second = [m.copy(deep=True) for m in second]
        if len({m.multiplicity for m in first}) != 1 or len({m.multiplicity for m in second}) != 1:
            return MALFORMED
        if first[0].multiplicity != second[0].multiplicity:
            return False
        if any(a.is_isomorphic(b, strict=True) for a in first for b in second):
            return True
        forms_first = [f for m in first for f in _forms(m)]
        forms_second = [f for m in second for f in _forms(m)]
        return any(a.multiplicity == b.multiplicity and a.is_isomorphic(b, strict=True)
                   for a in forms_first for b in forms_second)
    except _Undecidable:
        return UNDECIDABLE
    except Exception:
        return UNDECIDABLE
