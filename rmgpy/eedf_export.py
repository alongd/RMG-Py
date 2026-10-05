"""Qualification boundary for mechanism exports from EEDF reactor jobs."""

from rmgpy.exceptions import EEDFExportError


def require_non_eedf_reactor_mode(rmg, target):
    """Refuse a writer job whenever any attached reactor is in EEDF mode.

    Low-level serializers can recognize marker kinetics, but they cannot infer
    reactor mode from an ordinary reaction model.  Writer listeners retain that
    context and therefore enforce the mode-level qualification boundary before
    opening an output file.
    """
    if any(getattr(system, 'eedf_mode', False)
           for system in getattr(rmg, 'reaction_systems', ())):
        raise EEDFExportError(
            'EEDF export to {0} requires qualification into a supported '
            'non-EEDF rate, even when the current mechanism has no EEDFChannel '
            'markers.'.format(target))
