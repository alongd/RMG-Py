Plasma run provenance
====================

Modes are derived from each reactor, with literal labels shared by artifact
writers. A mixed run reports each distinct mode; no process-global mode is
inherited by a subsequent run.

* P: ``LEGACY: prescribed-Te transport diagnostic``.
* S: ``LEGACY: solved-Te energy-balance engineering intermediate``.
* E: ``EEDF energy-balance engineering intermediate``.
* M: ``MODEL-FORM COMPARATOR: Maxwellian-table energy-balance diagnostic``.

P and S use the existing electron-temperature-law kinetics when
``electronKinetics`` is absent. E uses ``loki-table``. M uses
``maxwellian-table`` and carries this qualification verbatim:

   Electron kinetics use a prescribed Maxwellian on the source-table grid.
   This is a model-form comparator; it is not a production EEDF result and
   does not establish EEDF self-consistency or a physical reduced electric field.

Mode labels describe the model choice, not acceptance or terminal
qualification. Requesting provenance for an EEDF reactor without its required
energy balance refuses. Constructor error identities can describe incomplete
configuration without publishing valid run provenance.

.. automodule:: rmgpy.diagnostic_mode
   :members:

.. autoexception:: rmgpy.exceptions.PrescribedTeDiagnosticError
