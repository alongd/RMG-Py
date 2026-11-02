EEDF tables
===========

``EEDFTable.load(..., provider_kind='loki-table')`` and
``EEDFProvider(..., provider_kind='loki-table')`` default to production LoKI
identity. Pass ``provider_kind='maxwellian-table'`` for a validated comparator.
The same provider class exposes rows, reaction rates, transport and mean
energy for both kinds. The kind is immutable after provider construction and
is checked against the artifact before rows become available.

Comparator generation is an offline operation. Runtime table/provider imports
do not import the generator. Comparator numerical validation is distinct from
production qualification: ``accepted`` remains false and production admission
refuses. Source field coordinates are returned as ``reference_EN_Td``;
``EN_Td`` and physical field diagnostics are unavailable for the comparator.

.. automodule:: rmgpy.solver.eedf
    :members:

.. automodule:: rmgpy.solver.eedf_provider
    :members: EEDFProvider

.. automodule:: rmgpy.tools.eedf.generate
    :members: generate, model_inputs

.. automodule:: rmgpy.tools.eedf.maxwellian
    :members: generate_maxwellian
