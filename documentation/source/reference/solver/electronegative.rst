Electronegative wall qualification
=================================

.. automodule:: rmgpy.solver.electronegative
    :members:

``PlasmaReactor`` exposes accepted-state qualification through
``monitor_electronegative_wall`` and ``electronegative_wall_manifest``,
transport decomposition through ``compute_ion_wall_components`` and
``compute_anion_transport_data``, accepted reaction classifications through
``compute_reference_reaction_data``, and the
analytic charged-wall derivative through
``compute_electronegative_wall_jacobian``. See the input-file documentation
for required external references and frozen thresholds.
