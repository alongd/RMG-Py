Operating branches
==================

The branch scanner records the steady states of a plasma reactor at one absorbed power, continues
the stable ones over power under a stated certificate, and matches a run's terminal state to one
recorded branch.  See :ref:`plasmaoperatingbranch` for the input option that uses it.

.. automodule:: rmgpy.tools.eedf.branches
    :members: scan, continue_branch, record_continuation, certificate_issue, tolerance_issue,
        reference_issue, matchable_issue, matching_branch, declared_branch, path_consistency,
        certified_coverage, aliasing_radius, PlasmaReactorAdapter, ReactorRun, Seed
