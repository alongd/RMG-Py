.. _plasma_reversal:

Resolved reactions and explicit reverses
=======================================

A reaction whose kinetics depends on electron temperature or density and whose
participants carry a resolved electronic or vibrational state must be declared
irreversible. Supply its reverse as a separate irreversible reaction with its own
kinetics. RMG raises ``NonEquilibriumReverseRateError`` with the equation and
resolved participant when an equilibrium reverse would otherwise be inferred.
Gas-temperature-only reactions retain equilibrium reversibility.

Chemkin duplicate import preserves the explicit irreversible declaration when it
combines terms. This also corrects the old unresolved-chemistry behavior: two
same-direction irreversible duplicate terms now reload as one irreversible
multi-rate reaction, rather than defaulting to a reversible reaction and
fabricating ``kb`` from ``kf``. Chemkin and Cantera writers check the original rate before
converting or expanding it; Chemkin readers check after attaching TDEP and again
before returning imported reactions. Arkane thermal fitting, direct Arkane
Chemkin output, pressure-dependence input saving, RMS output, and legacy library
output may refuse resolved nonthermal rates with the same error before conversion
or serialization. RMS writes an explicit
``reversible: false`` for irreversible thermal reactions.

The source inventory
--------------------

``scripts/generate_reversal_census.py`` tokenizes Python and Cython source in
``rmgpy/`` (including ``rmgpy/tools/``) and ``arkane/``. It discovers calls to
reaction direction and reverse-rate primitives, reaction builders, reader and
exporter entry points, and executable expressions in Jinja templates. It does
not import a built extension or scan comments and documentation as executable
calls. Structural Species/Molecule/Group calls appear in the inventory too; their
verdicts explain why they cannot decide reaction direction or produce kinetics.

``reversal_census_verdicts.json`` records the reviewed verdict, reason and test
selectors for each stable site identity. Discovery never assigns verdicts.
``reversal_census.json`` is the generated inventory joined to those verdicts.
Regenerate it from the repository root with::

    python scripts/generate_reversal_census.py \
        --verdicts documentation/source/users/rmg/reversal_census_verdicts.json \
        --output documentation/source/users/rmg/reversal_census.json

Run ``test/rmgpy/data/kinetics/superelasticCensusTest.py`` to check that discovery,
verdicts, generated output and test selectors agree. A new unclassified site, a
removed verdict, or a disappeared site fails the gate. Line numbers are regenerated
evidence; the site identity uses its file, enclosing function, primitive and
occurrence number. Tests cover shared policy primitives and reachable branch
witnesses; the verdicts identify delegation and structural exclusions explicitly.
