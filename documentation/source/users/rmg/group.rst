.. _group:

********************
Group Representation
********************

Group representations are used to represent molecular substructures within RMG.
These are commonly used for identifying functional groups for use in both the
thermo and kinetic databases.

For syntax of how to define groups, see :ref:`rmgpy.molecule.adjlist`.

State constraints
=================

Groups accept optional ``electronicstate`` and ``vibrationallevel`` headers.
Each is a list-valued constraint, independent of the other. A missing header
or an empty list admits only the corresponding unresolved molecule value
(``''`` for electronic state, ``-1`` for vibrational level). ``x`` admits every
value, including unresolved. Explicit lists admit exactly their listed values::

    electronicstate [A3Su+,B3Pg]
    vibrationallevel [0,1]
    1 *1 N u0 p1 c0 {2,T}
    2 *2 N u0 p1 c0 {1,T}

Use ``electronicstate x`` and ``vibrationallevel x`` for a template that admits
both resolved and unresolved molecules. Quoted Python lists can represent
an unresolved electronic value or a token containing a comma, for example
``electronicstate ['', 'a,b.c']``. An explicit vibrational list may include
``-1`` to admit its unresolved value. The molecular token grammar and integer
bounds also apply to group values. State constraints participate in Group
isomorphism and subgraph specificity; list order does not affect matching.
They survive adjacency-list I/O, pickle, copying, and deepcopy. Copies have
independent state lists. Logical tree nodes combine their components' declared
state domains; structural negation does not invert those domains. An empty
logical node admits unresolved species only.

A reaction family also needs ``allowExcitedReactants = True`` in ``groups.py``
to admit a resolved reactant. Its Python attribute is
``KineticsFamily.allow_excited_reactants`` and defaults to ``False``. A true
flag still requires matching group constraints; a wildcard template alone
cannot opt a family in. The flag applies in both generation directions.
Existing families require no changes to keep their unresolved behavior.

Product state assignment
========================

A family recipe can set both state fields on a product::

    allowExcitedReactants = True
    reversible = False
    recipe(actions=[
        ['SET_STATE', '*1', 'A3Su+', 1],
    ])

``SET_STATE`` targets the final product containing the uniquely labeled atom.
The assignment takes effect after graph actions, product splitting, and
updates, and applies only to that product. A later ``SET_STATE`` targeting
the same product replaces both values. Use ``''`` and ``-1`` for unresolved
fields; ``['SET_STATE', '*1', '', -1]`` clears both fields. Product Group
templates receive the same values as singleton lists (or empty lists for
unresolved values).

Reactant state is preserved on the original molecules. Products retaining exactly
one input species' atoms inherit its state, in either recipe direction. The
merged working graph carries no species state. Other products are unresolved
unless ``SET_STATE`` assigns them. Product state constraints are checked
independently of ``allowExcitedReactants`` and multiplicity constraints. Each
state-constrained product template requires a distinct matching product.

Automatic reversal of ``SET_STATE`` uses the input template's unique declared
state (an absent header means unresolved). When species merge or split, the
inverse restores each affected original species, including species without a
forward ``SET_STATE`` action. Wildcards or multiple possible original states
raise ``StateReversalError`` naming the family during loading; declare such a
family ``reversible = False`` and provide its reverse channel separately.
Reversible recipes that change species boundaries and cannot carry declared
state constraints are likewise refused. Direct ``ReactionRecipe.apply_forward``
can assign state to a connected structure; use ``KineticsFamily.apply_recipe``
for multiple products. ``x`` is reserved for Group wildcards and cannot be a
Molecule electronic-state token or a ``SET_STATE`` value. Groups containing
only state headers round-trip without atom lines. Constructor constraints must
be lists, with valid state tokens or integer vibrational levels; booleans are
not levels.


Data provenance and training
============================

Matching a wildcard parent does not authorize its children's data. Resolved
estimates check alias targets, ancestors, and the sources of cached averages.
Unmatched sources raise ``StateProvenanceError``; averaging provenance is kept
independently of comments and verbosity. Rate rules must cover the matched
template's structure and state constraints. When a bare template has no
reactant context, this check conservatively treats state-bearing groups as
potential resolved inputs. Templates selected from unresolved reactions keep
the existing averaging behavior.

Resolved ring averaging, polycyclic decomposition, radical saturation, halogen
replacement, surface desorption, adsorption group-tree estimation, and the transport
Lennard-Jones fallback are refused with ``StateProvenanceError``. These
transformations have no supported whole-species state provenance model.
Direct group values and aliases remain available when every source node matches.
All SIDT adsorption estimates are refused for resolved input, including fitted
roots with no group, because their values lack matched state provenance.
QM and ML evaluators refuse resolved input before any state-blind representation
can supply an estimate. The topology-only connectivity probe is not proof of
matched fitted data.

Training ingestion and automatic tree fitting refuse resolved reactants or
products with ``ResolvedStateTrainingError``, before modifying data or structures.
Read-only reaction matching continues to support opted-in state constraints.
Tree root merging retains their conjunctive state domain; incompatible domains
raise ``StateConstraintMergeError`` rather than silently clearing constraints.

Pickle retains derived-entry provenance. Modern and legacy text writers refuse
state-bearing group entries with derived data when the text format would erase
those sources. The writers inspect logical members as well as literal Group
headers; unresolved-only text output is unchanged.

Resolved frequency and electron data
------------------------------------

Tree selection for a resolved molecule checks the complete atom pattern as well
as both state fields, including an ``Others-*`` fallback. An unmatched fallback
raises ``StateProvenanceError`` before it can supply a rate template or data.
Unresolved tree selection retains the existing fallback convention.

Characteristic-frequency estimation refuses resolved molecules with
``StateProvenanceError``. This includes the fixed ring C-H frequencies, which
have no matched group provenance. Statmech libraries and depositories can still
supply a record with an exactly matching molecular state.

The electron thermo-library convention still ignores radical multiplicity, but
both declared state fields must match the library record. Records with a different
state are skipped, so the ordering of ground and resolved records cannot prevent
a matching record from supplying thermo. A resolved electron raises
``StateProvenanceError`` only when no matching record is found. Library searches
continue past a library that lacks the declared electron state.

Kinetics reconstruction from comment/source records also refuses resolved
reactions with ``StateProvenanceError`` because those supplied sources have no
original-molecule matching proof. Existing unresolved reconstruction is unchanged.

The static data-route census is checked by
``test/rmgpy/data/stateDataCensusTest.py``. It includes every function under
``rmgpy/`` and ``arkane/``, including Cython functions and property accessors,
as a conservative superset of data suppliers and their helpers. This covers QM,
ML, species/model assignment, energy transfer, thermo, statmech and ESS file
loaders as well as database routes, with explicit verdicts
in ``test/rmgpy/data/stateDataCensus.json``. Added or changed functions require
a reviewed classification; ``python scripts/state_data_census.py --report <path>``
generates a source-located table from the AST and the checked manifest.

QM and ML thermo estimation refuse resolved input with ``StateProvenanceError``
before initialization, cache reads, requests, or prediction. The refusal applies
to independently callable estimation boundaries as well as database wrappers.
Default energy-transfer estimation and Arkane statmech file loading also refuse
resolved inputs. Explicit species models remain explicit inputs; no fitted
state model or new database record is supplied by these refusals.

Census ``carried`` verdicts for numeric operators mean they evaluate or convert
explicitly supplied models without selecting an independent ground template.
Non-acquisition utilities are included in the inventory but supply no scientific
value for a species. Their classifications do not prove that generic graph or
format helpers preserve state. State-erasing representations can never justify
matched provenance at an acquisition boundary. New or changed sites require a
reviewed reason and evidence; the checker assigns no default verdict.

Heat-capacity limits and symmetry arithmetic remain reachable for resolved
molecular identities. These helpers select no database record; their thermo
callers provide state-matched data or an explicit numeric model.
