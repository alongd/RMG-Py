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
independent state lists.

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
updates, so other products stay unresolved. A later ``SET_STATE`` targeting
the same product replaces both values. Use ``''`` and ``-1`` for unresolved
fields; ``['SET_STATE', '*1', '', -1]`` clears both fields. Product Group
templates receive the same values as singleton lists (or empty lists for
unresolved values).

Reactant state is preserved on the original molecules. The temporary merged
recipe graph and unassigned products are unresolved. Automatic reversal of
``SET_STATE`` is refused because the action does not specify the previous
state; declare such a family ``reversible = False`` and provide any reverse
channel separately. Direct ``ReactionRecipe.apply_forward`` can assign state
to a connected structure; use ``KineticsFamily.apply_recipe`` for multiple
products.
