.. _kineticsDatabase:

*****************
Kinetics Database
*****************
This section describes the general usage of RMG's kinetic database. See :ref:`kinetic-database-modification` for 
instructions on modifying the database.

Pressure independent reaction rates in RMG are calculated using a modified 
Arrhenius equation, designating the reaction coefficient as :math:`k(T)` at 
temperature :math:`T`.

.. math:: k(T) = A\left(\frac{T}{T_0}\right)^ne^{-(E_0 + \alpha \Delta H_{rxn})/(RT)}

:math:`R` is the universal gas constant. The **kinetic parameters** determining 
the rate coefficient are:

* :math:`A`:	the pre-exponential A-factor 

* :math:`T_0`:	the reference temperature

* :math:`n`:	the temperature exponent

* :math:`E_0`:	the activation energy for a thermoneutral reaction (a barrier height intrinsic to the kinetics family)

* :math:`\alpha`:	the Evans-Polanyi coefficient (characterizes the position of the transition state along the reaction coordinate, :math:`0 \le \alpha \le 1`)

* :math:`\Delta H_{rxn}`: the enthalpy of reaction

When Evans-Polanyi corrections are used, :math:`\Delta H_{rxn}` is calculated
using RMG's thermo database, instead of being specified in the kinetic database.
When Evans-Polanyi corrections are not used, :math:`\Delta H_{rxn}` and :math:`\alpha`
are set to zero, and :math:`E_0` is the activation energy of the reaction.

Libraries
=========
Kinetic libraries delineate kinetic parameters for specific reactions. 
RMG always chooses to use kinetics from libraries over families. If multiple libraries
contain the same reaction, then precedence is given to whichever library is
listed first in the input.py file.

For combustion mechanisms, you should always use at least one small-molecule 
combustion library, such as the pre-packaged *BurkeH2O2* and/or *FFCM1*
for natural gas.
The reactions contained in these libraries are poorly estimated by kinetic 
families and are universally important to combustion systems.

Kinetic libraries should also be used in the cases where:

* A set of reaction rates were optimized together
* You know the reaction rate is not generalizable to similar species (perhaps due to catalysis or aromatic structures)
* No family exists for the class of reaction
* You are not confident about the accuracy of kinetic parameters

`A full list of the libraries can be found on the RMG website <https://rmg.mit.edu/database/kinetics/libraries/>`_ (please allow 1-2 minutes for the website to load).


.. _kinetics_families_db:

Families
========
Allowable reactions in RMG are divided up into classes called **reaction families**.
All reactions not listed in a kinetic library have their kinetic parameters 
estimated from the reaction families. 

Each reaction family contains the files:

* groups.py containing the recipe, group definitions, and hierarchical trees
* training.py containing a training set for the family
* rules.py containing kinetic parameters for rules

A family may optionally also contain:

* quarantine.py declaring that some of the family's data must not enter a
  quantitative mechanism (see :ref:`kinetics_quarantine`)

`A full list of the kinetic families can be found on the RMG website <https://rmg.mit.edu/database/kinetics/families/>`_ (please allow 1-2 minutes for the website to load).

Recipe
------
The recipe can be found near the top of groups.py and describes the changes in
bond order and radicals that occur during the reaction. Reacting atoms are
labelled with a starred number. Shown below is the recipe for the H-abstraction 
family.

.. image:: images/Recipe.png
	:scale: 65%
	:align: center

The table below shows the possible actions for recipes. The arguments are given 
in the curly braces as shown above. For the order of bond change in the 
Change_Bond action, a -1 could represent a triple bond changing to a double 
bond while a +1 could represent a single bond changing to a double bond. 

+------------+-----------------+---------------------+------------------+
|Action      |Argument1        |Argument2            |Argument3         |
+============+=================+=====================+==================+
|Break_Bond  |First bonded atom|Type of bond         |Second bonded atom|
+------------+-----------------+---------------------+------------------+
|Form_Bond   |First bonded atom|Type of bond         |Second bonded atom|
+------------+-----------------+---------------------+------------------+
|Change_Bond |First bonded atom|Order of bond change |Second bonded atom|
+------------+-----------------+---------------------+------------------+
|Gain_Radical|Specified atom   |Number of radicals   |                  |
+------------+-----------------+---------------------+------------------+
|Lose_Radical|Specified atom   |Number of radicals   |                  |
+------------+-----------------+---------------------+------------------+

Change_Bond order cannot be directly used on benzene bonds. During generation,
aromatic species are kekulized to alternating double and single bonds such that
reaction families can be applied. However, RMG cannot properly handle benzene bonds 
written in the kinetic group definitions.

Training Set vs Rules
---------------------
The training set and rules both contain trusted kinetics that are used to fill in
templates in a family. The **training set** contains kinetics for specific reactions,
which are then matched to a template. The kinetic **rules** contain kinetic 
parameters that do not necessarily correspond to a specific reaction, but have 
been generalized for a template.

When determining the kinetics for a reaction, a match for the template
is searched for in the kinetic database. The three cases in order
of decreasing reliability are:

#. Reaction match from training set
#. Node template exact match using either training set or rules
#. Node template estimate averaged from children nodes

Both training sets and reaction libraries use the observed rate, but rules must
first be divided by the degeneracy of the reaction. For example, the reaction
CH4 + OH --> H2O + CH3 has a reaction degeneracy of 4. If one performed an
experiment or obtained this reaction rate using Arkane (applying the correct
symmetry), the resultant rate parameters would be entered into libraries and
training sets unmodified. However a kinetic rule created for this reaction must
have its A-factor divided by 4 before being entered into the database. 

The reaction match from training set is accurate within the documented uncertainty for that
reaction. A template exact match is usually accurate within about one order
of magnitude. When there is no kinetics available for for the template in
either the training set or rules, the kinetics are averaged from the children
nodes as an estimate. In these cases, the kinetic parameters are much less reliable.
For more information on the estimation algorithm see :ref:`kinetics`. 

The training set can be modified in training.py and the rules can be modified in
rules.py. For more information on modification see :ref:`kinetic-training-set` and :ref:`kinetic-rules`.

.. _kinetics_quarantine:

Quarantined data
----------------
Some data in the database is real and worth keeping, yet its *evaluation* has lost
the semantics of its source model -- the number the rate law returns is not the
number the data means. Rates fitted under conditions that a mechanism does not
reproduce are the usual case: an electrode electron-transfer rate evaluated in a
gas-phase mechanism, where there is no electrode and no reference potential, is
computed correctly and means nothing.

Deleting such data destroys provenance; leaving it in circulation puts a
meaningless rate into a mechanism that then reports success. **Quarantine** is the
third option: keep the data, and refuse loudly at the boundary where it would
become a quantitative claim.

A family declares a quarantine by carrying a ``quarantine.py`` sidecar next to its
``groups.py`` and ``rules.py``::

    name = "Some_Family/quarantine"
    state = "QUARANTINED FOR QUANTITATIVE PLASMA USE"
    appliesToKineticsClass = "Marcus"
    reason = "electrochemical reference/domain unavailable"

``state`` and ``reason`` are free text and are quoted back in the error.
``appliesToKineticsClass`` names a class in :mod:`rmgpy.kinetics` and is the
criterion: **the manifest never names an entry.** RMG computes which entries are
affected by testing every rate rule and every training entry of the family against
that class, so the quarantine cannot fall out of step with the data. Adding an
entry of that class quarantines it automatically; refitting one to another kinetics
class releases it automatically; renumbering or relabelling entries changes
nothing. A misspelled class name raises ``DatabaseError`` when the family loads,
rather than producing a manifest that announces a quarantine while gating nothing.

When a reaction whose kinetics come from quarantined data is admitted to an RMG
reaction model, the run stops with ``QuarantinedKineticsError``, naming the family,
the rule or training-entry provenance, the generated reaction, the kinetics class,
the reason, and the path of the manifest. The reaction is **not** dropped,
averaged, zeroed, reversed, made irreversible, or given a substitute rate: each of
those would let the run report success on a mechanism that is quietly wrong.

The gate sits inside RMG's mechanism builder and nowhere else. The database still
loads, the family stays registered and still generates reactions, and the rate law
itself still evaluates, so the data remains usable for refitting, for provenance
audits, and for any consumer that supplies the missing reference. Lifting a
quarantine means deleting the ``quarantine.py`` file.

Families with no ``quarantine.py`` -- which is all of them by default -- are
completely unaffected.
