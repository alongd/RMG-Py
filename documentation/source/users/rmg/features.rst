.. _features:

********************
Overview of Features
********************

**Thermodynamics estimation using group additivity.**
	Group additivity based on Benson's groups provide fast and reliable thermochemistry estimates. A standalone utility for estimating heat of formation, entropy, and heat capacity is also included.

**Rate-based model enlargement.**
 	Reactions are added to the model based on their rate, fastest first.

**Rate-based termination.** 
	The model enlargement stops when all excluded reactions are slower than a given threshold.
	This provides a controllable error bound on the kinetic model that is generated.

**Extensible libraries.**
	Ability to include reaction models on top of the provided reaction families.

**Independent electron channels from libraries.**
    Opposite irreversible reactions from different libraries remain separate when
    an explicit electron participates and the two directions differ. Each retains
    its own rate law; other cross-library duplicates keep first-library priority.
    Every reaction with an electron participant is refused from pressure-dependent
    networks, whether the electron is explicit, a specific collider, a kinetics efficiency or
    coverage key, stored as signed reaction metadata, or required by an owner
    placement declaration. Nested and cached kinetics and every declared species
    physical record are checked too, including conformers, modes, thermo and
    transport. Object-dtype arrays anywhere in these records refuse. An
    unsupported exact type is refused. Only known RMG records and standard
    containers are inspected; user subclasses and
    foreign carriers cannot establish electron absence. Restored networks are
    checked at routing and computation entries after circular state is installed.
    This includes reversible and single-library entries and direct module-level
    network computations. Network admission validates the incoming channel;
    computation entries freshly validate current channels and inspect each
    shared physical record once within that call. No verdict survives a call.
    Electron rates remain explicit; electron-impact kinetics are not falloff chemistry.
    Disable the elementary high-pressure and pressure-dependent routing flags and
    cached network kinetics for these entries, or disable pressure dependence.
    Admission checks the entire batch before registration, so a refused seed or
    restart leaves the model unchanged and can be retried with corrected routing.
    If reconstruction or preflight raises, original entry comments are restored
    and the original exception propagates.

**Pressure-dependent reaction networks.**
	Dissociation, combination, and isomerization reactions have the potential to have rate coefficients that are dependent on both temperature and pressure, and RMG is able to estimate both for networks of arbitrary complexity with a bounded error.
	
**Simultaneous mechanism generation for several conditions.**
	Concurrent generation of a reaction mechanism over multiple temperature and pressure conditions. 
	Mechanisms generated this way are valid over a range of reaction conditions.

**Dynamic simulation to a target conversion or time.**
	Often the desired simulation time is not known *a priori*, so a target conversion is preferred.

**Transport properties estimation using group additivity.**
	The Lennard-Jones sigma and epsilon parameters are estimated using empirical correlations (based on a species' critical properties and acentric factor).
	The critical properties are estimated using a group-additivity approach; the acentric factor is also estimated using empirical correlations.
	A standalone application for estimating these parameters is provided, and the output is stored in CHEMKIN-readable format.

Resolved-state export boundaries
-------------------------------

Electronic-state and vibrational-level identities are retained in paired Chemkin
files, RMS YAML, Cantera YAML notes, and modern kinetics dictionaries. Paired
kinetics and dictionary files carry matching generation headers. Load the matching
dictionary with a stamped Chemkin file; mixed generations and missing dictionaries
raise a named identity error. Native Cantera conversion supports resolved species
with ``use_chemkin_identifier=True`` and a complete species inventory; its plain
label mode refuses resolved states. Flux diagrams and simulation profile CSVs
validate names against full molecular identities.

Label-only HTML reports, saved RMG input decks, model-merge output, observable
comparison plots, uncertainty reports, and native Julia RMS conversion refuse
resolved states with a named error. These restrictions affect output identity;
they do not add thermodynamic or transport rules.


Electronic and vibrational state identity is retained in supported mechanism
exports and database adjacency records. Reloading two different full molecular
identities under the same record label raises ``SpeciesIdentityError`` when
either record carries a resolved state; ground-only duplicate handling, distinct
labels and repeated identical records retain their existing behavior.

Molecule and reaction drawings cannot represent resolved state identity and raise
``SpeciesIdentityError`` before creating or replacing an image. QM geometry files
and Gaussian/Mopac inputs retain the augmented InChI in their title and use
state-specific file keys. Ground-only kinetics library coverage keys retain
their legacy indexed spelling, such as ``X(3)``.

**Resolved nonthermal reactions.**

    Electron-temperature- or density-dependent reactions involving resolved
    species require independent irreversible forward and reverse declarations.
    See :doc:`plasma_reversal` for import/export policy and its source census.

.. toctree::
    :hidden:

    plasma_reversal
