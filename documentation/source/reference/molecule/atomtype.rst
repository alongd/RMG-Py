***********************
rmgpy.molecule.AtomType
***********************

.. _atom-types:

.. autoclass:: rmgpy.molecule.AtomType
    :members:

.. autofunction:: rmgpy.molecule.getAtomType

The atom type of an atom describes the atom itself and (often) something about
the local bond structure around that atom. This is a useful semantic tool for
accelerating graph isomorphism queries, and a useful shorthand when specifying
molecular substructure patterns via an RMG-style adjacency list.

We define the following basic atom types:

=============== ============================================================
Atom type       Description
=============== ============================================================
*General atom types*
----------------------------------------------------------------------------
``R``           any atom with any local bond structure
``R!H``         any non-hydrogen atom with any local bond structure
*Hydrogen atom types*
----------------------------------------------------------------------------
``H``           hydrogen atom with one single bond
*Carbon atom types*
----------------------------------------------------------------------------
``C``           carbon atom with any local bond structure
``Ca``          carbon atom with two lone pairs and no bonds
``Cs``          carbon atom with up to four single bonds
``Csc``         charged (-1) carbon atom with up to five single bonds
``Cd``          carbon atom with one double bond (not to O or S) and up to two single bonds
``CO``          carbon atom with one double bond to oxygen and up to two single bonds
``CS``          carbon atom with one double bond to sulfur and up to two single bonds
``Cdc``         charged (-1) carbon atom with one double bond and up to three single bonds
``Cdd``         carbon atom with two double bonds
``Ct``          carbon atom with one triple bond and up to one single bond
``Ctc``         charged (-1) carbon atom with one triple bond and up to two single bonds
``Cb``          carbon atom with up to two benzene bonds and up to one single bond
``Cbf``         carbon atom with three benzene bonds
``C2s``         carbon atom with one lone pair (valance 2) and up to two single bonds
``C2sc``        charged (-1) carbon atom with one lone pair (valance 2) and up to three single bonds
``C2d``         carbon atom with one lone pair (valance 2) and one double bond
``C2dc``        charged (-1) carbon atom with one lone pair (valance 2), one double bond and up to one single bond
``C2tc``        charged (-1) carbon atom with one lone pair (valance 2) and one triple bond
*Nitrogen atom types*
----------------------------------------------------------------------------
``N``           nitrogen atom with any local bond structure
``N0sc``        charged (-2) nitrogen atom with three lone pairs (valance 0) with up to one single bond
``N1s``         nitrogen atom with two lone pairs (valance 1) and up to one single bond
``N1sc``        charged (-1) nitrogen atom with two lone pairs (valance 1) up to two single bonds
``N1dc``        charged (-1) nitrogen atom with two lone pairs (valance 1) and one double bond
``N3s``         nitrogen atom with one lone pair (valance 3) with up to three single bonds
``N3d``         nitrogen atom with one lone pair (valance 3), one double bond and up to one single bond
``N3t``         nitrogen atom with one lone pair (valance 3) and one triple bond
``N3b``         nitrogen atom with one lone pair (valance 3) and two benzene bonds
``N5sc``        charged (+1 or +2) nitrogen atom with no lone pairs (valance 5) with up to four single bonds
``N5dc``        charged (+1) nitrogen atom with no lone pairs (valance 5), one double bond and up to two single bonds
``N5dc2``       charged (+2) nitrogen atom with no lone pairs (valance 5), one double bond and up to one single bond
``N5dd``        nitrogen atom with with no lone pairs (valance 5), two double bonds and up to one single bond
``N5ddc``       charged (+1) nitrogen atom with with no lone pairs (valance 5) and two double bonds
``N5tc``        charged (+1) nitrogen atom with with no lone pairs (valance 5), one triple bond and up to one single bond
``N5td``        nitrogen atom with with no lone pairs (valance 5), one triple bond and one double bond
``N5b``         nitrogen atom with with no lone pairs (valance 5) and 2 benzene bonds (one of the lone pairs also participates in the aromatic bond) and up to one single bond
*Oxygen atom types*
----------------------------------------------------------------------------
``O``           oxygen atom with any local bond structure
``Oa``          oxygen atom with three lone pairs and no bonds
``O0sc``        charged (-1) oxygen with three lone pairs (valance 0) and up to one single bond
``O0dc``        charged (-2) oxygen atom with three lone pairs (valance 0) and one double bond
``O2s``         oxygen atom with two lone pairs (valance 2) and up to two single bonds
``O2sp``        charged (+1) oxygen atom with two lone pairs (valance 2) and up to one single bond
``O2sn``        charged (-1) oxygen atom with two lone pairs (valance 2) and up to three single bonds
``O2d``         oxygen atom with two lone pairs (valance 2) and one doubel bond
``O2dc``        charged (-1) oxygen atom with two lone pairs (valance 2), one double bond and up to one single bond
``O2tc``        charged (-1) oxygen atom with two lone pairs (valance 2) and one triple bond
``O4sc``        charged (+1) oxygen atom with one one pair (valance 4) and up to three single bonds
``O4dc``        charged (+1) oxygen atom with one one pair (valance 4), one double bond and up to one single bond
``O4tc``        charged (+1) oxygen atom with one one pair (valance 4) and one triple bond
*Silicon atom types*
----------------------------------------------------------------------------
``Si``          silicon atom with any local bond structure
``Sis``         silicon atom with four single bonds
``Sid``         silicon atom with one double bond (to carbon) and two single bonds
``SiO``         silicon atom with one double bond (to oxygen) and two single bonds
``Sidd``        silicon atom with two double bonds
``Sit``         silicon atom with one triple bond and one single bond
``Sib``         silicon atom with two benzene bonds and one single bond
``Sibf``        silicon atom with three benzene bonds
*Sulfur atom types*
----------------------------------------------------------------------------
``S``           sulfur atom with any local bond structure
``Sa``          sulfur atom with three lone pairs and no bonds
``S0s``         charged (-1) sulfur atom with three lone pairs (valance 0) and up to one single bond
``S2s``         sulfur atom with two lone pairs (valance 2) and up to two single bonds
``S2sp``        charged (+1) sulfur atom with two lone pairs (valance 2) and up to one single bond
``S2sn``        charged (-1) sulfur atom with two lone pairs (valance 2) and up to three single bonds
``S2d``         sulfur atom with two lone pairs (valance 2) and one double bond
``S2dc``        charged (-1) sulfur atom with two lone pairs (valance 2), one double bond and up to one single bond
``S4s``         sulfur atom with one lone pair (valance 4) and up to four single bonds
``S4sc``        charged (+1) sulfur atom with one lone pair (valance 4) and up to three single bonds
``S4d``         sulfur atom with one lone pair (valance 4), one double bond and up to two single bonds
``S4dc``        charged (+1) sulfur atom with one lone pair (valance 4), one double bond and up to one single bond
``S4b``         sulfur atom with one lone pair (valance 4) and two benzene bonds (one of the lone pairs also participates in the aromatic bond)
``S4dd``        sulfur atom with one lone pair (valance 4) and two double bonds
``S4ddc``       charged (-1) sulfur atom with one lone pair (valance 4), two double bonds and up to one single bond
``S4t``         sulfur atom with one lone pair (valance 4), one triple bond and up to one single bond
``S4tc``        charged (+1) sulfur atom with one lone pair (valance 4) and one triple bond
``S6s``         sulfur atom with no lone pairs (valance 6) and up to six single bonds
``S6d``         sulfur atom with no lone pairs (valance 6), one double bond and up to four single bonds
``S6dc``        charged (+1 or +2) sulfur atom with no lone pairs (valance 6), one double bond and up to three single bonds
``S6dd``        sulfur atom with no lone pairs (valance 6), two double bonds and up to two single bonds
``S6ddc``       charged (+1) sulfur atom with no lone pairs (valance 6), two double bonds and up to one single bond
``S6ddd``       sulfur atom with no lone pairs (valance 6) and three double bonds
``S6t``         sulfur atom with no lone pairs (valance 6), one triple bond and up to three single bonds
``S6tc``        charged (+1) sulfur atom with no lone pairs (valance 6), one triple bond and up to two single bonds
``S6td``        sulfur atom with no lone pairs (valance 6), one triple bond, one double bond and up to one single bond
``S6tdc``       charged (+1) sulfur atom with no lone pairs (valance 6), one triple bond and one double bond
``S6tt``        sulfur atom with no lone pairs (valance 6) and two triple bonds
=============== ============================================================
