#!/usr/bin/env python
# encoding: utf-8

################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2017 Prof. William H. Green (whgreen@mit.edu), 
#   Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
This module defines the atom types that are available for representing
molecular functional groups and substructure patterns. Each available atom type
is defined as an instance of the :class:`AtomType` class. The atom types 
themselves are available in the ``atomTypes`` module-level variable, or as
the return value from the :meth:`getAtomType()` method.

If you want to change which atom types are available in RMG and/or what they
represent, this should be the only module you need to change to do so.
"""

import cython
from rmgpy.exceptions import AtomTypeError

################################################################################

class AtomType:
    """
    A class for internal representation of atom types. Using unique objects
    rather than strings allows us to use fast pointer comparisons instead of
    slow string comparisons, as well as store extra metadata. In particular,
    we store metadata describing the atom type's hierarchy with regard to other
    atom types, and the atom types that can result when various actions
    involving this atom type are taken. The attributes are:

    =================== =================== ====================================
    Attribute           Type                Description
    =================== =================== ====================================
    `label`             ``str``             A unique label for the atom type
    `generic`           ``list``            The atom types that are more generic than this one
    `specific`          ``list``            The atom types that are more specific than this one
    `incrementBond`     ``list``            The atom type(s) that result when an adjacent bond's order is incremented
    `decrementBond`     ``list``            The atom type(s) that result when an adjacent bond's order is decremented
    `formBond`          ``list``            The atom type(s) that result when a new single bond is formed to this atom type
    `breakBond`         ``list``            The atom type(s) that result when an existing single bond to this atom type is broken
    `incrementRadical`  ``list``            The atom type(s) that result when the number of radical electrons is incremented
    `decrementRadical`  ``list``            The atom type(s) that result when the number of radical electrons is decremented
    `incrementLonePair` ``list``            The atom type(s) that result when the number of lone electron pairs is incremented
    `decrementLonePair` ``list``            The atom type(s) that result when the number of lone electron pairs is decremented

    The following features are what are required in a given atomtype. Any int in the list is acceptable.
    An empty list is a wildcard
    'single'            ''list''            The total number of single bonds on the atom
    'allDouble'         ''list''            The total number of double bonds on the atom
    'rDouble'           ''list''            The number of double bonds to any non-oxygen, nonsulfur
    'oDouble'           ''list''            The number of double bonds to oxygen
    'sDouble'           ''list''            The number of double bonds to sulfur
    'triple'            ''list''            The total number of triple bonds on the atom
    'benzene'           ''list''            The total number of benzene bonds on the atom
    'lonePairs'         ''list''            The number of lone pairs on the atom
    'charge'            ''list''            The partial charge of the atom
    =================== =================== ====================================

    """
    

    def __init__(self, label='', generic=None, specific=None,
                 single=None,
                 allDouble=None,
                 rDouble=None,
                 oDouble=None,
                 sDouble=None,
                 triple=None,
                 benzene=None,
                 lonePairs=None,
                 charge=None):
        self.label = label
        self.generic = generic or []
        self.specific = specific or []
        self.incrementBond = []
        self.decrementBond = []
        self.formBond = []
        self.breakBond = []
        self.incrementRadical = []
        self.decrementRadical = []
        self.incrementLonePair = []
        self.decrementLonePair = []
        self.single = single or []
        self.allDouble = allDouble or []
        self.rDouble = rDouble or []
        self.oDouble = oDouble or []
        self.sDouble = sDouble or []
        self.triple = triple or []
        self.benzene = benzene or []
        self.lonePairs = lonePairs or []
        self.charge = charge or []

    def __repr__(self):
        return '<AtomType "%s">' % self.label

    def __reduce__(self):
        """
        A helper function used when pickling an AtomType object.
        """
        d = {
            'label': self.label,
            'generic': self.generic,
            'specific': self.specific,
            'incrementBond': self.incrementBond,
            'decrementBond': self.decrementBond,
            'formBond': self.formBond,
            'breakBond': self.breakBond,
            'incrementRadical': self.incrementRadical,
            'decrementRadical': self.decrementRadical,
            'incrementLonePair': self.incrementLonePair,
            'decrementLonePair': self.decrementLonePair,
            'single': self.single,
            'allDouble': self.allDouble,
            'rDouble': self.rDouble,
            'oDouble': self.oDouble,
            'sDouble': self.sDouble,
            'triple': self.triple,
            'benzene': self.benzene,
            'lonePairs': self.lonePairs,
            'charge': self.charge
        }
        return (AtomType, (), d)

    def __setstate__(self, d):
        """
        A helper function used when unpickling an AtomType object.
        """
        self.label = d['label']
        self.generic = d['generic']
        self.specific = d['specific']
        self.incrementBond = d['incrementBond']
        self.decrementBond = d['decrementBond']
        self.formBond = d['formBond']
        self.breakBond = d['breakBond']
        self.incrementRadical = d['incrementRadical']
        self.decrementRadical = d['decrementRadical']
        self.incrementLonePair = d['incrementLonePair']
        self.decrementLonePair = d['decrementLonePair']
        self.single = d['single']
        self.allDouble = d['allDouble']
        self.rDouble = d['rDouble']
        self.oDouble = d['oDouble']
        self.sDouble = d['sDouble']
        self.triple = d['triple']
        self.benzene = d['benzene']
        self.lonePairs = d['lonePairs']
        self.charge = d['charge']

    def setActions(self, incrementBond, decrementBond, formBond, breakBond, incrementRadical, decrementRadical, incrementLonePair, decrementLonePair):
        self.incrementBond = incrementBond
        self.decrementBond = decrementBond
        self.formBond = formBond
        self.breakBond = breakBond
        self.incrementRadical = incrementRadical
        self.decrementRadical = decrementRadical
        self.incrementLonePair = incrementLonePair
        self.decrementLonePair = decrementLonePair

    def equivalent(self, other):
        """
        Returns ``True`` if two atom types `atomType1` and `atomType2` are
        equivalent or ``False``  otherwise. This function respects wildcards,
        e.g. ``R!H`` is equivalent to ``C``.
        """
        return self is other or self in other.specific or other in self.specific

    def isSpecificCaseOf(self, other):
        """
        Returns ``True`` if atom type `atomType1` is a specific case of
        atom type `atomType2` or ``False``  otherwise.
        """
        return self is other or self in other.specific
    
    def getFeatures(self):
        """
        Returns a list of the features that are checked to determine atomtype
        """
        features=[self.single,
                  self.allDouble,
                  self.rDouble,
                  self.oDouble,
                  self.sDouble,
                  self.triple,
                  self.benzene,
                  self.lonePairs,
                  self.charge]
        return features

################################################################################


"""
Note: function to read adjacency lists assumes that all atom types begin
with a capital letter [A-Z]

For making sample atoms, we use the first atomtype under specific,
therefore the first one in the list should always be an element.

The atomTypes naming convention is:
<element> <valence> <characteristic bonds> <charge(optional)>
For example:
- N3d is nitrogen with valence=3 (i.e., 3 electronce are able to form bonds or remain as radicals) with one double bond
- S2tc is a charged sulful with valence=2 with a triple bonds
- Oa is atomic oxygen, i.e., a closed shell atom
- O2sn is a negatively charged oxygen with valence=2
Some charged atom types were merged together, and are marked as '*Composite atomType'
"""

atomTypes = {}
atomTypes['R']    = AtomType(label='R', generic=[], specific=[
    'H',
    'R!H',
    'Val4','Val5','Val6','Val7',    
    'He',
    'C','Ca','Cs','Csc','Cd','CO','CS','Cdd','Cdc','Ct','Ctc','Cb','Cbf','C2s','C2sc','C2d','C2dc','C2tc',
    'N','N0sc','N0dc','N1s','N1sc','N1dc','N1tc','N3s','N3d','N3dc','N3t','N3tdc','N3b','N5sc','N5d','N5dc','N5dd','N5ddc','N5dddc','N5t','N5tc','N5td','N5tdc','N5b',
    'O','Oa','O0sc','O0dc','O2s','O2sp','O2sn','O2d','O2dc','O2tc','O4sc','O4dc','O4tc',
    'Ne',
    'Si','Sis','Sid','Sidd','Sit','SiO','Sib','Sibf',
    'S','Sa','S0c','S2s','S2sc','S2d','S2dc','S2tc','S4s','S4sc','S4d','S4dd','S4dc','S4b','S4t','S4tdc','S6s','S6sc','S6d','S6dd','S6ddd','S6dc','S6t','S6td','S6tt','S6tdc',
    'Cl','Ar']
)
atomTypes['R!H']  = AtomType(label='R!H', generic=['R'], specific=[
    'He',
    'Val4','Val5','Val6','Val7',
    'C','Ca','Cs','Csc','Cd','CO','CS','Cdd','Cdc','Ct','Ctc','Cb','Cbf','C2s','C2sc','C2d','C2dc','C2tc',
    'N','N0sc','N0dc','N1s','N1sc','N1dc','N1tc','N3s','N3d','N3dc','N3t','N3tdc','N3b','N5sc','N5d','N5dc','N5dd','N5ddc','N5dddc','N5t','N5tc','N5td','N5tdc','N5b',
    'O','Oa','O0sc','O0dc','O2s','O2sp','O2sn','O2d','O2dc','O2tc','O4sc','O4dc','O4tc',
    'Ne',
    'Si','Sis','Sid','Sidd','Sit','SiO','Sib','Sibf',
    'S','Sa','S0c','S2s','S2sc','S2d','S2dc','S2tc','S4s','S4sc','S4d','S4dd','S4dc','S4b','S4t','S4tdc','S6s','S6sc','S6d','S6dd','S6ddd','S6dc','S6t','S6td','S6tt','S6tdc',
    'Cl','Ar'])

atomTypes['Val4'] = AtomType(label='Val4', generic=['R','R!H'], specific=[
    'C','Ca','Cs','Csc','Cd','CO','CS','Cdd','Cdc','Ct','Ctc','Cb','Cbf','C2s','C2sc','C2d','C2dc','C2tc',
    'Si','Sis','Sid','Sidd','Sit','SiO','Sib','Sibf'])

atomTypes['Val5'] = AtomType(label='Val5', generic=['R','R!H'], specific=[
    'N','N0sc','N0dc','N1s','N1sc','N1dc','N1tc','N3s','N3d','N3dc','N3t','N3tdc','N3b','N5sc','N5d','N5dc','N5dd','N5ddc','N5dddc','N5t','N5tc','N5td','N5tdc','N5b'])

atomTypes['Val6'] = AtomType(label='Val6', generic=['R','R!H'], specific=[
    'O','Oa','O0sc','O0dc','O2s','O2sp','O2sn','O2d','O2dc','O2tc','O4sc','O4dc','O4tc',
    'S','Sa','S0c','S2s','S2sc','S2d','S2dc','S2tc','S4s','S4sc','S4d','S4dd','S4dc','S4b','S4t','S4tdc','S6s','S6sc','S6d','S6dd','S6ddd','S6dc','S6t','S6td','S6tt','S6tdc'])

atomTypes['Val7'] = AtomType(label='Val7', generic=['R','R!H'], specific=[
    'Cl'])

atomTypes['H'   ] = AtomType('H',    generic=['R'],            specific=[])

atomTypes['He'   ] = AtomType('He',  generic=['R','R!H'],      specific=[])

atomTypes['C'   ] = AtomType('C',    generic=['R','R!H','Val4'],      specific=['Ca','Cs','Csc','Cd','CO','CS','Cdd','Cdc','Ct','Ctc','Cb','Cbf','C2s','C2sc','C2d','C2dc','C2tc'],
                             single=[], allDouble=[], rDouble=[], oDouble=[], sDouble=[], triple=[], benzene=[], lonePairs=[], charge=[])
atomTypes['Ca'  ] = AtomType('Ca',   generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[0])
# examples for Ca: atomic carbon (closed shell)
atomTypes['Cs'  ] = AtomType('Cs',   generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0,1,2,3,4], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for Cs: C, CC, C1=CCC([O-])[CH+]1
atomTypes['Csc' ] = AtomType('Csc',  generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0,1,2,3,4,5,6], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[-2,-1,+1,+2])
# examples for Csc: [CH4-][N+]#N
atomTypes['Cd'  ] = AtomType('Cd',   generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0,1,2], allDouble=[1], rDouble=[1], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[0], charge=[])
# examples for Cd: C=C, C=N
atomTypes['CO'  ] = AtomType('CO',   generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0,1,2], allDouble=[1], rDouble=[0], oDouble=[1], sDouble=[0], triple=[0], benzene=[0], lonePairs=[0], charge=[])
# examples for CO: C=O
atomTypes['CS'  ] = AtomType('CS',   generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0,1,2], allDouble=[1], rDouble=[0], oDouble=[0], sDouble=[1], triple=[0], benzene=[0], lonePairs=[0], charge=[])
# examples for CS: C=S
atomTypes['Cdd' ] = AtomType('Cdd',  generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0], allDouble=[2], rDouble=[0,1,2], oDouble=[0,1,2], sDouble=[0,1,2], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for Cdd: O=C=O, C=C=C
atomTypes['Cdc' ] = AtomType('Cdc',  generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0,1,2,3,4], allDouble=[1,2], rDouble=[0,1,2], oDouble=[0,1,2], sDouble=[0,1,2], triple=[0], benzene=[0], lonePairs=[0], charge=[-2,-1])
# examples for Cdc: O=[CH2-][N+]#N
atomTypes['Ct'  ] = AtomType('Ct',   generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[0], charge=[])
# examples for Ct: C#C, C#N
atomTypes['Ctc' ] = AtomType('Ctc',  generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0,1,2,3], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[0], charge=[-2,-1])
# examples for Ctc: [CH2-]#[NH+]
atomTypes['Cb'  ] = AtomType('Cb',   generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[1,2], lonePairs=[], charge=[])
# examples for Cb: benzene (C6H6)
atomTypes['Cbf' ] = AtomType('Cbf',  generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[3], lonePairs=[], charge=[])
# examples for Cbf: Naphthalene
atomTypes['C2s' ] = AtomType('C2s',  generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0,1,2], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[0])
# examples for C2s: singlet[CH2]
atomTypes['C2sc'] = AtomType('C2sc', generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0,1,2,3,4], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[-2,-1])
# examples for C2sc: [CH2-][N+]#N
atomTypes['C2d' ] = AtomType('C2d',  generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[0])
# examples for C2d: singlet[C]=C
atomTypes['C2dc'] = AtomType('C2dc', generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0,1,2], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[-2,-1])
# examples for C2dc: C=[C-][N+]#N, [CH-]=[N+]=O
atomTypes['C2tc' ] = AtomType('C2tc',  generic=['R','R!H','C','Val4'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[1], charge=[-2,-1])
# examples for C2tc: [C-]#[O+]

atomTypes['N'   ] = AtomType('N',    generic=['R','R!H','Val5'],      specific=['N0sc','N0dc','N1s','N1sc','N1dc','N1tc','N3s','N3d','N3dc','N3t','N3tdc','N3b','N5sc','N5d','N5dc','N5dd','N5ddc','N5dddc','N5t','N5tc','N5td','N5tdc','N5b'],
                             single=[], allDouble=[], rDouble=[], oDouble=[], sDouble=[], triple=[], benzene=[], lonePairs=[], charge=[])
atomTypes['N0sc'] = AtomType('N0sc',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1,2], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[3], charge=[-3,-2])
# examples for N0sc: [NH+]#[N+][N-2] with adjList 1 N u0 p0 c+1 {2,S} {3,T}; 2 H u0 p0 c0 {1,S}; 3 N u0 p0 c+1 {1,T} {4,S}; 4 N u0 p3 c-2 {3,S}
atomTypes['N0dc'] = AtomType('N0dc',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[3], charge=[-3])
# examples for N0dc: [O-]=[N+]=[O+][N-3]#[N+2]
atomTypes['N1s' ] = AtomType('N1s',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[2], charge=[0])
# examples for N1s: closed shell N-N, closed shell NH
atomTypes['N1sc' ] = AtomType('N1sc',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1,2,3,4], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[-3,-2,-1])
# examples for N1sc: [NH-][S+]=C, [NH-][N+]#C
atomTypes['N1dc' ] = AtomType('N1dc',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1,2], allDouble=[1,2], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[-3,-2,-1])
# examples for N1dc: [N-]=[N+]=N terminal nitrogen on azide (two lone pairs), [N-]=[NH+]
atomTypes['N1tc' ] = AtomType('N1tc',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[1], benzene=[0], lonePairs=[2], charge=[-3,-2])
atomTypes['N3s' ] = AtomType('N3s',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1,2,3], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[1], charge=[])
# examples for N3s: NH3, NH2, NH, N, C[NH]...
atomTypes['N3d' ] = AtomType('N3d',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[0])
# examples for N3d: N=O, N=N, C=N, [O]N=O, [N]=O, [N]=C
atomTypes['N3dc'] = AtomType('N3dc', generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1,2,3,4], allDouble=[1,2], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[-3,-2,-1,+1])
# *Composite atomType; examples for N3dc: [NH-]=[SH+], [N+]=[S-]=O, C=[N-]=[SH+]
atomTypes['N3t' ] = AtomType('N3t',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[1], benzene=[0], lonePairs=[1], charge=[0])
# examples for N3t: N2, N#C, N#[C], N#CC
atomTypes['N3tdc'] = AtomType('N3tdc',generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1,2,3], allDouble=[0,1], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[1], charge=[-3,-2,-1])
# *Composite atomType; examples for N3tdc: [NH-]#[S+], [S+]#[N-2]=[SH+]
atomTypes['N3b' ] = AtomType('N3b',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[2], lonePairs=[1], charge=[0])
# examples for N3b: Oxazole, Pyradine, Pyrazine, 1,3,5-Triazine, Benzimidazole, Purine
atomTypes['N5sc' ] = AtomType('N5sc',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1,2,3,4], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[0], charge=[+1,+2,+3])
# examples for N5sc: [NH4+], [NH2+][O-] {N has u1 p0}, [NH3+][O-]
atomTypes['N5d'  ] = AtomType('N5d',   generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1,2,3], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for N5d: [NH2]=N
atomTypes['N5dc' ] = AtomType('N5dc',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1,2], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[+1,+2,+3])
# examples for N5dc: O[N+](=O)(O-) nitrate group, [N+](=O)(O)[O-], O=[N+][O-], [N+](=O)(O[N+](=O)[O-])[O-], [O-2]=[NH+2]
atomTypes['N5dd' ] = AtomType('N5dd',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1], allDouble=[2], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for N5dd: ON(=O)=O
atomTypes['N5ddc'] = AtomType('N5ddc', generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0], allDouble=[2], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[+1])
# examples for N5ddc: N=[N+]=[N-] center nitrogen on azide, [N-]=[N+]=O
atomTypes['N5dddc'] = AtomType('N5dddc', generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0], allDouble=[3], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[-1])
atomTypes['N5t'  ] = AtomType('N5t',   generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1,2], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[1], benzene=[0], lonePairs=[0], charge=[0])
# examples for N5t: C#[NH2]
atomTypes['N5tc' ] = AtomType('N5tc',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[1], benzene=[0], lonePairs=[0], charge=[+1,+2])
# examples for N5tc: C[N+]#[C-] isocyano group, N#[N+][O-], [NH+]#[C-] (note that C- has p1 here), [N+]#[C-] (note that C- has p1 here), [O-][N+]#C, C[N+]#[C-] (note that C- has p1 here)
atomTypes['N5td' ] = AtomType('N5td',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[0], charge=[0])
# examples for N5td: N#N=O
atomTypes['N5tdc'] = AtomType('N5tdc', generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1,2], allDouble=[0,1,2], rDouble=[], oDouble=[], sDouble=[], triple=[1,2], benzene=[0], lonePairs=[0], charge=[-3,-2,-1])
# *Composite atomType; examples for N5tdc: [NH+]#[NH-]=O
atomTypes['N5b' ] = AtomType('N5b',  generic=['R','R!H','N','Val5'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[2], lonePairs=[0], charge=[0,+1])
# examples for N5b: Pyrrole, Indole, Benzimidazole, Purine; Note that this is the only N atomType with valence 5 which isn't necessarily charged.

atomTypes['O'   ] = AtomType('O',    generic=['R','R!H','Val6'],      specific=['Oa','O0sc','O0dc','O2s','O2sp','O2sn','O2d','O2dc','O2tc','O4sc','O4dc','O4tc'],
                             single=[], allDouble=[], rDouble=[], oDouble=[], sDouble=[], triple=[], benzene=[], lonePairs=[], charge=[])
atomTypes['Oa'  ] = AtomType('Oa',   generic=['R','R!H','O','Val6'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[3], charge=[0])
# examples for Oa: atomic oxygen (closed shell)
atomTypes['O0sc'] = AtomType('O0sc', generic=['R','R!H','O','Val6'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[3], charge=[-2,-1])
# examples for O0sc: Nitric acid O[N+](=O)([O-])
atomTypes['O0dc'] = AtomType('O0dc', generic=['R','R!H','O','Val6'],  specific=[],
                             single=[0], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[3], charge=[-2])
# examples for O0dc: [O-2][NH+2]
atomTypes['O2s' ] = AtomType('O2s',  generic=['R','R!H','O','Val6'],  specific=[],
                             single=[0,1,2], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[0])
# examples for O2s: H2O, OH, CH3OH
atomTypes['O2sp'] = AtomType('O2sp', generic=['R','R!H','O','Val6'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[+1])
# examples for O2sp: C=[S-][O+]
atomTypes['O2sn'] = AtomType('O2sn', generic=['R','R!H','O','Val6'],  specific=[],
                             single=[0,1,2,3], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[-1])
# examples for O2sn: C=[NH+][OH2-]
atomTypes['O2d' ] = AtomType('O2d',  generic=['R','R!H','O','Val6'],  specific=[],
                             single=[0], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[0])
# examples for O2d: CO2, CH2O
atomTypes['O2dc'] = AtomType('O2dc', generic=['R','R!H','O','Val6'],  specific=[],
                             single=[0,1], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[-1])
# examples for O2dc: O=[NH+][O-]=C
atomTypes['O2tc'] = AtomType('O2tc', generic=['R','R!H','O','Val6'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[2], charge=[-1])
# examples for O2tc: [O-]#[NH+]
atomTypes['O4sc'] = AtomType('O4sc', generic=['R','R!H','O','Val6'],  specific=[],
                             single=[0,1,2,3], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[+1,+2,+3])
# examples for O4sc: [O-][OH+]C
atomTypes['O4dc'] = AtomType('O4dc', generic=['R','R!H','O','Val6'],  specific=[],
                             single=[0,1], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[+1,+2])
# examples for O4dc: the positively charged O in ozone [O-][O+]=O
atomTypes['O4tc' ] = AtomType('O4tc',  generic=['R','R!H','O','Val6'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[1], charge=[+1])
# examples for O4tc: [C-]#[O+]

atomTypes['Ne'  ] = AtomType('Ne',   generic=['R','R!H'],      specific=[])
atomTypes['Si'  ] = AtomType('Si',   generic=['R','R!H','Val4'],      specific=['Sis','Sid','Sidd','Sit','SiO','Sib','Sibf'],
                             single=[], allDouble=[], rDouble=[], oDouble=[], sDouble=[], triple=[], benzene=[], lonePairs=[], charge=[])
atomTypes['Sis' ] = AtomType('Sis',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[], charge=[])
atomTypes['SiO' ] = AtomType('SiO',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[1], rDouble=[], oDouble=[1], sDouble=[], triple=[0], benzene=[0], lonePairs=[], charge=[])
atomTypes['Sid' ] = AtomType('Sid',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[1], rDouble=[], oDouble=[0], sDouble=[], triple=[0], benzene=[0], lonePairs=[], charge=[])
atomTypes['Sidd'] = AtomType('Sidd', generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[2], rDouble=[0,1,2], oDouble=[0,1,2], sDouble=[0,1,2], triple=[0], benzene=[0], lonePairs=[], charge=[])
atomTypes['Sit' ] = AtomType('Sit',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[], charge=[])
atomTypes['Sib' ] = AtomType('Sib',  generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[2], lonePairs=[], charge=[])
atomTypes['Sibf'] = AtomType('Sibf', generic=['R','R!H','Si','Val4'], specific=[],
                             single=[], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[3], lonePairs=[], charge=[])

atomTypes['S'   ] = AtomType('S',    generic=['R','R!H','Val6'],      specific=['Sa','S0c','S2s','S2sc','S2d','S2dc','S2tc','S4s','S4sc','S4d','S4dd','S4dc','S4b','S4t','S4tdc','S6s','S6sc','S6d','S6dd','S6ddd','S6dc','S6t','S6td','S6tt','S6tdc'],
                             single=[], allDouble=[], rDouble=[], oDouble=[], sDouble=[], triple=[], benzene=[], lonePairs=[], charge=[])
atomTypes['Sa'  ] = AtomType('Sa',   generic=['R','R!H','S','Val6'],  specific=[],
                            single=[0], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[3], charge=[0])
# examples for Sa: atomic sulfur (closed shell)
atomTypes['S0c'] = AtomType('S0c', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2], allDouble=[0,1], rDouble=[], oDouble=[], sDouble=[], triple=[0,1], benzene=[0], lonePairs=[3], charge=[-3,-2,-1])
# examples for S0c: [S-][S+], S3 {[S-][S+]=S}
atomTypes['S2s' ] = AtomType('S2s',  generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[2], charge=[0])
# examples for S2s: [S], [SH], S {H2S}, [S][S], SS {H2S2}, SSC, CSSC, SO {HSOH}...
atomTypes['S2sc'] = AtomType('S2sc', generic=['R','R!H','S','Val6'],  specific=[],
                            single=[0,1,2,3,4,5], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[2], charge=[-3,-2,-1,+1])
# examples for S2sc: [S-][S+], N#[N+][S-](O)O
atomTypes['S2d' ] = AtomType('S2d',  generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[0])
# examples for S2d: S=S, C=S, S=O, S=N, S=C=S, S=C=O, S=C=S...
atomTypes['S2dc'] = AtomType('S2dc', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2,3], allDouble=[1,2], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[2], charge=[-3,-2,-1])
# *Composite atomType; examples for S2dc: [SH-]=[N+]
atomTypes['S2tc'] = AtomType('S2tc', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[2], charge=[-3,-2,-1])
# examples for S2tc: [S-]#[NH+]
atomTypes['S4s' ] = AtomType('S4s',  generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2,3,4], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[1], charge=[0])
# examples for S4s: H4S, SH3CH3...
atomTypes['S4sc'] = AtomType('S4sc', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2,3], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[1], charge=[+1,+2,+3])
# examples for S4sc: CS[S+]([O-])C, [O-][S+2][O-]
atomTypes['S4d' ] = AtomType('S4d',  generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[0])
# examples for S4d: O=S(O)O {Sulfurous acid}
atomTypes['S4dd'] = AtomType('S4dd', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0], allDouble=[2], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[0])
# examples for S4dd: O=S=O
atomTypes['S4dc'] = AtomType('S4dc', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2,3,4], allDouble=[1,2,3], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[1], charge=[-3,-2,-1,+1,+2])
# *Composite atomType; examples for S4dc: [CH2-][S+]=C {where the [CH2-] has a lone pair}, [O+][S-](=O)=O, [O-][S+]=C, [NH-][S+]=C {where the [NH-] has two lone pairs}, [O-][S+]=O
atomTypes['S4b' ] = AtomType('S4b',  generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[2], lonePairs=[1], charge=[0])
# examples for S4b: Thiophene, Benzothiophene, Benzo[c]thiophene, Thiazole, Benzothiazole...
atomTypes['S4t' ] = AtomType('S4t',  generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[1], benzene=[0], lonePairs=[1], charge=[0])
# examples for S4t: C#S, C#SO, C#[S]
atomTypes['S4tdc'] = AtomType('S4tdc',generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2,3], allDouble=[0,1,2], rDouble=[], oDouble=[], sDouble=[], triple=[1,2], benzene=[0], lonePairs=[1], charge=[-3,-2,-1,+1,+2])
# *Composite atomType; examples for S4tdc: [C-]#[S+]
atomTypes['S6s' ] = AtomType('S6s',  generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2,3,4,5,6], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for S6s: H6S, F6S
atomTypes['S6sc'] = AtomType('S6sc', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2,3,4,5], allDouble=[0], rDouble=[0], oDouble=[0], sDouble=[0], triple=[0], benzene=[0], lonePairs=[0], charge=[+1,+2,+3])
# examples for S6sc: O[SH3+][O-], [O-][S+3]([O-])[O-]
atomTypes['S6d' ] = AtomType('S6d',  generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2,3,4], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for S6d: [SH4]=O, SF4=O, [SH4]=C, C[SH3]=C...
atomTypes['S6dd'] = AtomType('S6dd', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2], allDouble=[2], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for S6dd: S(=O)(=O)(O)O {H2SO4, Sulfuric acid}, Perfluorooctanesulfonic acid, Pyrosulfuric acid, Thiosulfuric acid {middle S}, OS(=O)(=O)OOS(=O)(=O)O
atomTypes['S6ddd'] = AtomType('S6ddd', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0], allDouble=[3], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[0])
# examples for S6ddd: O=S(=O)(=O)
atomTypes['S6dc'] = AtomType('S6dc', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2,3,4,5], allDouble=[1,2,3], rDouble=[], oDouble=[], sDouble=[], triple=[0], benzene=[0], lonePairs=[0], charge=[-3,-2,-1,+1,+2,+3])
# *Composite atomType; examples for S6dc: O=[S+2]([O-])[O-], [CH-]=[SH3+], [CH-]=[SH2+]O, [CH-][SH2+], O=[S+](=O)[O-], [OH+]=[S-](=O)=O
atomTypes['S6t'  ] = AtomType('S6t', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2,3], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[0], charge=[0])
# examples for S6t: H3S#N
atomTypes['S6td' ] = AtomType('S6td', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1], allDouble=[1], rDouble=[], oDouble=[], sDouble=[], triple=[1], benzene=[0], lonePairs=[0], charge=[0])
# examples for S6td: HS(=O)#N
atomTypes['S6tt' ] = AtomType('S6tt', generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0], allDouble=[0], rDouble=[], oDouble=[], sDouble=[], triple=[2], benzene=[0], lonePairs=[0], charge=[0])
# examples for S6tt: N#S#N
atomTypes['S6tdc'] = AtomType('S6tdc',generic=['R','R!H','S','Val6'],  specific=[],
                             single=[0,1,2,3,4], allDouble=[0,1,2], rDouble=[], oDouble=[], sDouble=[], triple=[1,2], benzene=[0], lonePairs=[0], charge=[-3,-2,-1,+1,+2,+3])
# *Composite atomType; examples for S6tdc: [N-2][N+]#[S+]=O, [SH2+]#[C-]


atomTypes['Cl'  ] = AtomType('Cl',   generic=['R','R!H','Val7'],      specific=[])

atomTypes['Ar'  ] = AtomType('Ar',   generic=['R','R!H'],      specific=[])

atomTypes['R'   ].setActions(incrementBond=['R'],            decrementBond=['R'],            formBond=['R'],         breakBond=['R'],         incrementRadical=['R'],    decrementRadical=['R'],    incrementLonePair=['R'],   decrementLonePair=['R'])
atomTypes['R!H' ].setActions(incrementBond=['R!H'],          decrementBond=['R!H'],          formBond=['R!H'],       breakBond=['R!H'],       incrementRadical=['R!H'],  decrementRadical=['R!H'],  incrementLonePair=['R!H'], decrementLonePair=['R!H'])
atomTypes['Val4'].setActions(incrementBond=['Val4'],         decrementBond=['Val4'],         formBond=['Val4'],      breakBond=['Val4'],      incrementRadical=['Val4'], decrementRadical=['Val4'], incrementLonePair=['Val4'],decrementLonePair=['Val4'])
atomTypes['Val5'].setActions(incrementBond=['Val5'],         decrementBond=['Val5'],         formBond=['Val5'],      breakBond=['Val5'],      incrementRadical=['Val5'], decrementRadical=['Val5'], incrementLonePair=['Val5'],decrementLonePair=['Val5'])
atomTypes['Val6'].setActions(incrementBond=['Val6'],         decrementBond=['Val6'],         formBond=['Val6'],      breakBond=['Val6'],      incrementRadical=['Val6'], decrementRadical=['Val6'], incrementLonePair=['Val6'],decrementLonePair=['Val6'])
atomTypes['Val7'].setActions(incrementBond=['Val7'],         decrementBond=['Val7'],         formBond=['Val7'],      breakBond=['Val7'],      incrementRadical=['Val7'], decrementRadical=['Val7'], incrementLonePair=['Val7'],decrementLonePair=['Val7'])


atomTypes['H'   ].setActions(incrementBond=[],               decrementBond=[],               formBond=['H'],         breakBond=['H'],         incrementRadical=['H'],    decrementRadical=['H'],    incrementLonePair=[],      decrementLonePair=[])

atomTypes['He'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=['He'],   decrementRadical=['He'],   incrementLonePair=[],      decrementLonePair=[])

atomTypes['C'   ].setActions(incrementBond=['C'],            decrementBond=['C'],            formBond=['C'],         breakBond=['C'],         incrementRadical=['C'],    decrementRadical=['C'],    incrementLonePair=['C'],   decrementLonePair=['C'])
atomTypes['Ca'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['C2s'])
atomTypes['Cs'  ].setActions(incrementBond=['Cd','CO','CS'], decrementBond=[],               formBond=['Cs','Csc'],  breakBond=['Cs'],        incrementRadical=['Cs'],   decrementRadical=['Cs'],   incrementLonePair=['C2s'], decrementLonePair=['C2s'])
atomTypes['Csc' ].setActions(incrementBond=['Cdc'],          decrementBond=[],               formBond=['Csc'],       breakBond=['Csc','Cs'],  incrementRadical=['Csc'],  decrementRadical=['Csc'],  incrementLonePair=['C2sc'],decrementLonePair=['C2sc'])
atomTypes['Cd'  ].setActions(incrementBond=['Cdd','Ct','C2tc'],decrementBond=['Cs'],         formBond=['Cd','Cdc'],  breakBond=['Cd'],        incrementRadical=['Cd'],   decrementRadical=['Cd'],   incrementLonePair=['C2d'], decrementLonePair=[])
atomTypes['CO'  ].setActions(incrementBond=['Cdd','C2tc'],   decrementBond=['Cs'],           formBond=['CO','Cdc'],  breakBond=['CO'],        incrementRadical=['CO'],   decrementRadical=['CO'],   incrementLonePair=['C2d'], decrementLonePair=[])
atomTypes['CS'  ].setActions(incrementBond=['Cdd','C2tc'],   decrementBond=['Cs'],           formBond=['CS','Cdc'],  breakBond=['CS'],        incrementRadical=['CS'],   decrementRadical=['CS'],   incrementLonePair=['C2d'], decrementLonePair=[])
atomTypes['Cdc' ].setActions(incrementBond=['Ctc'],          decrementBond=['Csc'],          formBond=['Cdc'],       breakBond=['Cdc','Cd','CO','CS'],incrementRadical=['Cdc'],decrementRadical=['Cdc'],incrementLonePair=['C2dc'],decrementLonePair=[])
atomTypes['Cdd' ].setActions(incrementBond=[],               decrementBond=['Cd','CO','CS'], formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['Ct'  ].setActions(incrementBond=[],               decrementBond=['Cd','CO','CS'], formBond=['Ct','Ctc'],  breakBond=['Ct'],        incrementRadical=['Ct'],   decrementRadical=['Ct'],   incrementLonePair=['C2tc'],decrementLonePair=[])
atomTypes['Ctc' ].setActions(incrementBond=[],               decrementBond=['Cdc'],          formBond=['Ctc'],       breakBond=['Ctc','Ct'],   incrementRadical=['Ctc'],  decrementRadical=['Ctc'],  incrementLonePair=['C2tc'],decrementLonePair=[])
atomTypes['Cb'  ].setActions(incrementBond=['Cbf'],          decrementBond=[],               formBond=['Cb'],        breakBond=['Cb'],        incrementRadical=['Cb'],   decrementRadical=['Cb'],   incrementLonePair=[],      decrementLonePair=[])
atomTypes['Cbf' ].setActions(incrementBond=[],               decrementBond=['Cb'],           formBond=[],            breakBond=['Cb'],        incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['C2s' ].setActions(incrementBond=['C2d'],          decrementBond=[],               formBond=['C2s'],       breakBond=['C2s'],       incrementRadical=['C2s'],  decrementRadical=['C2s'],  incrementLonePair=['Ca'],  decrementLonePair=['Cs'])
atomTypes['C2sc'].setActions(incrementBond=['C2dc'],         decrementBond=[],               formBond=['C2sc'],      breakBond=['C2sc'],      incrementRadical=['C2sc'], decrementRadical=['C2sc'], incrementLonePair=[],      decrementLonePair=['Cs'])
atomTypes['C2d' ].setActions(incrementBond=['C2tc'],         decrementBond=['C2s'],          formBond=['C2dc'],      breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['Cd','CO','CS'])
atomTypes['C2dc'].setActions(incrementBond=[],               decrementBond=['C2sc'],         formBond=[],            breakBond=['C2d'],       incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['Cdc'])
atomTypes['C2tc'].setActions(incrementBond=[],               decrementBond=['C2d'],          formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['Ct','Ctc'])

atomTypes['N'   ].setActions(incrementBond=['N'],            decrementBond=['N'],            formBond=['N'],         breakBond=['N'],         incrementRadical=['N'],    decrementRadical=['N'],    incrementLonePair=['N'],   decrementLonePair=['N'])
atomTypes['N0sc'].setActions(incrementBond=['N0dc'],         decrementBond=[],               formBond=['N0sc'],      breakBond=['N0sc'],      incrementRadical=['N0sc'], decrementRadical=['N0sc'], incrementLonePair=[],      decrementLonePair=['N1s','N1sc'])
atomTypes['N0dc'].setActions(incrementBond=[],               decrementBond=['N0sc'],         formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['N1dc'])
atomTypes['N1s' ].setActions(incrementBond=['N1dc'],         decrementBond=[],               formBond=['N1s'],       breakBond=['N1s'],       incrementRadical=['N1s'],  decrementRadical=['N1s'],  incrementLonePair=['N0sc'],decrementLonePair=['N3s'])
atomTypes['N1sc'].setActions(incrementBond=[],               decrementBond=[],               formBond=['N1sc'],      breakBond=['N1sc'],      incrementRadical=['N1sc'], decrementRadical=['N1sc'], incrementLonePair=[],      decrementLonePair=['N3s'])
atomTypes['N1dc'].setActions(incrementBond=['N1tc','N1dc'],  decrementBond=['N1s','N1dc'],   formBond=['N1dc'],      breakBond=['N1dc'],      incrementRadical=['N1dc'], decrementRadical=['N1dc'], incrementLonePair=['N0dc'],decrementLonePair=['N3d','N3dc'])
atomTypes['N1tc'].setActions(incrementBond=[],               decrementBond=['N1dc'],         formBond=['N1tc'],      breakBond=['N1tc'],      incrementRadical=['N1tc'], decrementRadical=['N1tc'], incrementLonePair=[],      decrementLonePair=['N3tdc'])
atomTypes['N3s' ].setActions(incrementBond=['N3d'],          decrementBond=[],               formBond=['N3s'],       breakBond=['N3s'],       incrementRadical=['N3s'],  decrementRadical=['N3s'],  incrementLonePair=['N1s','N1sc'],decrementLonePair=['N5sc'])
atomTypes['N3d' ].setActions(incrementBond=['N3t'],          decrementBond=['N3s'],          formBond=['N3d'],       breakBond=['N3d'],       incrementRadical=['N3d'],  decrementRadical=['N3d'],  incrementLonePair=['N1dc'],decrementLonePair=['N5dc'])
atomTypes['N3dc'].setActions(incrementBond=['N3t','N3dc','N3tdc'],decrementBond=['N3s','N3dc'],formBond=['N3dc'],    breakBond=['N3dc'],      incrementRadical=['N3dc'], decrementRadical=['N3dc'], incrementLonePair=['N1dc'],decrementLonePair=['N5dc','N5ddc'])
atomTypes['N3t' ].setActions(incrementBond=[],               decrementBond=['N3d'],          formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['N5t'])
atomTypes['N3tdc'].setActions(incrementBond=['N3tdc'],       decrementBond=['N3dc','N3tdc'], formBond=['N3tdc'],     breakBond=['N3tdc'],     incrementRadical=['N3tdc'],decrementRadical=['N3tdc'],incrementLonePair=['N1tc'],decrementLonePair=['N5tc','N5tdc'])
atomTypes['N3b' ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['N5sc'].setActions(incrementBond=['N5dc'],         decrementBond=[],               formBond=['N5sc'],      breakBond=['N5sc'],      incrementRadical=['N5sc'], decrementRadical=['N5sc'], incrementLonePair=['N3s'], decrementLonePair=[])
atomTypes['N5d' ].setActions(incrementBond=['N5ddc','N5tc','N5dd'],decrementBond=['N5sc'],   formBond=['N5d'],       breakBond=['N5d'],       incrementRadical=['N5d'],  decrementRadical=['N5d'],  incrementLonePair=['N3d'], decrementLonePair=[])
atomTypes['N5dc'].setActions(incrementBond=['N5ddc','N5tc','N5dd'],decrementBond=['N5sc'],   formBond=['N5dc'],      breakBond=['N5dc'],      incrementRadical=['N5dc'], decrementRadical=['N5dc'], incrementLonePair=['N3d','N3dc'],decrementLonePair=[])
atomTypes['N5dd'].setActions(incrementBond=['N5td'],         decrementBond=['N5dc'],         formBond=['N5dd'],      breakBond=['N5dd'],      incrementRadical=['N5dd'], decrementRadical=['N5dd'], incrementLonePair=[],      decrementLonePair=[])
atomTypes['N5ddc'].setActions(incrementBond=['N5tdc','N5dddc'],decrementBond=['N5dc'],       formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=['N3dc'],decrementLonePair=[])
atomTypes['N5dddc'].setActions(incrementBond=['N5tdc'],      decrementBond=['N5ddc'],        formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['N5t' ].setActions(incrementBond=['N5td'],         decrementBond=['N5d'],          formBond=['N5t'],       breakBond=['N5t'],       incrementRadical=['N5t'],  decrementRadical=['N5t'],  incrementLonePair=['N3t'], decrementLonePair=[])
atomTypes['N5tc'].setActions(incrementBond=[],               decrementBond=['N5dc'],         formBond=['N5tc'],      breakBond=['N5tc'],      incrementRadical=['N5tc'], decrementRadical=['N5tc'], incrementLonePair=['N3tdc'],decrementLonePair=[])
atomTypes['N5td'].setActions(incrementBond=[],               decrementBond=['N5dd','N5t'],   formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['N5tdc'].setActions(incrementBond=['N5tdc'],       decrementBond=['N5ddc','N5dddc','N5tdc'],formBond=['N5tdc'],breakBond=['N5tdc'], incrementRadical=['N5tdc'],decrementRadical=['N5tdc'],incrementLonePair=['N3tdc'],decrementLonePair=[])
atomTypes['N5b' ].setActions(incrementBond=[],               decrementBond=[],               formBond=['N5b'],       breakBond=['N5b'],       incrementRadical=['N5b'],  decrementRadical=['N5b'],  incrementLonePair=[],      decrementLonePair=[])

atomTypes['O'   ].setActions(incrementBond=['O'],            decrementBond=['O'],            formBond=['O'],         breakBond=['O'],         incrementRadical=['O'],    decrementRadical=['O'],    incrementLonePair=['O'],   decrementLonePair=['O'])
atomTypes['Oa'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=['O0sc'],      breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['O2s'])
atomTypes['O0sc'].setActions(incrementBond=['O0dc'],         decrementBond=[],               formBond=['O0sc'],      breakBond=['Oa','O0sc'], incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['O2s','O2sn'])
atomTypes['O0dc'].setActions(incrementBond=[],               decrementBond=['O0sc'],         formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['O2d','O2dc'])
atomTypes['O2s' ].setActions(incrementBond=['O2d'],          decrementBond=[],               formBond=['O2s','O2sp','O2sn'],breakBond=['O2s'],incrementRadical=['O2s'],  decrementRadical=['O2s'],  incrementLonePair=['Oa','O0sc'],decrementLonePair=['O4sc'])
atomTypes['O2sp'].setActions(incrementBond=['O2d'],          decrementBond=[],               formBond=[],            breakBond=['O2s'],       incrementRadical=['O2sp'], decrementRadical=['O2sp'], incrementLonePair=[],      decrementLonePair=['O4sc'])
atomTypes['O2sn'].setActions(incrementBond=['O2dc'],         decrementBond=[],               formBond=['O2sn'],      breakBond=['O2s','O2sn'],incrementRadical=['O2sn'], decrementRadical=['O2sn'], incrementLonePair=['O0sc'],decrementLonePair=['O4sc'])
atomTypes['O2d' ].setActions(incrementBond=['O4tc'],         decrementBond=['O2s','O2sp','O2d'],formBond=['O2dc'],   breakBond=[],            incrementRadical=['O2d'],  decrementRadical=['O2d'],  incrementLonePair=['O0dc'],decrementLonePair=['O4dc','O4tc'])
atomTypes['O2dc'].setActions(incrementBond=['O2tc'],         decrementBond=['O2sn'],         formBond=['O2dc'],      breakBond=['O2d','O2dc'],incrementRadical=['O2dc'], decrementRadical=['O2dc'], incrementLonePair=['O0dc'],decrementLonePair=['O4dc','O4tc'])
atomTypes['O2tc'].setActions(incrementBond=[],               decrementBond=['O2dc'],         formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['O4tc'])
atomTypes['O4sc'].setActions(incrementBond=['O4dc'],         decrementBond=[],               formBond=['O4sc'],      breakBond=['O4sc'],      incrementRadical=['O4sc'], decrementRadical=['O4sc'], incrementLonePair=['O2s','O2sp'],decrementLonePair=[])
atomTypes['O4dc'].setActions(incrementBond=['O4dc'],         decrementBond=['O4dc','O4sc'],  formBond=['O4dc'],      breakBond=['O4dc'],      incrementRadical=['O4dc'], decrementRadical=['O4dc'], incrementLonePair=['O2d','O2dc'],decrementLonePair=[])
atomTypes['O4tc'].setActions(incrementBond=[],               decrementBond=['O2d'],          formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=['O2d','O2dc','O2tc'],decrementLonePair=[])

atomTypes['Ne'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=['Ne'],   decrementRadical=['Ne'],   incrementLonePair=[],      decrementLonePair=[])

atomTypes['Si'  ].setActions(incrementBond=['Si'],           decrementBond=['Si'],           formBond=['Si'],        breakBond=['Si'],        incrementRadical=['Si'],   decrementRadical=['Si'],   incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sis' ].setActions(incrementBond=['Sid','SiO'],    decrementBond=[],               formBond=['Sis'],       breakBond=['Sis'],       incrementRadical=['Sis'],  decrementRadical=['Sis'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sid' ].setActions(incrementBond=['Sidd','Sit'],   decrementBond=['Sis'],          formBond=['Sid'],       breakBond=['Sid'],       incrementRadical=['Sid'],  decrementRadical=['Sid'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sidd'].setActions(incrementBond=[],               decrementBond=['Sid','SiO'],    formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sit' ].setActions(incrementBond=[],               decrementBond=['Sid'],          formBond=['Sit'],       breakBond=['Sit'],       incrementRadical=['Sit'],  decrementRadical=['Sit'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['SiO' ].setActions(incrementBond=['Sidd'],         decrementBond=['Sis'],          formBond=['SiO'],       breakBond=['SiO'],       incrementRadical=['SiO'],  decrementRadical=['SiO'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sib' ].setActions(incrementBond=[],               decrementBond=[],               formBond=['Sib'],       breakBond=['Sib'],       incrementRadical=['Sib'],  decrementRadical=['Sib'],  incrementLonePair=[],      decrementLonePair=[])
atomTypes['Sibf'].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])

atomTypes['S'   ].setActions(incrementBond=['S'],            decrementBond=['S'],            formBond=['S'],         breakBond=['S'],         incrementRadical=['S'],    decrementRadical=['S'],    incrementLonePair=['S'],   decrementLonePair=['S'])
atomTypes['S0c'].setActions(incrementBond=['S0c'],           decrementBond=['S0c'],          formBond=['S0c'],       breakBond=['Sa','S0c'],  incrementRadical=['S0c'],  decrementRadical=['S0c'],  incrementLonePair=[],      decrementLonePair=['S2s','S2sc','S2dc','S2tc'])
atomTypes['Sa'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=['S0c'],       breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['S2s'])
atomTypes['S2s' ].setActions(incrementBond=['S2d','S2dc'],   decrementBond=[],               formBond=['S2s','S2sc'],breakBond=['S2s'],       incrementRadical=['S2s'],  decrementRadical=['S2s'],  incrementLonePair=['Sa','S0c'],decrementLonePair=['S4s','S4sc'])
atomTypes['S2sc'].setActions(incrementBond=['S2dc'],         decrementBond=[],               formBond=['S2sc'],      breakBond=['S2sc','S2s'],incrementRadical=['S2sc'], decrementRadical=['S2sc'], incrementLonePair=['S0c'], decrementLonePair=['S4s','S4sc'])
atomTypes['S2d' ].setActions(incrementBond=['S2tc'],         decrementBond=['S2s'],          formBond=['S2d'],       breakBond=['S2d'],       incrementRadical=['S2d'],  decrementRadical=['S2d'],  incrementLonePair=[],      decrementLonePair=['S4dc','S4d'])
atomTypes['S2dc'].setActions(incrementBond=['S2tc','S2dc'],  decrementBond=['S2sc','S2s','S2dc'],formBond=['S2dc'],  breakBond=['S2dc'],      incrementRadical=['S2dc'], decrementRadical=['S2dc'], incrementLonePair=['S0c'], decrementLonePair=['S4d','S4dc'])
atomTypes['S2tc'].setActions(incrementBond=[],               decrementBond=['S2d','S2dc'],   formBond=['S2tc'],      breakBond=['S2tc'],      incrementRadical=['S2tc'], decrementRadical=['S2tc'], incrementLonePair=['S0c'], decrementLonePair=['S4t'])
atomTypes['S4s' ].setActions(incrementBond=['S4d','S4dc'],   decrementBond=[],               formBond=['S4s'],       breakBond=['S4s'],       incrementRadical=['S4s'],  decrementRadical=['S4s'],  incrementLonePair=['S2s','S2sc'],decrementLonePair=['S6s'])
atomTypes['S4sc'].setActions(incrementBond=['S4d','S4dc'],   decrementBond=[],               formBond=['S4s','S4sc'],breakBond=['S4sc'],      incrementRadical=['S4sc'], decrementRadical=['S4sc'], incrementLonePair=['S2s','S2sc'],decrementLonePair=['S6s'])
atomTypes['S4d' ].setActions(incrementBond=['S4dd','S4dc','S4t','S4tdc'],decrementBond=['S4s','S4sc'],formBond=['S4dc','S4d'],breakBond=['S4d','S4dc'],incrementRadical=['S4d'],decrementRadical=['S4d'],incrementLonePair=['S2d','S2dc'],decrementLonePair=['S6d','S6dc'])
atomTypes['S4dc'].setActions(incrementBond=['S4dd','S4dc','S4tdc'],decrementBond=['S4sc','S4dc'],formBond=['S4d','S4dc'],breakBond=['S4d','S4dc'],incrementRadical=['S4dc'],decrementRadical=['S4dc'],incrementLonePair=['S2d','S2dc'],decrementLonePair=['S6d','S6dc'])
atomTypes['S4b' ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['S4dd'].setActions(incrementBond=['S4dc'],         decrementBond=['S4dc','S4d'],   formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=['S6dd'])
atomTypes['S4t' ].setActions(incrementBond=[],               decrementBond=['S4d'],          formBond=['S4t'],       breakBond=['S4t'],       incrementRadical=['S4t'],  decrementRadical=['S4t'],  incrementLonePair=['S2tc'],decrementLonePair=['S6t','S6tdc'])
atomTypes['S4tdc'].setActions(incrementBond=['S4tdc'],       decrementBond=['S4d','S4tdc'],  formBond=['S4tdc'],     breakBond=['S4tdc'],     incrementRadical=['S4tdc'],decrementRadical=['S4tdc'],incrementLonePair=['S6tdc'],decrementLonePair=['S6td','S6tdc'])
atomTypes['S6s' ].setActions(incrementBond=['S6d','S6dc'],   decrementBond=[],               formBond=['S6s'],       breakBond=['S6s'],       incrementRadical=['S6s'],  decrementRadical=['S6s'],  incrementLonePair=['S4s','S4sc'],decrementLonePair=[])
atomTypes['S6sc'].setActions(incrementBond=['S6dc'],         decrementBond=[],               formBond=['S6sc'],      breakBond=['S6sc'],      incrementRadical=['S6sc'], decrementRadical=['S6sc'], incrementLonePair=['S4s','S4sc'],decrementLonePair=[])
atomTypes['S6d' ].setActions(incrementBond=['S6dd','S6t','S6tdc'],decrementBond=['S6s'],     formBond=['S6d','S6dc'],breakBond=['S6d','S6dc'],incrementRadical=['S6d'],  decrementRadical=['S6d'],  incrementLonePair=['S4d','S4dc'], decrementLonePair=[])
atomTypes['S6dc'].setActions(incrementBond=['S6dd','S6ddd','S6dc','S6t','S6td','S6tdc'],decrementBond=['S6sc','S6dc'],formBond=['S6d','S6dc'],breakBond=['S6d','S6dc'],incrementRadical=['S6dc'],decrementRadical=['S6dc'],incrementLonePair=['S4d','S4dc'],decrementLonePair=[])
atomTypes['S6dd'].setActions(incrementBond=['S6ddd','S6td'], decrementBond=['S6d','S6dc'],   formBond=['S6dd','S6dc'],breakBond=['S6dd'],     incrementRadical=['S6dd'], decrementRadical=['S6dd'], incrementLonePair=['S4dd'],decrementLonePair=[])
atomTypes['S6ddd'].setActions(incrementBond=[],              decrementBond=['S6dd','S6dc'],  formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['S6t' ].setActions(incrementBond=['S6td'],         decrementBond=['S6d','S6dc'],   formBond=['S6t'],       breakBond=['S6t'],       incrementRadical=['S6t'],  decrementRadical=['S6t'],  incrementLonePair=['S4t'], decrementLonePair=[])
atomTypes['S6td'].setActions(incrementBond=['S6tt','S6tdc'], decrementBond=['S6dc','S6t','S6dd','S6tdc'],formBond=['S6td'],breakBond=['S6td'],incrementRadical=['S6td'], decrementRadical=['S6td'], incrementLonePair=['S4tdc'],decrementLonePair=[])
atomTypes['S6tt'].setActions(incrementBond=[],               decrementBond=['S6td','S6tdc'], formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])
atomTypes['S6tdc'].setActions(incrementBond=['S6td','S6tdc','S6tt'],decrementBond=['S6dc','S6tdc'],formBond=['S6tdc'],breakBond=['S6tdc'],    incrementRadical=['S6tdc'],decrementRadical=['S6tdc'],incrementLonePair=['S4t','S4tdc'],decrementLonePair=[])


atomTypes['Cl'  ].setActions(incrementBond=[],               decrementBond=['Cl'],           formBond=['Cl'],        breakBond=['Cl'],        incrementRadical=['Cl'],   decrementRadical=['Cl'],   incrementLonePair=[],      decrementLonePair=[])

atomTypes['Ar'  ].setActions(incrementBond=[],               decrementBond=[],               formBond=[],            breakBond=[],            incrementRadical=[],       decrementRadical=[],       incrementLonePair=[],      decrementLonePair=[])

#list of elements that do not have more specific atomTypes
#these are ordered on priority of picking if we encounter a more general atomType for make
allElements=['H', 'C', 'O', 'N', 'S', 'Si', 'Cl', 'Ne', 'Ar', 'He',]
nonSpecifics=['H', 'He', 'Ne', 'Cl', 'Ar',]

for atomType in atomTypes.values():
    for items in [atomType.generic, atomType.specific,
      atomType.incrementBond, atomType.decrementBond, atomType.formBond,
      atomType.breakBond, atomType.incrementRadical, atomType.decrementRadical, atomType.incrementLonePair, atomType.decrementLonePair]:
        for index in range(len(items)):
            items[index] = atomTypes[items[index]]


def getFeatures(atom, bonds):
    """
    Returns a list of features needed to determine atomType for :class:'Atom'
    or :class:'GroupAtom' object 'atom and with local bond structure `bonds`,
    a ``dict`` containing atom-bond pairs.

    """
    cython.declare(single=cython.int, allDouble=cython.int, rDouble=cython.int,
                   sDouble=cython.int, oDouble=cython.int, triple=cython.int,
                   benzene=cython.int)
    cython.declare(features=cython.list)

    # Count numbers of each higher-order bond type
    single = 0; rDouble = 0; oDouble = 0; sDouble = 0; triple = 0; benzene = 0
    for atom2, bond12 in bonds.iteritems():
        if bond12.isSingle():
            single += 1
        elif bond12.isDouble():
            if atom2.isOxygen():
                oDouble += 1
            elif atom2.isSulfur():
                sDouble += 1
            else:
                # rDouble is for double bonds NOT to oxygen or Sulfur
                rDouble += 1
        elif bond12.isTriple(): triple += 1
        elif bond12.isBenzene(): benzene += 1

    # allDouble is for all double bonds, to anything
    allDouble = rDouble + oDouble + sDouble
    features = [single, allDouble, rDouble, oDouble, sDouble, triple, benzene, atom.lonePairs, atom.charge]

    return features

def getAtomType(atom, bonds):
    """
    Determine the appropriate atom type for an :class:`Atom` object `atom`
    with local bond structure `bonds`, a ``dict`` containing atom-bond pairs.
    """

    cython.declare(atomSymbol=str)
    cython.declare(molFeatureList=cython.list, atomTypeFeatureList=cython.list)

    # Use element and counts to determine proper atom type
    atomSymbol = atom.symbol
    #These elements do not do not have a more specific atomType
    if atomSymbol in nonSpecifics:
        return atomTypes[atomSymbol]

    molFeatureList = getFeatures(atom, bonds)
    for specificAtomType in atomTypes[atomSymbol].specific:
        atomtypeFeatureList = specificAtomType.getFeatures()
        for molFeature, atomtypeFeature in zip(molFeatureList, atomtypeFeatureList):
            if atomtypeFeature == []:
                continue
            elif molFeature not in atomtypeFeature:
                break
        else:
            return specificAtomType
    else:
        single = molFeatureList[0]
        rDouble = molFeatureList[2]
        oDouble = molFeatureList[3]
        sDouble = molFeatureList[4]
        triple = molFeatureList[5]
        benzene = molFeatureList[6]
        lonePairs = molFeatureList[7]
        charge = molFeatureList[8]

        raise AtomTypeError(
            'Unable to determine atom type for atom {0}, which has {1:d} single bonds, {2:d} double bonds to C, {3:d} double bonds to O, {4:d} double bonds to S, {5:d} triple bonds, {6:d} benzene bonds, {7:d} lone pairs, and {8:d} charge.'.format(
                atom, single, rDouble, oDouble, sDouble, triple, benzene, lonePairs, charge))

