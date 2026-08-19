#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
This module provides methods for translating to and from common molecule
representation formats, e.g. SMILES, InChI, SMARTS.
"""

import logging
import re

import cython
# Assume that rdkit is installed
from rdkit import Chem
# Test if openbabel is installed
try:
    from openbabel import openbabel
except ImportError:
    BACKENDS = ['rdkit']
else:
    BACKENDS = ['openbabel', 'rdkit']

import rmgpy.molecule.inchi as inchiutil
import rmgpy.molecule.molecule as mm
import rmgpy.molecule.util as util
from rmgpy.exceptions import DependencyError
from rmgpy.molecule.converter import to_rdkit_mol, from_rdkit_mol, to_ob_mol, from_ob_mol

# constants

INCHI_LOOKUPS = {
    'H': '[H]',  # RDkit was improperly handling the Hydrogen radical from InChI
    'He': '[He]',
}
SMILES_LOOKUPS = {
    '[He]':  # RDKit improperly handles helium and returns it in a triplet state
        """
        He
        multiplicity 1
        1 He u0 p1
        """,
    '[Ar]':  # RDKit improperly handles argon
        """
        Ar
        multiplicity 1
        1 Ar u0 p4
        """,
    '[C]':  # We'd return the quintuplet without this
        """
        multiplicity 3
        1 C u2 p1 c0
        """,
    '[CH]':  # We'd return the quartet without this
        """
        multiplicity 2
        1 C u1 p1 c0 {2,S}
        2 H u0 p0 c0 {1,S}
        """,
    '[C]F':  # We'd return the quartet without this
        """
        multiplicity 2
        1 C u1 p1 c0 {2,S}
        2 F u0 p3 c0 {1,S}
        """,
    '[C]Cl':  # We'd return the quartet without this
        """
        multiplicity 2
        1 C u1 p1 c0 {2,S}
        2 Cl u0 p3 c0 {1,S}
        """,
    '[C]Br':  # We'd return the quartet without this
        """
        multiplicity 2
        1 C u1 p1 c0 {2,S}
        2 Br u0 p3 c0 {1,S}
        """,
    '[X]':  # Surface site
        """
        multiplicity 1
        1 X u0
        """,
    'e':
        """
        multiplicity 1
        1 e u0 p0 c-1
        """,
    '[H+]':
        """
        multiplicity 1
        1 H u0 p0 c+1
        """,
}

#: This dictionary is used to shortcut lookups of a molecule's SMILES string from its chemical formula.
#: The key is the charge signature returned by :func:`get_charge_signature` - the formula together
#: with the charges its atoms carry - because a formula alone does not identify a species.
MOLECULE_LOOKUPS = {
    ('N2', ()): 'N#N',
    ('CH4', ()): 'C',
    ('H2O', ()): 'O',
    ('C2H6', ()): 'CC',
    ('H2', ()): '[H][H]',
    ('H2O2', ()): 'OO',
    ('C3H8', ()): 'CCC',
    ('Ar', ()): '[Ar]',
    ('He', ()): '[He]',
    ('CH4O', ()): 'CO',
    ('CO', (-1, 1)): '[C-]#[O+]',
    ('O2', ()): 'O=O',
    ('C', ()): '[C]',  # for this to be in the "molecule" list it must be singlet with 2 lone pairs
    ('H2S', ()): 'S',
    ('NH3', ()): 'N',
    ('O3', (-1, 1)): '[O-][O+]=O',
    ('Cl2', ()): '[Cl][Cl]',
    ('ClH', ()): 'Cl',
    ('I2', ()): '[I][I]',
    ('HI', ()): 'I',
    ('H', (1,)): '[H+]',  # the proton; bare 'H+' is not parseable SMILES
    ('e', (-1,)): 'e',
}

RADICAL_LOOKUPS = {
    ('CH3', ()): '[CH3]',
    ('HO', ()): '[OH]',
    ('C2H5', ()): 'C[CH2]',
    ('O', ()): '[O]',
    ('S', ()): '[S]',
    ('N', ()): '[N]',
    ('HO2', ()): '[O]O',
    ('CH', ()): '[CH]',
    ('CH2', ()): '[CH2]',
    ('H', ()): '[H]',
    ('C', ()): '[C]',  # this, in the radical list, could be triplet or quintet.
    ('O2', ()): '[O][O]',
    ('S2', ()): '[S][S]',
    ('OS', ()): '[S][O]',
    ('HS', ()): '[SH]',
    ('H2N', ()): '[NH2]',
    ('HN', ()): '[NH]',
    ('NO', ()): '[N]=O',
    ('F', ()): '[F]',
    ('Cl', ()): '[Cl]',
    ('Br', ()): '[Br]',
    ('I', ()): '[I]',
    ('CF', ()): '[C]F',
    ('CCl', ()): '[C]Cl',
    ('CBr', ()): '[C]Br',
    ('e', (-1,)): 'e',
}

#: Matches a bracket atom in a SMILES string; SMILES cannot write a charge anywhere else.
_BRACKET_ATOM = re.compile(r'\[([^\]]*)\]')
#: Matches the charge at the end of a bracket atom's contents: ``-``, ``++``, ``+2``, ...
_BRACKET_CHARGE = re.compile(r'([+-])(\1*)(\d*)$')


def to_inchi(mol, backend='rdkit-first', aug_level=0):
    """
    Convert a molecular structure to an InChI string.
    For aug_level=0, generates the canonical InChI.
    For aug_level=1, appends the molecule multiplicity.
    For aug_level=2, appends positions of unpaired and paired electrons.

    Uses RDKit or OpenBabel for conversion.

    Args:
        backend     choice of backend, 'rdkit-first' (default), 'openbabel-first', 'rdkit', or 'openbabel'
        aug_level   level of augmentation, 0, 1, or 2
    """
    cython.declare(inchi=str, ulayer=str, player=str, mlayer=str)

    if aug_level == 0:
        return _write(mol, 'inchi', backend)

    elif aug_level == 1:
        inchi = to_inchi(mol, backend=backend)

        mlayer = '/mult{0}'.format(mol.multiplicity) if mol.multiplicity != 0 else ''

        return inchi + mlayer

    elif aug_level == 2:
        inchi = to_inchi(mol, backend=backend)

        ulayer, player = inchiutil.create_augmented_layers(mol)

        return inchiutil.compose_aug_inchi(inchi, ulayer, player)

    else:
        raise ValueError("Implemented values for aug_level are 0, 1, or 2.")


def to_inchi_key(mol, backend='rdkit-first', aug_level=0):
    """
    Convert a molecular structure to an InChI Key string.
    For aug_level=0, generates the canonical InChI.
    For aug_level=1, appends the molecule multiplicity.
    For aug_level=2, appends positions of unpaired and paired electrons.

    Uses RDKit or OpenBabel for conversion.

    Args:
        backend     choice of backend, 'rdkit-first' (default), 'openbabel-first', 'rdkit', or 'openbabel'
        aug_level   level of augmentation, 0, 1, or 2
    """
    cython.declare(key=str, ulayer=str, player=str, mlayer=str)

    if aug_level == 0:
        return _write(mol, 'inchikey', backend)

    elif aug_level == 1:
        key = to_inchi_key(mol, backend=backend)

        mlayer = '-mult{0}'.format(mol.multiplicity) if mol.multiplicity != 0 else ''

        return key + mlayer

    elif aug_level == 2:
        key = to_inchi_key(mol, backend=backend)

        ulayer, player = inchiutil.create_augmented_layers(mol)

        return inchiutil.compose_aug_inchi_key(key, ulayer, player)

    else:
        raise ValueError("Implemented values for aug_level are 0, 1, or 2.")


def to_smarts(mol, backend='rdkit'):
    """
    Convert a molecular structure to an SMARTS string. Uses
    `RDKit <https://rdkit.org/>`_ to perform the conversion.
    Perceives aromaticity and removes Hydrogen atoms.
    """
    return _write(mol, 'sma', backend)


def to_smiles(mol, backend='default'):
    """
    Convert a molecular structure to an SMILES string.

    If there is a Nitrogen/Sulfur atom present it uses
    `OpenBabel <https://openbabel.org/>`_ to perform the conversion,
    and the SMILES may or may not be canonical.

    Otherwise, it uses `RDKit <https://rdkit.org/>`_ to perform the
    conversion, so it will be canonical SMILES.
    While converting to an RDMolecule it will perceive aromaticity
    and removes Hydrogen atoms.

    Raises a `ValueError` if the molecule carries a charge that the resulting
    SMILES does not, rather than returning a string for some other species.
    """
    # If we're going to have to check the formula anyway,
    # we may as well shortcut a few small known molecules.
    # Dictionary lookups are O(1) so this should be fast.
    # The dictionary is defined at the top of this file.
    key = get_charge_signature(mol)
    try:
        if mol.is_radical():
            output = RADICAL_LOOKUPS[key]
        else:
            output = MOLECULE_LOOKUPS[key]
    except KeyError:
        if backend == 'default':
            for atom in mol.atoms:
                if atom.is_nitrogen() or atom.is_sulfur():
                    return _write(mol, 'smi', backend='openbabel')
            return _write(mol, 'smi', backend='rdkit')
        else:
            return _write(mol, 'smi', backend=backend)
    else:
        return output


def get_charge_signature(mol):
    """
    Return the key under which `mol` is looked up in the SMILES shortcut tables:
    its formula, together with the charges its atoms carry.

    A formula on its own does not identify a species. It says nothing about the
    charge, so keying on it alone hands an ion the neutral species' SMILES - which
    is how O(-.) used to be written ``[O]`` and H(-) written ``H+``. The net charge
    is not enough either: an ion pair such as H(+).OH(-) is net neutral and would
    still collide with water's entry. The charges themselves are what distinguish
    them, so they are what the key carries.
    """
    cython.declare(charges=list, atom=mm.Atom)

    charges = []
    for atom in mol.atoms:
        if atom.charge != 0:
            charges.append(atom.charge)
    charges.sort()

    return mol.get_formula(), tuple(charges)


def get_smiles_net_charge(smiles):
    """
    Return the sum of the formal charges written in the SMILES string `smiles`.

    A SMILES string can only carry a charge inside a bracket atom, so the bracket
    contents are the whole of it.
    """
    cython.declare(total=cython.int, body=str, match=object)

    total = 0
    for body in _BRACKET_ATOM.findall(smiles):
        # An atom class (the ':1' of '[CH3:1]') trails the charge and is not one.
        body = body.split(':', 1)[0]
        match = _BRACKET_CHARGE.search(body)
        if match is None:
            continue
        # A charge is written either as repeated signs ('--') or as a count ('-2').
        total += (1 if match.group(1) == '+' else -1) * (int(match.group(3)) if match.group(3)
                                                         else 1 + len(match.group(2)))
    return total


def from_inchi(mol, inchistr, backend='openbabel-first', raise_atomtype_exception=True):
    """
    Convert an InChI string `inchistr` to a molecular structure. Uses
    a user-specified backend for conversion, currently supporting 'openbabel-first' (default), rdkit-first,
    rdkit, and openbabel.
    """
    if inchiutil.INCHI_PREFIX in inchistr:
        return _read(mol, inchistr, 'inchi', backend, raise_atomtype_exception=raise_atomtype_exception)
    else:
        return _read(mol, inchiutil.INCHI_PREFIX + '/' + inchistr, 'inchi', backend,
                     raise_atomtype_exception=raise_atomtype_exception)


def from_augmented_inchi(mol, aug_inchi, raise_atomtype_exception=True):
    """
    Creates a Molecule object from the augmented inchi.

    First, the inchi is converted into a Molecule using
    the backend parsers.

    Next, the multiplicity and unpaired electron information
    is used to fix a number of parsing errors made by the backends.

    Finally, the atom types of the corrected molecule are perceived.

    Returns a Molecule object
    """

    if not isinstance(aug_inchi, inchiutil.AugmentedInChI):
        aug_inchi = inchiutil.AugmentedInChI(aug_inchi)

    mol = from_inchi(mol, aug_inchi.inchi)

    mol.multiplicity = len(aug_inchi.u_indices) + 1 if aug_inchi.u_indices else 1

    inchiutil.fix_molecule(mol, aug_inchi)

    mol.update_atomtypes(log_species=True, raise_exception=raise_atomtype_exception)

    return mol


def from_smarts(mol, smartsstr, backend='rdkit', raise_atomtype_exception=True):
    """
    Convert a SMARTS string `smartsstr` to a molecular structure. Uses
    `RDKit <https://rdkit.org/>`_ to perform the conversion.
    This Kekulizes everything, removing all aromatic atom types.
    """
    return _read(mol, smartsstr, 'sma', backend, raise_atomtype_exception=raise_atomtype_exception)


def from_smiles(mol, smilesstr, backend='openbabel-first', raise_atomtype_exception=True):
    """
    Convert a SMILES string `smilesstr` to a molecular structure. Uses
    a user-specified backend for conversion, currently supporting openbabel-first (default), rdkit-first,
    rdkit and openbabel.
    """
    return _read(mol, smilesstr, 'smi', backend, raise_atomtype_exception=raise_atomtype_exception)


def _rdkit_translator(input_object, identifier_type, mol=None):
    """
    Converts between formats using RDKit. If input is a :class:`Molecule`,
    the identifier_type is used to determine the output type. If the input is
    a `str`, then the identifier_type is used to identify the input, and the
    desired output is assumed to be a :class:`Molecule` object.

    Args:
        input_object: either molecule or string identifier
        identifier_type: format of string identifier
            'inchi'    -> InChI
            'inchikey' -> InChI Key
            'sma'      -> SMARTS
            'smi'      -> SMILES
        mol: molecule object for output (optional)
    """
    if identifier_type == 'inchi' and not Chem.inchi.INCHI_AVAILABLE:
        raise DependencyError("RDKit installed without InChI. Please reinstall to read and write InChI strings.")

    if isinstance(input_object, str):
        # We are converting from a string identifier to a molecule
        if identifier_type == 'inchi':
            rdkitmol = Chem.inchi.MolFromInchi(input_object, removeHs=False)
        elif identifier_type == 'sma':
            rdkitmol = Chem.MolFromSmarts(input_object)
        elif identifier_type == 'smi':
            rdkitmol = Chem.MolFromSmiles(input_object)
        else:
            raise ValueError('Identifier type {0} is not supported for reading using RDKit.'.format(identifier_type))
        if rdkitmol is None:
            raise ValueError("Could not interpret the identifier {0!r}".format(input_object))
        if mol is None:
            mol = mm.Molecule()
        output = from_rdkit_mol(mol, rdkitmol)
    elif isinstance(input_object, mm.Molecule):
        # We are converting from a molecule to a string identifier
        if identifier_type == 'smi':
            rdkitmol = to_rdkit_mol(input_object, sanitize=False)
        else:
            rdkitmol = to_rdkit_mol(input_object, sanitize=True)
        if identifier_type == 'inchi':
            output = Chem.inchi.MolToInchi(rdkitmol, options='-SNon')
        elif identifier_type == 'inchikey':
            inchi = to_inchi(input_object)
            output = Chem.inchi.InchiToInchiKey(inchi)
        elif identifier_type == 'sma':
            output = Chem.MolToSmarts(rdkitmol)
        elif identifier_type == 'smi':
            if input_object.is_aromatic():
                output = Chem.MolToSmiles(rdkitmol)
            else:
                output = Chem.MolToSmiles(rdkitmol, kekuleSmiles=True)
        else:
            raise ValueError('Identifier type {0} is not supported for writing using RDKit.'.format(identifier_type))
    else:
        raise ValueError('Unexpected input format. Should be a Molecule or a string.')

    return output


def _openbabel_translator(input_object, identifier_type, mol=None):
    """
    Converts between formats using OpenBabel. If input is a :class:`Molecule`,
    the identifier_type is used to determine the output type. If the input is
    a `str`, then the identifier_type is used to identify the input, and the
    desired output is assumed to be a :class:`Molecule` object.

    Args:
        input_object: either molecule or string identifier
        identifier_type: format of string identifier
            'inchi'    -> InChI
            'inchikey' -> InChI Key
            'smi'      -> SMILES
        mol: molecule object for output (optional)
    """
    ob_conversion = openbabel.OBConversion()

    if isinstance(input_object, str):
        # We are converting from a string identifier to a Molecule
        ob_conversion.SetInFormat(identifier_type)
        obmol = openbabel.OBMol()
        ob_conversion.ReadString(obmol, input_object)
        obmol.AddHydrogens()
        # In OpenBabel 3+ the function obmol.AssignSpinMultiplicity(True) does nothing.
        # We could write our own method here and call obatom.SetSpinMultiplicity on
        # each atom, but instead we will leave them blank for now and fix them 
        # in the from_ob_mol() method.
        if mol is None:
            mol = mm.Molecule()
        output = from_ob_mol(mol, obmol)
    elif isinstance(input_object, mm.Molecule):
        # We are converting from a Molecule to a string identifier
        if identifier_type == 'inchi':
            ob_conversion.SetOutFormat('inchi')
            ob_conversion.AddOption('w')
        elif identifier_type == 'inchikey':
            ob_conversion.SetOutFormat('inchi')
            ob_conversion.AddOption('w')
            ob_conversion.AddOption('K')
        elif identifier_type == 'smi':
            ob_conversion.SetOutFormat('can')
            # turn off isomer and stereochemistry information
            ob_conversion.AddOption('i')
        else:
            raise ValueError('Unexpected identifier type {0}.'.format(identifier_type))
        obmol = to_ob_mol(input_object)
        output = ob_conversion.WriteString(obmol).strip()
    else:
        raise ValueError('Unexpected input format. Should be a Molecule or a string.')

    return output


def _lookup(mol, identifier, identifier_type):
    """
    Looks up the identifier and parses it the way we think is best.

    For troublesome inchis, we look up the smiles, and parse smiles.
    For troublesome smiles, we look up the adj list, and parse the adj list.

    """
    if identifier_type.lower() == 'inchi':
        try:
            smi = INCHI_LOOKUPS[identifier.split('/', 1)[1]]
            return mol.from_smiles(smi)
        except KeyError:
            return None
    elif identifier_type.lower() == 'smi':
        try:
            adjList = SMILES_LOOKUPS[identifier]
            return mol.from_adjacency_list(adjList)
        except KeyError:
            return None


def _check_smiles_charge(mol, identifier, identifier_type):
    """
    Check that a written SMILES carries the charge the molecule has.

    A wrong-but-parseable SMILES is worse than no SMILES at all: everything keyed on it
    - Chemkin and YAML export, species dictionaries, library matching - then silently
    treats the ion as whatever neutral species the string names. Returning False here
    lets :func:`_write` try the next backend, and if none can carry the charge it raises.

    Only charged molecules are checked, so neutral chemistry pays nothing for this.
    """
    cython.declare(charge=cython.int)

    if identifier_type != 'smi':
        return True

    charge = mol.get_net_charge()
    if charge == 0:
        return True

    if get_smiles_net_charge(identifier) != charge:
        logging.error('The SMILES {0!r} does not carry the net charge {1:+d} of this molecule:\n'
                      '{2}'.format(identifier, charge, mol.to_adjacency_list()))
        return False

    return True


def _check_output(mol, identifier):
    """Check if molecule object has been correctly parsed."""
    conditions = []

    # Check that the molecule has atoms
    conditions.append(bool(mol.atoms))
    # Check that the identifier is not blank
    conditions.append(bool(identifier.strip()))

    # Check that the InChI element count matches the molecule
    if 'InChI=1' in identifier:
        inchi_elementcount = util.get_element_count(identifier)
        mol_elementcount = util.get_element_count(mol)
        conditions.append(inchi_elementcount == mol_elementcount)

    return all(conditions)


def _read(mol, identifier, identifier_type, backend, raise_atomtype_exception=True):
    """
    Parses the identifier based on the type of identifier (inchi/smi/sma)
    and the backend used.

    First, look up the identifier in a dictionary to see if it can be processed
    this way.

    If not in the dictionary, parse it through the specified backed,
    or try all backends.
    """
    # Check for potential mistakes in input arguments
    if 'InChIKey' in identifier:
        raise ValueError('InChIKey is a write-only format and cannot be parsed.')
    elif 'InChI' in identifier and identifier_type != 'inchi':
        raise ValueError('Improper identifier type "{0}". The provided identifier appears to be an InChI.'.format(identifier_type))

    if _lookup(mol, identifier, identifier_type) is not None:
        if _check_output(mol, identifier):
            mol.update(log_species=True, raise_atomtype_exception=raise_atomtype_exception, sort_atoms=False)
            return mol

    for option in _get_backend_list(backend):
        if option == 'rdkit':
            mol = _rdkit_translator(identifier, identifier_type, mol)
        elif option == 'openbabel':
            mol = _openbabel_translator(identifier, identifier_type, mol)
        else:
            raise NotImplementedError("Unrecognized backend {0}".format(option))

        if _check_output(mol, identifier):
            mol.update(log_species=True, raise_atomtype_exception=raise_atomtype_exception, sort_atoms=False)
            return mol
        else:
            logging.debug('Backend {0} is not able to parse identifier {1}'.format(option, identifier))

    raise ValueError("Unable to correctly parse {0} with backend {1}.".format(identifier, backend))


def _write(mol, identifier_type, backend):
    """
    Converts the input molecule to the specified identifier type.

    Uses backends as specified by the `backend` argument.

    Returns a string identifier of the requested type.
    """
    # Check that the molecule is not empty
    if not mol.atoms:
        return ''

    for option in _get_backend_list(backend):
        if option == 'rdkit':
            try:
                output = _rdkit_translator(mol, identifier_type)
            except ValueError:
                continue
        elif option == 'openbabel':
            try:
                output = _openbabel_translator(mol, identifier_type)
            except ValueError:
                continue
        else:
            raise NotImplementedError("Unrecognized backend {0}".format(option))

        if _check_output(mol, output) and _check_smiles_charge(mol, output, identifier_type):
            return output
        else:
            logging.debug('Backend {0} is not able to generate {1} for this molecule:\n'
                          '{2}'.format(option, identifier_type, mol.to_adjacency_list()))

    logging.error('Unable to generate identifier for this molecule:\n{0}'.format(mol.to_adjacency_list()))

    raise ValueError("Unable to generate identifier type {0} with backend {1}.".format(identifier_type, backend))


def _get_backend_list(backend):
    """
    Returns the appropriate list or iterator of backends given the provided keyword.
    """
    if not isinstance(backend, str):
        raise ValueError("The backend argument should be a string. "
                         "Accepted values are 'openbabel-first', 'rdkit-first', 'rdkit', and 'openbabel'")
    backend = backend.strip().lower()
    if backend == 'openbabel-first':
        return BACKENDS
    elif backend == 'rdkit-first':
        return reversed(BACKENDS)
    elif backend in ['rdkit', 'openbabel']:
        return [backend]
    else:
        raise ValueError("Unrecognized value for backend argument. "
                         "Accepted values are 'openbabel-first', 'rdkit-first', 'rdkit', and 'openbabel'")
