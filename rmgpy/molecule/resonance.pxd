from .graph cimport Vertex, Edge, Graph
from .molecule cimport Atom, Bond, Molecule

cpdef list populateResonanceAlgorithms(dict features=?)

cpdef dict analyzeMolecule(Molecule mol)

cpdef list generateResonanceStructures(Molecule mol, bint clarStructures=?, bint keepIsomorphic=?, bint filterStructures=?)

cpdef list _generateResonanceStructures(list molList, list methodList, bint keepIsomorphic=?, bint copy=?)

cpdef list filter_resonance_structures(list molList)

cpdef list generateAdjacentResonanceStructures(Molecule mol)

cpdef list generateLonePairRadicalResonanceStructures(Molecule mol)

cpdef list generateLonePairMultipleBondResonanceStructures(Molecule mol)

cpdef list generateIsomorphicResonanceStructures(Molecule mol, bint saturateH=?)

cpdef list generateAromaticResonanceStructures(Molecule mol, dict features=?)

cpdef list generateKekuleStructure(Molecule mol)

cpdef list generateOppositeKekuleStructure(Molecule mol)

cpdef list generateClarStructures(Molecule mol)

cpdef list _clarOptimization(Molecule mol, list constraints=?, maxNum=?)

cpdef list _clarTransformation(Molecule mol, list ring)
