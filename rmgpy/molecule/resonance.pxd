from .graph cimport Vertex, Edge, Graph
from .molecule cimport Atom, Bond, Molecule

cpdef list populateResonanceAlgorithms(dict features=?)

cpdef dict analyzeMolecule(Molecule mol)

cpdef list generateResonanceStructures(Molecule mol, bint clarStructures=?, bint keepIsomorphic=?)

cpdef list _generateResonanceStructures(list molList, list methodList, bint keepIsomorphic=?, bint copy=?)

cpdef list generateAdjacentResonanceStructures(Molecule mol)
