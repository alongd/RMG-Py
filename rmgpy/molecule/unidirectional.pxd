from .graph cimport Vertex, Edge, Graph
from .molecule cimport Atom, Bond, Molecule

cpdef list generate_birad_multiple_bond_unidirectional_transitions(Molecule mol)

cpdef list generate_OS_unidirectional_transitions(Molecule mol)
