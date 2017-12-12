from .graph cimport Vertex, Edge, Graph
from .molecule cimport Atom, Bond, Molecule

cpdef list populate_resonance_algorithms(dict features=?)

cpdef dict analyze_molecule(Molecule mol)

cpdef list generate_resonance_structures(Molecule mol, bint clar_structures=?, bint keep_isomorphic=?, bint filter_structures=?)

cpdef list _generate_resonance_structures(list mol_list, list method_list, bint keep_isomorphic=?, bint copy=?)

cpdef list generate_ally_delocalization_resonance_structures(Molecule mol)

cpdef list generate_lone_pair_radical_resonance_structures(Molecule mol)

cpdef list generate_lone_pair_multiple_bond_resonance_structures(Molecule mol)

cpdef list generate_lone_pair_radical_multiple_bond_resonance_structures(Molecule mol)

cpdef list generate_N5ddc_N5tc_resonance_structures(Molecule mol)

cpdef list generate_isomorphic_resonance_structures(Molecule mol)

cpdef list generate_aromatic_resonance_structures(Molecule mol, dict features=?)

cpdef list generate_kekule_structure(Molecule mol)

cpdef list generate_opposite_kekule_structure(Molecule mol)

cpdef list generate_clar_structures(Molecule mol)

cpdef list _clar_optimization(Molecule mol, list constraints=?, maxNum=?)

cpdef list _clar_transformation(Molecule mol, list ring)
