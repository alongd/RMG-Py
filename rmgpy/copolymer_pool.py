"""One method-of-moments pool for a linear statistical copolymer (I-240,
slice S5 of the copolymer design).

A ``PolymerSpec`` (schema v2) declares residues, their junctions, the
sequence support, the composition and the sequence statistics.  The method
of moments carries it as one pool whose unit is one residue:

- **Mass.**  ``M_i`` is residue i's molar mass and ``Mbar0 = sum x_i M_i``
  the mean residue molar mass, both exact rationals from the declared graphs
  (``melt_sequence``: the kMC melt's own numbers).  A chain of n residues
  weighs ``n Mbar0 + e`` on average, ``e`` the two end groups.
- **Moments.**  ``mu_k = s sum_d c_d d^k`` (mol, d in residues) over the
  declared chain-length ensemble ``c_d``, with ``s`` the factor that makes
  ``mu_1 Mbar0 + mu_0 e`` the declared initial mass.
- **Period.**  Every present residue must hold the same number of backbone
  atoms; that number is the pool's periods per unit (the kernels' length
  scale), so every backbone cut leaves whole periods.
- **Proxy.**  The pool's RMG species (its proxy oligomers) is stitched from
  one unit: the shortest closed walk of the sequence support through the
  most abundant residue (one residue when it may follow itself; E-P for an
  alternating pair), entered at its first residue's entry site (``*2``) and
  left at its last residue's exit site (``*1``).  The proxy names the pool
  (with its unit mass, polymer-mom's ``unit_molar_mass``, in its
  fingerprint); it never enters the mass or the moments.  A larger unit (every junction once) made a
  three-residue declaration's proxy trimer 27 residues and its construction
  minutes long.

What one pool cannot represent is refused or reported by name: a network
(``refused:network_topology``; its Mw diverges at the gel point), a
branched copolymer (``refused:branched_copolymer``; it needs branch-point
moments), residues of different backbone lengths
(``refused:unequal_residue_backbones``).  Composition drift under
residue-selective chemistry is the kernels' limitation
(``rmgpy.kmc.copolymer_kernels``)."""

from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, Tuple

from rmgpy.molecule import Atom, Bond, Molecule
from rmgpy.molecule.element import get_element

REFUSED_NETWORK = "refused:network_topology"
REFUSED_BRANCHED = "refused:branched_copolymer"
REFUSED_BACKBONES = "refused:unequal_residue_backbones"


class CopolymerPoolRefusal(ValueError):
    """A declaration one copolymer pool cannot carry."""


@dataclass(frozen=True)
class CopolymerPoolSeed:
    label: str
    composition: Tuple[Tuple[str, Fraction], ...]          # present residues, x > 0
    residue_molar_mass: Tuple[Tuple[str, Fraction], ...]   # g/mol, exact
    mean_residue_molar_mass: Fraction                      # Mbar0, g/mol
    end_molar_mass: Fraction                               # head + tail end groups, g/mol
    backbone_atoms_per_residue: int
    dp_counts: Tuple[Tuple[int, int], ...]                 # residues -> representative chains
    moments_mol: Tuple[float, float, float]                # residues
    proxy_word: Tuple[str, ...]


def copolymer_pool_seed(spec) -> CopolymerPoolSeed:
    """The one-pool seed of a linear statistical copolymer (module doc)."""
    from rmgpy.kmc import melt_sequence as msq
    from rmgpy.polymer_spec import PolymerSpec
    if not isinstance(spec, PolymerSpec):
        raise CopolymerPoolRefusal(f"a copolymer pool needs a PolymerSpec, got {type(spec).__name__}")
    if spec.architecture == "network":
        raise CopolymerPoolRefusal(f"{REFUSED_NETWORK}: a network has no moment seed (Mw diverges at the gel "
                                   "point); one pool carries linear chains only")
    if spec.architecture != "linear":
        raise CopolymerPoolRefusal(f"{REFUSED_BRANCHED}: a {spec.architecture} copolymer needs branch-point "
                                   "moments; one pool carries linear chains only")
    composition = tuple((name, x) for name, x in spec.initial_ensemble.composition if x > 0)
    backbone = {name: len(msq.residue_backbone_keys(spec, name)) for name, _x in composition}
    if len(set(backbone.values())) != 1:
        raise CopolymerPoolRefusal(f"{REFUSED_BACKBONES}: residues {backbone} backbone atoms; one pool counts "
                                   "periods only when every residue holds the same number")
    masses = tuple((name, msq.residue_molar_mass(spec, name)) for name, _x in composition)
    mbar = msq.mean_residue_molar_mass(spec)
    ends = msq._graph_mass(spec.head_end_group) + msq._graph_mass(spec.tail_end_group)
    counts = tuple(sorted(spec.initial_ensemble.chains.counts.items()))
    mass = sum((c * (d * mbar + ends) for d, c in counts), Fraction(0))
    scale = Fraction(spec.initial_mass_kg) * 1000 / mass
    moments = tuple(float(scale * sum((c * Fraction(d) ** k for d, c in counts), Fraction(0))) for k in range(3))
    return CopolymerPoolSeed(spec.label, composition, masses, mbar, ends, backbone[composition[0][0]], counts,
                             moments, _proxy_word(spec, [name for name, _x in composition]))


def mean_residue_atoms(spec) -> Tuple[Tuple[str, Fraction], ...]:
    """((element, atoms), ...) of the mean residue, sum x_i (residue i's
    element counts, hydrogens included), exact, sorted by element."""
    residues = dict(spec.residues)
    total: Dict[str, Fraction] = {}
    for name, x in spec.initial_ensemble.composition:
        if not x:
            continue
        for atom in residues[name].graph.atoms:
            total[atom.element] = total.get(atom.element, Fraction(0)) + x
            if atom.hydrogens:
                total["H"] = total.get("H", Fraction(0)) + x * atom.hydrogens
    return tuple(sorted(total.items()))


def _proxy_word(spec, present) -> Tuple[str, ...]:
    """The shortest closed walk of the support (over the present residues)
    through the most abundant residue (ties: name order), from it."""
    x = dict(spec.initial_ensemble.composition)
    start = sorted(present, key=lambda name: (-x[name], name))[0]
    support = dict(spec.sequence_support)
    previous, frontier = {start: None}, [start]
    while frontier:
        reached = []
        for a in frontier:
            for b in sorted(support[a]):
                if b == start:
                    word = [a]
                    while previous[word[-1]] is not None:
                        word.append(previous[word[-1]])
                    return tuple(reversed(word))
                if b in present and b not in previous:
                    previous[b] = a
                    reached.append(b)
        frontier = reached
    raise CopolymerPoolRefusal(f"the sequence support has no closed walk through {start!r}")


def proxy_unit(spec) -> Molecule:
    """The legacy stitching view of the proxy unit (module doc): the proxy
    word's residues joined by their junctions, every declared hydrogen
    explicit, ``*2`` (a stitching radical) on the first residue's entry
    atom and ``*1`` on the last residue's exit atom."""
    from rmgpy.kmc import melt_sequence as msq
    word = copolymer_pool_seed(spec).proxy_word
    residues = dict(spec.residues)
    molecule, maps = Molecule(), []
    entry, exit_ = msq.residue_sites(spec, word[0])[0].atom, msq.residue_sites(spec, word[-1])[1].atom
    for i, name in enumerate(word):
        by_key: Dict[str, Atom] = {}
        for declaration in residues[name].graph.atoms:
            label = "*2" if (i == 0 and declaration.key == entry) else (
                "*1" if (i == len(word) - 1 and declaration.key == exit_) else "")
            atom = Atom(element=get_element(declaration.element, declaration.isotope),
                        radical_electrons=declaration.radical_electrons + (1 if label else 0),
                        charge=declaration.charge, lone_pairs=declaration.lone_pairs, label=label)
            molecule.add_atom(atom)
            by_key[declaration.key] = atom
            for _h in range(declaration.hydrogens):
                hydrogen = Atom(element=get_element("H"), lone_pairs=0)
                molecule.add_atom(hydrogen)
                molecule.add_bond(Bond(atom, hydrogen, 1.0))
        for bond in residues[name].graph.bonds:
            molecule.add_bond(Bond(by_key[bond.atoms[0]], by_key[bond.atoms[1]], bond.order))
        if maps:
            left_atom, right_atom, order = spec.junctions[(word[i - 1], name)]
            molecule.add_bond(Bond(maps[-1][left_atom], by_key[right_atom], order))
        maps.append(by_key)
    molecule.update_multiplicity()
    molecule.update_atomtypes()
    return molecule


def legacy_end_groups(spec) -> Tuple[Molecule, Molecule]:
    """The head cap (``*1``) and tail cap (``*2``), as v1's stitch adapter
    gives them (``LinearHomopolymerSpec.legacy_molecules``)."""
    head, tail = spec.head_end_group, spec.tail_end_group
    return (head.to_legacy_molecule({head.ports[0][0]: "*1"}),
            tail.to_legacy_molecule({tail.ports[0][0]: "*2"}))
