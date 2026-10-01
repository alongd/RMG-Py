"""Vapour-pressure routing for whole polymer-kMC components.

The constants and correlations in this module are the frozen R-010 v13
volatility law.  The input is an RMG :class:`~rmgpy.molecule.Molecule`; no
cached formula or group data are trusted.
"""

from collections import Counter, deque
from dataclasses import dataclass
import hashlib
import math
from typing import Optional, Tuple

import numpy as np

from rmgpy.molecule import Molecule


R = 8.314462618
ATM_TO_BAR = 1.01325
U = 2.0 ** -53
M_OFF = 1500.0

MASS_TEXT = "C 12.01064\nH 1.007971\n"
MASS_SHA256 = "2e81fc038a470c52af7975ff69e298118ac89b02b54a372a9d2e0c3cd36d9b28"
assert hashlib.sha256(MASS_TEXT.encode()).hexdigest() == MASS_SHA256
ATOMIC_MASS = {line.split()[0]: float(line.split()[1]) for line in MASS_TEXT.splitlines()}
M_C36 = 36 * ATOMIC_MASS["C"] + 74 * ATOMIC_MASS["H"]
assert M_C36 == 506.972894

GROUP_MAP_TEXT = """terminal_CH3|1
aromatic_methyl|3
chain_CH2|4
chain_CH|5
chain_C|6
benzylic_CHx|8
cyclic_CH2|9
cyclic_CH|10
cyclic_C|11
aromatic_CH|15
aromatic_substituted_C|16
internal_CeqC|58
aryl_internal_CeqC|59
terminal_CeqC|61
cyclic_CeqC|62
cyclic_conjugated_diene|88
chain_conjugated_diene|89
ortho_pair|127
meta_pair|128
para_pair|129
nu5|130
n4_22|131
n5|132
n6|133
"""
GROUP_MAP_SHA256 = "5060826e86189b9fa5f4082fcb2acb681b6ff3e891bda3216f0cdd69e8cd4b90"
assert hashlib.sha256(GROUP_MAP_TEXT.encode()).hexdigest() == GROUP_MAP_SHA256

CRIT_COEFF_TEXT = """eq 0.988948 0.699003 0.86074 0.0093898 -0.140414
id tc pc
1 0.04186824 0.0008162027
3 -0.00107095 0.0004165987
4 0.04009768 0.0005262280
5 0.03020686 0.0002300874
6 -0.00387783 -0.0002992479
8 0.00944222 0.0002366517
9 0.02128983 0.0003402668
10 0.02635125 0.0003616188
11 -0.01704586 -0.0005129896
15 0.01611540 0.0002106379
16 0.06820448 0.0004182580
58 0.04515305 0.0007158083
59 0.00000000 0.0000000000
61 0.04544056 0.0009641269
62 0.05640592 0.0003473101
88 0.02473020 -0.0010245100
89 0.00000000 0.0000000000
127 0.00128230 0.0000706118
128 0.00670988 -0.0000724556
129 0.00000000 0.0000000000
130 -0.03382014 -0.0008845745
131 -0.01848152 -0.0002254193
132 -0.02360240 -0.0003245959
133 -0.02458015 -0.0005311340
"""
CRIT_COEFF_SHA256 = "1cb234337365900ba32080d7661207ed0fe11ec99fc54e110a763cb421ba13fa"
assert hashlib.sha256(CRIT_COEFF_TEXT.encode()).hexdigest() == CRIT_COEFF_SHA256

SECONDARY_SMARTS_TEXT = """130 [CX4][CX4]([CX4])([CX4])[cX3]
131 [CX4,CX3][CX4;H1]([CX4,CX3])[CX4;H1]([CX4,CX3])[CX4,CX3]
1311 [CX4,CX3][CX4]([!#6])([CX4,CX3])[CX4]([!#6])([CX4,CX3])[CX4,CX3]
1312 [CX4,CX3][CX4]([!#6])([CX4,CX3])[CX4;H1]([CX4,CX3])[CX4,CX3]
1313 [CX4,CX3][CX4;$([CX4;r3][OX2;r3])]([CX4,CX3])[CX4;$([CX4;r3][OX2;r3])]([CX4,CX3])[CX4,CX3]
132 [CX4,CX3][CX4;H1]([CX4,CX3])[CX4]([CX4,CX3])([CX4,CX3])[CX4,CX3]
1321 [CX4,CX3][CX4]([!#6])([CX4,CX3])[CX4]([CX4,CX3])([CX4,CX3])[CX4,CX3]
133 [CX4,CX3][CX4]([CX4,CX3])([CX4,CX3])[CX4]([CX4,CX3])([CX4,CX3])[CX4,CX3]
"""
SECONDARY_SMARTS_SHA256 = "6c102e33d6874d6342127a54759049bf1dcf55e7ed4f7915f92202487a233b54"
assert hashlib.sha256(SECONDARY_SMARTS_TEXT.encode()).hexdigest() == SECONDARY_SMARTS_SHA256


_TC = {}
_PC = {}
for _line in CRIT_COEFF_TEXT.splitlines()[2:]:
    _group_id, _tc, _pc = _line.split()
    _TC[int(_group_id)] = float(_tc)
    _PC[int(_group_id)] = float(_pc)


@dataclass(frozen=True)
class VolatilityOutcome:
    """Immutable result of a volatility-law evaluation."""

    outcome: str
    error: str
    M: Optional[float] = None
    n_ar: Optional[int] = None
    groups: Tuple[Tuple[int, int], ...] = ()
    Tb: Optional[float] = None
    Tc: Optional[float] = None
    Pc: Optional[float] = None
    log10P: Optional[float] = None
    regime: Optional[str] = None
    flags: Tuple[str, ...] = ()


class _VolatilityError(Exception):
    pass


def _reject(code):
    raise _VolatilityError(code)


def _error(code):
    return VolatilityOutcome("error", code)


def _below_tb(T, Tb):
    """Single mutation point for the frozen ``T <= Tb`` branch boundary."""

    return T <= Tb


def _bond_order(bond):
    order = bond.order
    if isinstance(order, (bool, np.bool_)) or not isinstance(
        order, (int, float, np.integer, np.floating)
    ):
        _reject("E_TYPE_BOND_ORDER")
    order = float(order)
    if not math.isfinite(order):
        _reject("E_NONFINITE_BOND_ORDER")
    if order not in {1.0, 1.5, 2.0}:
        _reject("E_RANGE_BOND_ORDER")
    return order


def _connected(atoms):
    reached = {atoms[0]}
    stack = [atoms[0]]
    while stack:
        for neighbour in stack.pop().bonds:
            if neighbour not in reached:
                reached.add(neighbour)
                stack.append(neighbour)
    return len(reached) == len(atoms)


def _aromatic_projection(molecule, carbon_set):
    try:
        rings, _ = molecule.get_aromatic_rings(save_order=True)
    except (ValueError, RuntimeError, TypeError):
        rings = []
    aromatic_rings = [tuple(ring) for ring in rings if len(ring) == 6 and set(ring) <= carbon_set]
    aromatic_atoms = set()
    aromatic_edges = set()
    for ring in aromatic_rings:
        aromatic_atoms.update(ring)
        for index, atom in enumerate(ring):
            aromatic_edges.add(frozenset((atom, ring[(index + 1) % 6])))
    return aromatic_rings, aromatic_atoms, aromatic_edges


def _compile_graph(molecule, r1_class):
    if not isinstance(molecule, Molecule):
        _reject("E_TYPE_RECORD")
    atoms = tuple(molecule.atoms)
    if not atoms:
        _reject("E_TYPE_ATOMS")
    if not _connected(atoms):
        _reject("E_DOMAIN_DISCONNECTED")
    if any(atom.symbol not in {"C", "H"} for atom in atoms):
        _reject("E_DOMAIN_ELEMENT")
    carbons = tuple(atom for atom in atoms if atom.symbol == "C")
    if not carbons:
        _reject("E_DOMAIN_ELEMENT")
    carbon_set = set(carbons)
    hydrogens = tuple(atom for atom in atoms if atom.symbol == "H")
    if any(len(atom.bonds) != 1 for atom in hydrogens):
        _reject("E_DOMAIN_VALENCE")

    aromatic_rings, aromatic_atoms, aromatic_edges = _aromatic_projection(molecule, carbon_set)
    n_ar = len(aromatic_rings)
    adjacency = {atom: [] for atom in carbons}
    edges = []
    for atom in carbons:
        for neighbour, bond in atom.bonds.items():
            if neighbour.symbol == "H":
                if _bond_order(bond) != 1.0:
                    _reject("E_RANGE_BOND_ORDER")
                continue
            if neighbour not in carbon_set:
                _reject("E_DOMAIN_ELEMENT")
            if id(atom) >= id(neighbour):
                continue
            order = 1.5 if frozenset((atom, neighbour)) in aromatic_edges else _bond_order(bond)
            adjacency[atom].append((neighbour, order))
            adjacency[neighbour].append((atom, order))
            edges.append((atom, neighbour, order))

    explicit_h = Counter()
    for atom in hydrogens:
        neighbour, bond = next(iter(atom.bonds.items()))
        if neighbour not in carbon_set or _bond_order(bond) != 1.0:
            _reject("E_DOMAIN_ELEMENT")
        explicit_h[neighbour] += 1

    h_count = {}
    radicals = 0
    all_h_explicit = bool(hydrogens)
    for atom in carbons:
        radical = atom.radical_electrons
        if isinstance(radical, (bool, np.bool_)) or not isinstance(radical, (int, np.integer)):
            _reject("E_TYPE_ATOM_RADICAL")
        radical = int(radical)
        if radical < 0 or radical > 4:
            _reject("E_RANGE_ATOM_RADICAL")
        heavy_valence = math.fsum(order for _, order in adjacency[atom])
        implicit = 4.0 - radical - heavy_valence
        rounded = round(implicit)
        if abs(implicit - rounded) > 8.0 * U or not 0 <= rounded <= 4:
            _reject("E_DOMAIN_VALENCE")
        h_count[atom] = int(rounded)
        if all_h_explicit and explicit_h[atom] != rounded:
            _reject("E_DOMAIN_VALENCE")
        if atom.charge != 0 or atom.lone_pairs not in {0, -100}:
            _reject("E_DOMAIN_VALENCE")
        radicals += radical

    mass = math.fsum(
        (len(carbons) * ATOMIC_MASS["C"], sum(h_count.values()) * ATOMIC_MASS["H"])
    )

    if radicals:
        # Routing precedes correlation-domain rejection, but the retained
        # record still carries the source groups of every closed-shell site.
        retained_groups = Counter()
        retained_alkenes = {atom for a, b, order in edges if order == 2.0 for atom in (a, b)}
        for a, b, order in edges:
            if order != 2.0:
                continue
            aryl = sum(
                neighbour in aromatic_atoms
                for atom, other in ((a, b), (b, a))
                for neighbour, _ in adjacency[atom]
                if neighbour is not other
            )
            if len(adjacency[a]) == 1 or len(adjacency[b]) == 1:
                retained_groups[61] += 1
            elif aryl:
                retained_groups[59] += aryl
            else:
                retained_groups[58] += 1
        for atom in carbons:
            if atom in retained_alkenes or atom.radical_electrons:
                continue
            h = h_count[atom]
            aromatic_attachments = sum(neighbour in aromatic_atoms for neighbour, _ in adjacency[atom])
            if atom in aromatic_atoms:
                retained_groups[15 if h == 1 else 16] += 1
            elif h == 3 and aromatic_attachments:
                retained_groups[3] += 1
            elif aromatic_attachments:
                retained_groups[8] += aromatic_attachments
            elif h in {3, 4}:
                retained_groups[1] += 1
            elif h == 2:
                retained_groups[4] += 1
            elif h == 1:
                retained_groups[5] += 1
            elif h == 0:
                retained_groups[6] += 1
        canonical_retained = tuple(
            sorted((group, count) for group, count in retained_groups.items() if count)
        )
        return mass, n_ar, canonical_retained, radicals, (), adjacency, edges, h_count, aromatic_atoms
    if r1_class in {"J_para", "J_ortho"}:
        return mass, n_ar, (), radicals, (), adjacency, edges, h_count, aromatic_atoms

    # Aromatic systems in the correlation domain are isolated six-rings.  This
    # check deliberately follows radical and J routing.
    aromatic_adjacency = {atom: set() for atom in aromatic_atoms}
    for edge in aromatic_edges:
        a, b = tuple(edge)
        aromatic_adjacency[a].add(b)
        aromatic_adjacency[b].add(a)
    unseen = set(aromatic_atoms)
    n_ar = 0
    while unseen:
        seed = unseen.pop()
        component = {seed}
        stack = [seed]
        while stack:
            for neighbour in aromatic_adjacency[stack.pop()]:
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    component.add(neighbour)
                    stack.append(neighbour)
        if len(component) != 6 or any(len(aromatic_adjacency[atom]) != 2 for atom in component):
            _reject("E_DOMAIN_FUSED_BIARYL")
        n_ar += 1
    for a, b, _ in edges:
        if a in aromatic_atoms and b in aromatic_atoms and frozenset((a, b)) not in aromatic_edges:
            _reject("E_DOMAIN_FUSED_BIARYL")

    try:
        cycles = [tuple(ring) for ring in molecule.get_smallest_set_of_smallest_rings()]
    except (ValueError, RuntimeError, TypeError):
        cycles = []
    nonaromatic_cycles = [ring for ring in cycles if not set(ring) <= aromatic_atoms]
    if r1_class == "quinoid_free":
        if any(len(ring) != 6 or not set(ring) <= carbon_set for ring in nonaromatic_cycles):
            _reject("E_DOMAIN_ALIPHATIC_RING")
        r1_ring = frozenset(atom for ring in nonaromatic_cycles for atom in ring)
    else:
        r1_ring = frozenset()
        if nonaromatic_cycles:
            _reject("E_DOMAIN_ALIPHATIC_RING")
    if len(edges) - len(carbons) + 1 != n_ar + len(nonaromatic_cycles):
        _reject("E_DOMAIN_ALIPHATIC_RING")

    double_edges = [(a, b) for a, b, order in edges if order == 2.0]
    double_by_atom = Counter(atom for edge in double_edges for atom in edge)
    if any(count > 1 for count in double_by_atom.values()) and r1_class != "quinoid_free":
        _reject("E_DOMAIN_CUMULATED_ALKENE")
    conjugated_bonds = [
        (a, b)
        for a, b, order in edges
        if order == 1.0 and a in double_by_atom and b in double_by_atom
    ]
    if conjugated_bonds and r1_class != "quinoid_free":
        _reject("E_DOMAIN_CONJUGATED_ALKENE")

    groups = Counter()
    conjugated_double_atoms = set()
    if r1_class == "quinoid_free":
        cycle_sets = [set(ring) for ring in nonaromatic_cycles]
        for a, b in conjugated_bonds:
            left = next(edge for edge in double_edges if a in edge)
            right = next(edge for edge in double_edges if b in edge)
            conjugated_double_atoms.update(left)
            conjugated_double_atoms.update(right)
            four_atoms = set(left + right)
            cyclic = any(four_atoms <= ring for ring in cycle_sets)
            groups[88 if cyclic else 89] += 1

    alkene_atoms = set(atom for edge in double_edges for atom in edge)
    for a, b in double_edges:
        if a in conjugated_double_atoms or b in conjugated_double_atoms:
            continue
        if a in r1_ring and b in r1_ring:
            groups[62] += 1
            continue
        heavy_a = len(adjacency[a])
        heavy_b = len(adjacency[b])
        aryl_attachments = sum(neighbour in aromatic_atoms for neighbour, _ in adjacency[a] if neighbour is not b)
        aryl_attachments += sum(neighbour in aromatic_atoms for neighbour, _ in adjacency[b] if neighbour is not a)
        if heavy_a == 1 or heavy_b == 1:
            groups[61] += 1
        elif aryl_attachments:
            groups[59] += aryl_attachments
        else:
            groups[58] += 1

    for atom in carbons:
        if atom in alkene_atoms or atom.radical_electrons:
            continue
        h = h_count[atom]
        if atom in r1_ring:
            if h not in {0, 1, 2}:
                _reject("E_DOMAIN_PRIMARY_GROUP")
            groups[{2: 9, 1: 10, 0: 11}[h]] += 1
            continue
        if atom in aromatic_atoms:
            groups[15 if h == 1 else 16] += 1
            continue
        aromatic_attachments = sum(neighbour in aromatic_atoms for neighbour, _ in adjacency[atom])
        if h == 3 and aromatic_attachments:
            groups[3] += 1
        elif aromatic_attachments:
            groups[8] += aromatic_attachments
        elif h in {3, 4}:
            groups[1] += 1
        elif h == 2:
            groups[4] += 1
        elif h == 1:
            groups[5] += 1
        elif h == 0:
            groups[6] += 1
        else:
            _reject("E_DOMAIN_PRIMARY_GROUP")

    def cx4(atom):
        return (
            atom not in aromatic_atoms
            and atom.radical_electrons == 0
            and all(order == 1.0 for _, order in adjacency[atom])
            and h_count[atom] + len(adjacency[atom]) == 4
        )

    def cx3(atom):
        return (
            atom not in aromatic_atoms
            and atom.radical_electrons == 0
            and sum(order == 2.0 for _, order in adjacency[atom]) == 1
            and h_count[atom] + len(adjacency[atom]) == 3
        )

    def q(atom):
        return cx4(atom) and h_count[atom] == 0

    def ch(atom):
        return cx4(atom) and h_count[atom] == 1

    for atom in carbons:
        if not q(atom):
            continue
        aromatic_neighbours = [neighbour for neighbour, _ in adjacency[atom] if neighbour in aromatic_atoms]
        aliphatic_neighbours = [neighbour for neighbour, _ in adjacency[atom] if neighbour not in aromatic_atoms]
        groups[130] += int(
            len(aromatic_neighbours) == 1
            and len(aliphatic_neighbours) == 3
            and all(cx4(neighbour) for neighbour in aliphatic_neighbours)
        )

    for a, b, order in edges:
        if order != 1.0:
            continue
        other_a = [neighbour for neighbour, _ in adjacency[a] if neighbour is not b]
        other_b = [neighbour for neighbour, _ in adjacency[b] if neighbour is not a]
        source_allowed = all(cx4(atom) or cx3(atom) for atom in other_a + other_b)
        groups[131] += int(ch(a) and ch(b) and source_allowed)
        groups[132] += int(((ch(a) and q(b)) or (q(a) and ch(b))) and source_allowed)
        groups[133] += int(q(a) and q(b) and source_allowed)

    for ring in aromatic_rings:
        ring_set = set(ring)
        substituted = [atom for atom in ring if any(neighbour not in ring_set for neighbour, _ in adjacency[atom])]
        pair_counts = {1: 0, 2: 0, 3: 0}
        for index, a in enumerate(substituted):
            distance = {a: 0}
            queue = deque([a])
            while queue:
                atom = queue.popleft()
                for neighbour in aromatic_adjacency[atom]:
                    if neighbour in ring_set and neighbour not in distance:
                        distance[neighbour] = distance[atom] + 1
                        queue.append(neighbour)
            for b in substituted[index + 1 :]:
                d = min(distance[b], 6 - distance[b])
                pair_counts[d] += 1
        m127, m128, m129 = pair_counts[1], pair_counts[2], 2 * pair_counts[3]
        groups[127] += int(m127 == 1 and m128 == 0 and m129 == 0)
        groups[128] += int(m127 == 0 and m128 in {1, 3} and m129 == 0)
        groups[129] += int(m127 == 0 and m128 == 0 and m129 == 2)

    canonical_groups = tuple(sorted((group, count) for group, count in groups.items() if count))
    return mass, n_ar, canonical_groups, radicals, r1_ring, adjacency, edges, h_count, aromatic_atoms


def evaluate(molecule, T, *, off_above_1500=False, r1_class=None):
    """Evaluate the frozen volatility law for one connected RMG molecule."""

    try:
        if not isinstance(off_above_1500, (bool, np.bool_)):
            _reject("E_TYPE_OPTION")
        if r1_class is not None and not isinstance(r1_class, str):
            _reject("E_TYPE_R1_CLASS")
        if r1_class not in {None, "J_para", "J_ortho", "quinoid_free"}:
            _reject("E_RANGE_R1_CLASS")
        if isinstance(T, (bool, np.bool_)) or not isinstance(
            T, (int, float, np.integer, np.floating)
        ):
            _reject("E_TYPE_T")
        T = float(T)
        if not math.isfinite(T):
            _reject("E_NONFINITE_T")
        if T <= 0.0 or not math.isfinite(R * T) or R * T == 0.0:
            _reject("E_RANGE_T")

        M, n_ar, groups_tuple, radicals, _, _, _, _, _ = _compile_graph(molecule, r1_class)
        flags = ["extrapolated_M"] if M > M_C36 else []
        if radicals:
            return VolatilityOutcome(
                "retained_radical",
                "-",
                M,
                n_ar,
                groups_tuple,
                regime="radical_prohibited",
                flags=tuple(flags + ["radical_prohibited"]),
            )
        if r1_class in {"J_para", "J_ortho"}:
            return VolatilityOutcome(
                "retained_J_ring",
                "-",
                M,
                n_ar,
                (),
                regime=r1_class,
                flags=tuple(flags + ["R1_J_ring"]),
            )
        if r1_class == "quinoid_free":
            flags += ["R1_quinoid", "structure_extrapolated"]

        effective_mass = M + 14.2 * n_ar
        power = effective_mass ** (2.0 / 3.0)
        Tb = 1070.0 - math.exp(6.98291 - 0.02013 * power)
        sum_tc = math.fsum(_TC[group] * count for group, count in groups_tuple)
        sum_pc = math.fsum(_PC[group] * count for group, count in groups_tuple)
        if sum_tc <= 0.0 or 0.0093898 + sum_pc <= 0.0:
            _reject("E_CRITICAL_CORRELATION_DOMAIN")
        Tc = Tb * (1.0 / (0.988948 + sum_tc ** 0.86074) + 0.699003)
        Pc = (M ** -0.140414) / 100.0 / (0.0093898 + sum_pc) ** 2
        if not all(math.isfinite(value) and value > 0.0 for value in (Tb, Tc, Pc)):
            _reject("E_CRITICAL_CORRELATION_DOMAIN")
        if T >= Tc:
            return VolatilityOutcome(
                "route_gas",
                "-",
                M,
                n_ar,
                groups_tuple,
                Tb,
                Tc,
                Pc,
                regime="critical",
                flags=tuple(flags + ["above_critical"]),
            )
        if off_above_1500 and M > M_OFF:
            return VolatilityOutcome(
                "evaporation_disabled",
                "-",
                M,
                n_ar,
                groups_tuple,
                Tb,
                Tc,
                Pc,
                regime="off_arm",
                flags=tuple(flags + ["evaporation_disabled"]),
            )

        tau = max(0.0, (M - 76.1 * n_ar) / 14.03 - 3.0 + 0.5 * n_ar)
        delta_s = 86.0 + 0.4 * tau
        if _below_tb(T, Tb):
            delta_cp = -90.0 - 2.1 * tau
            log_p_atm = math.fsum(
                (
                    -delta_s * (Tb - T) / (R * T),
                    delta_cp / R * ((Tb - T) / T - math.log(Tb / T)),
                )
            )
            regime = "below_Tb"
        else:
            log_p_atm = delta_s / R * (1.0 - Tb / T)
            regime = "above_Tb"
        log10p = (log_p_atm + math.log(ATM_TO_BAR)) / math.log(10.0)
        if not math.isfinite(log10p):
            _reject("E_NUMERIC_RANGE")
        if log10p > math.log10(Pc):
            _reject("E_PSAT_GT_PC")
        return VolatilityOutcome(
            "P_sat",
            "-",
            M,
            n_ar,
            groups_tuple,
            Tb,
            Tc,
            Pc,
            log10p,
            regime,
            tuple(flags),
        )
    except _VolatilityError as exc:
        return _error(str(exc))
    except (OverflowError, ValueError, ZeroDivisionError):
        return _error("E_NUMERIC_RANGE")


__all__ = ["VolatilityOutcome", "evaluate"]
