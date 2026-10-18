"""The declared v1 polystyrene domain, and the separated-radical spin policy.

``check(molecule)`` accepts a molecule (one connected component, RMG
``Molecule``; atom ids are optional) only if it is one of

* a **linear PS chain**: a carbon backbone that is a simple path, phenyl
  (C6H5) pendants on backbone carbons (at most one per carbon, attached by a
  single bond), every other backbone substituent H;
* the declared small species ``H`` (atom) and ``H2``.

On a chain it requires

* a **saturated backbone**: every backbone-backbone bond single, except the
  bond at a chain end, which may be double (``END_GROUPS``: vinylidene
  ``CH2=C(Ph)-`` and styryl ``PhCH=CH-``, the unsaturated ends beta-scission
  makes; no other unsaturation, no triple bonds, no cumulated ends);
* **declared end groups** (``END_GROUPS``) at both termini;
* **declared radical sites** (``RADICAL_KINDS``): benzylic (a backbone carbon
  bearing a phenyl), primary end (a terminal CH2), secondary backbone (an
  interior CH without phenyl), methyl (CH3, a one-carbon chain); never on a
  ring carbon, a double-bonded carbon or next to one; one electron per site;
* at most ``MAX_RADICAL_SITES`` (R) sites per chain, every pair farther apart
  than ``RADICAL_SEPARATION`` (S) heavy-atom bonds and not in one pi system;
* the declared multiplicity convention (``MULTIPLICITY_CONVENTION``: maximum)
  on the input AS GIVEN, before any normalization: a separated biradical
  given as a singlet is refused, never silently rebuilt as a triplet.

An unsaturated end is read on both alkene carbons: styryl ``PhCH=CH-`` is a
terminal CH(Ph) double-bonded to a CH without phenyl; vinylidene
``CH2=C(Ph)-`` a terminal CH2 double-bonded to a phenyl-bearing carbon
(H counts and no radical on either carbon).

``check`` is a verdict on one species.  Certificates consume it as the
channel's ``eligibility`` disposition (``rmgpy.kmc.certificates``), distinct
from the spin/pathway ``support`` of ``rmgpy.kmc.channels``.

Everything else is ``refused(reason)``.  The chain is read on the given form
and, if that is refused, on its most aromatic reactive resonance form (RMG's
``generate_resonance_structures``; the least canonical SMILES among the forms
with the most benzene bonds), so a Kekule or quinoid input of a declared
species is accepted.

Constants (reported, see the I-135 return):

* ``MAX_RADICAL_SITES = 2``: the S2 acceptance states carry chains with two
  radical sites (two mid-chain benzylic radicals; a primary head and a
  benzylic tail); no reference has three, so a triradical is refused rather
  than certified by extrapolation.
* ``RADICAL_SEPARATION = 4``: a pair of sites four or fewer heavy bonds apart
  is refused.  Four heavy bonds (C-C-C-C-C) is the reach of the 1,5 H shift
  (intra_H_migration, the longest-range unimolecular S2 class), so two sites
  that close can react with each other (intramolecular recombination,
  disproportionation, H shift between them) and a doublet + doublet picture
  is not a declared approximation there; it is also beyond the through-bond
  exchange range of saturated sp3 linkers.  It is NOT the reach of the RMG
  tree predicates: those read the whole-molecule multiplicity, and
  R_Recombination's list-labelled predicates search the whole merged
  structure, which the rate-selection certificate records as global guards
  and remote candidates instead of hiding them behind a distance.

Separated-radical policy (``spin_verdict``; a declared approximation, not
physics): a species with two or more radical electrons is supported if no
atom carries two (``u2``) and every pair of radical sites is more than
``RADICAL_SEPARATION`` heavy bonds apart and not in one pi system; the
molecule is then built with the MAXIMUM multiplicity convention
(``MULTIPLICITY_CONVENTION``), which the certificate records as a global
guard.  Closer or conjugated pairs and any ``u2`` atom are
``unsupported:spin``.  The verdict is per species: two doublets on different
participants are never refused for that reason.
"""

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Optional, Tuple

MAX_RADICAL_SITES = 2       # R
RADICAL_SEPARATION = 4      # S: a pair of sites must be more than S heavy bonds apart
MULTIPLICITY_CONVENTION = ("maximum: multiplicity = number of unpaired electrons + 1 "
                           "(rmg_moves.molecule_from_rows -> Molecule.update_multiplicity; "
                           "RMG's own SMILES/adjacency default for a radical without a "
                           "multiplicity line)")

END_GROUPS = {
    "methyl": "CH3- (H-capped saturated end; ethyl-type when its neighbour is CH2)",
    "benzyl": "PhCH2- (H-capped benzylic end)",
    "primary_radical": ".CH2- (primary end radical)",
    "benzylic_radical": "PhCH.- (benzylic end radical)",
    "vinylidene": "CH2=C(Ph)- (unsaturated end from beta-scission of a benzylic radical)",
    "styryl": "PhCH=CH- (unsaturated end from beta-scission of a secondary radical)",
}
RADICAL_KINDS = ("benzylic", "primary_end", "secondary_backbone", "methyl", "h_atom")
SMALL_SPECIES = ("H", "H2")
SUPPORTED = "supported"
SPIN_UNSUPPORTED = "unsupported:spin"


@dataclass(frozen=True)
class DomainVerdict:
    accepted: bool
    reason: str
    kind: str = ""
    details: Dict = field(default_factory=dict, compare=False)

    def as_dict(self):
        return {"accepted": self.accepted, "reason": self.reason, "kind": self.kind,
                "details": self.details}


def refused(reason: str, **details) -> DomainVerdict:
    return DomainVerdict(False, reason, "refused", details)


def _accepted(kind, **details) -> DomainVerdict:
    return DomainVerdict(True, "", kind, details)


@dataclass(frozen=True)
class SpinVerdict:
    status: str                       # SUPPORTED | SPIN_UNSUPPORTED
    reason: str
    radical_electrons: int
    sites: Tuple[int, ...] = ()       # atom positions (ids when present) of the radical sites
    pairs: Tuple[Tuple, ...] = ()     # (site, site, heavy distance, same pi system)
    convention: Optional[str] = None  # the multiplicity convention when the species is a biradical+

    @property
    def supported(self):
        return self.status == SUPPORTED


# ---------------------------------------------------------------------------
# representation-independent form
# ---------------------------------------------------------------------------
def _key(atom):
    return atom.id if isinstance(getattr(atom, "id", None), int) and atom.id >= 0 else None


def _benzene_bonds(mol) -> int:
    return sum(1 for b in mol.get_all_edges() if b.is_benzene())


def canonical_form(molecule):
    """The most aromatic reactive resonance form (ties: least SMILES), with the
    input atom ids copied when RMG preserved them."""
    from rmgpy.molecule.resonance import generate_resonance_structures
    base = molecule.copy(deep=True)
    for atom, original in zip(base.atoms, molecule.atoms):
        atom.id = original.id
    if len(base.atoms) < 3:
        return base
    forms = generate_resonance_structures(base, keep_isomorphic=False, filter_structures=True)
    forms = [f for f in forms if f.reactive] or list(forms)
    best = min(forms, key=lambda f: (-_benzene_bonds(f), f.to_smiles()))
    return best


def _heavy_distances(mol, source):
    dist = {source: 0}
    frontier = [source]
    while frontier:
        nxt = []
        for a in frontier:
            for b in a.bonds:
                if b not in dist and not b.is_hydrogen():
                    dist[b] = dist[a] + 1
                    nxt.append(b)
        frontier = nxt
    return dist


def _pi_components(mol):
    """Atom -> component index of the pi subgraph: radical atoms and atoms with
    a bond of order > 1 (double, triple, benzene), joined by any bond."""
    pi = {a for a in mol.atoms if a.radical_electrons > 0 or
          any(b.order > 1.0 for b in a.bonds.values())}
    comp = {}
    k = 0
    for start in mol.atoms:
        if start not in pi or start in comp:
            continue
        stack = [start]
        comp[start] = k
        while stack:
            a = stack.pop()
            for b in a.bonds:
                if b in pi and b not in comp:
                    comp[b] = k
                    stack.append(b)
        k += 1
    return comp


# ---------------------------------------------------------------------------
# spin policy
# ---------------------------------------------------------------------------
def spin_verdict(molecule, separation: int = RADICAL_SEPARATION) -> SpinVerdict:
    """Separated-radical policy for one species (see the module docstring).

    Representation independent without enumerating resonance forms: a
    radical electron moves only within its conjugated (pi) system, which is
    the same in every resonance form (``_pi_components``: radical atoms and
    atoms with a multiple or benzene bond, joined by bonds).  Two sites in one
    pi system are refused.  Otherwise the separation of a pair is the least
    heavy-bond distance between their two pi systems (a lower bound of the
    site distance in every form, so the test is conservative)."""
    total = sum(int(a.radical_electrons) for a in molecule.atoms)
    if any(a.radical_electrons >= 2 for a in molecule.atoms):
        return SpinVerdict(SPIN_UNSUPPORTED, "an atom carries two unpaired electrons (u2)", total)
    if total <= 1:
        return SpinVerdict(SUPPORTED, "", total)
    sites = [a for a in molecule.atoms if a.radical_electrons]
    comp = _pi_components(molecule)
    members = defaultdict(list)
    for atom, k in comp.items():
        members[k].append(atom)
    pairs = []
    reason = ""
    for i, a in enumerate(sites):
        for b in sites[i + 1:]:
            same_pi = comp[a] == comp[b]
            if same_pi:
                d = 0
            else:
                d = None
                for x in members[comp[a]]:
                    if x.is_hydrogen():
                        continue
                    dist = _heavy_distances(molecule, x)
                    for y in members[comp[b]]:
                        if y in dist and (d is None or dist[y] < d):
                            d = dist[y]
            pairs.append((_site_key(molecule, a), _site_key(molecule, b), d, same_pi))
            if not reason and same_pi:
                reason = "radical sites in one pi system"
            if not reason and (d is None or d <= separation):
                reason = (f"radical pi systems {d} heavy bonds apart "
                          f"(must be > {separation})")
    site_keys = tuple(_site_key(molecule, a) for a in sites)
    if reason:
        return SpinVerdict(SPIN_UNSUPPORTED, reason, total, site_keys, tuple(pairs))
    return SpinVerdict(SUPPORTED, "", total, site_keys, tuple(pairs), MULTIPLICITY_CONVENTION)


def _site_key(molecule, atom):
    k = _key(atom)
    return k if k is not None else molecule.atoms.index(atom)


# ---------------------------------------------------------------------------
# domain checker
# ---------------------------------------------------------------------------
def _phenyl_rings(form):
    """(ring atom sets, ipso atom -> backbone neighbour) or a refusal string."""
    rings = form.get_smallest_set_of_smallest_rings() if any(
        form.is_atom_in_cycle(a) for a in form.atoms) else []
    seen = set()
    ipso_of = {}
    for ring in rings:
        ring = set(ring)
        if len(ring) != 6 or any(not a.is_carbon() for a in ring):
            return None, "a ring that is not a six-carbon ring"
        if ring & seen:
            return None, "fused rings"
        seen |= ring
        orders = []
        for a in ring:
            for b, bond in a.bonds.items():
                if b in ring and form.atoms.index(a) < form.atoms.index(b):
                    orders.append(bond.order)
        if not (all(o == 1.5 for o in orders) or sorted(orders) == [1.0] * 3 + [2.0] * 3):
            return None, "a six-ring that is not benzenoid"
        outside = []
        for a in ring:
            if a.radical_electrons:
                return None, "a radical on a ring carbon in the most aromatic form"
            subs = [b for b in a.bonds if b not in ring]
            if len(subs) != 1:
                return None, "a ring carbon without exactly one substituent"
            if not subs[0].is_hydrogen():
                outside.append((a, subs[0]))
        if len(outside) != 1:
            return None, "a ring that is not a monosubstituted phenyl"
        ipso, sub = outside[0]
        if not sub.is_carbon():
            return None, "a phenyl not attached to carbon"
        if ipso.bonds[sub].order != 1.0:
            return None, "a phenyl not attached by a single bond"
        ipso_of[ipso] = sub
    return (seen, ipso_of), ""


# neutral valence states: element -> {(lone pairs, bond-order sum + unpaired
# electrons)}; hypervalent S (sulfoxide S4, sulfone S6) and P (P5) included.
# A singlet carbene (C, 1 pair, valence 2) or a divalent S with one pair has
# no state; charges and u2 are refused by their own rules.
NEUTRAL_STATES = {"H": {(0, 1)}, "C": {(0, 4)}, "Si": {(0, 4)}, "N": {(1, 3)},
                  "P": {(1, 3), (0, 5)}, "O": {(2, 2)}, "S": {(2, 2), (1, 4), (0, 6)},
                  "F": {(3, 1)}, "Cl": {(3, 1)}, "Br": {(3, 1)}, "I": {(3, 1)}}


def in_neutral_state(symbol: str, lone_pairs: int, valence: float) -> bool:
    """Whether (``lone_pairs``, ``valence``) is a ``NEUTRAL_STATES`` state of
    ``symbol`` (an element without an entry has none)."""
    return (int(lone_pairs), valence) in NEUTRAL_STATES.get(symbol, ())


# aromatic atoms (one or more benzene bonds, B): the neutral aromatic atom
# types of ``rmgpy.molecule.atomtype`` (``get_atomtype`` decides the type;
# no type, or a type not listed, is a refusal), RMG's own charge consistency
# (``adjlist.ConsistencyChecker.check_partial_charge``: valence electrons -
# total bond order - unpaired electrons - 2 lone pairs equals the stored
# charge, 0; a B bond counts 3/2, or 4/3 on an atom with three) and, per type,
# the bond counts allowed: (lone pairs, B bonds, single, double, triple bonds,
# unpaired electrons), never summed.  Each type's set is exactly its neutral,
# charge-consistent states (RMG's types do not read unpaired electrons, so
# the counts carry them):
#   Cb, Sib     two B + one single bond, or two B + one unpaired electron
#   Cbf, Sibf   three B
#   N5b, P5b    two B + one single bond + one unpaired electron (with the
#               single bond alone the atom is formally +1)
#   N3b, P3b    two B + one lone pair
#   N5bd, P5bd  two B + one double bond
# Not listed: O4b and S4b have no neutral charge-consistent state (two B and
# one lone pair is formally +1: RMG's adjacency-list reader refuses it, and
# RMG's resonance never writes a five-membered ring with B bonds, so furan
# and thiophene reach this check in their Kekule form); Cb with one B bond
# (a single B bond is no aromatic ring).  A bond on an aromatic atom whose
# order is not single, double, triple or B (a quadruple bond, which Cb
# leaves free) is refused before the type is read.
_RING_ONE = frozenset({(0, 2, 1, 0, 0, 0), (0, 2, 0, 0, 0, 1)})
_RING_RADICAL = frozenset({(0, 2, 1, 0, 0, 1)})
_FUSED = frozenset({(0, 3, 0, 0, 0, 0)})
_PAIR = frozenset({(1, 2, 0, 0, 0, 0)})
_DOUBLE = frozenset({(0, 2, 0, 1, 0, 0)})
AROMATIC_TYPES = {"Cb": _RING_ONE, "Cbf": _FUSED, "Sib": _RING_ONE, "Sibf": _FUSED,
                  "N3b": _PAIR, "N5b": _RING_RADICAL, "N5bd": _DOUBLE,
                  "P3b": _PAIR, "P5b": _RING_RADICAL, "P5bd": _DOUBLE}


def valence_state(atom) -> tuple:
    """(lone pairs, bond-order sum + unpaired electrons) of an atom without a
    benzene bond; the sum rounded to a half."""
    total = sum(float(bond.order) for bond in atom.bonds.values()) + atom.radical_electrons
    return int(atom.lone_pairs), round(total * 2) / 2


def aromatic_state(atom) -> tuple:
    """(lone pairs, B bonds, single, double, triple bonds, unpaired electrons)
    of an atom with a benzene bond (``AROMATIC_TYPES``); counts, never a
    sum."""
    bonds = list(atom.bonds.values())
    return (int(atom.lone_pairs), sum(1 for b in bonds if b.is_benzene()),
            sum(1 for b in bonds if b.is_single()), sum(1 for b in bonds if b.is_double()),
            sum(1 for b in bonds if b.is_triple()), int(atom.radical_electrons))


def aromatic_refusal(atom) -> Optional[str]:
    """Why an atom with a benzene bond is outside ``AROMATIC_TYPES`` (None:
    inside): only single, double, triple and B bonds; RMG's atom type
    (``get_atomtype``) listed; RMG's charge consistency
    (``check_partial_charge``); and the atom's bond counts allowed for its
    type."""
    from rmgpy.molecule.adjlist import ConsistencyChecker, InvalidAdjacencyListError
    from rmgpy.molecule.atomtype import AtomTypeError, get_atomtype
    symbol = atom.element.symbol
    other = sorted({float(b.order) for b in atom.bonds.values()
                    if not (b.is_single() or b.is_double() or b.is_triple() or b.is_benzene())})
    if other:
        return (f"an aromatic atom ({symbol}) with a bond of order {other[0]:g}: only single, "
                f"double, triple and B bonds are read")
    try:
        label = get_atomtype(atom, atom.edges).label
    except AtomTypeError as error:
        return f"an aromatic atom ({symbol}) with no RMG atom type: {error}"
    try:
        ConsistencyChecker.check_partial_charge(atom)
    except InvalidAdjacencyListError as error:
        return (f"an aromatic atom ({symbol}) of RMG atom type {label} fails RMG's charge "
                f"consistency (formal charge != stored charge {atom.charge}): {error}")
    state = aromatic_state(atom)
    if label not in AROMATIC_TYPES:
        return f"an aromatic atom ({symbol}) of RMG atom type {label}: not a neutral " \
               f"aromatic type"
    if state not in AROMATIC_TYPES[label]:
        return (f"an aromatic atom ({symbol}) of RMG atom type {label} with bond counts "
                f"{state} (lone pairs, B, single, double, triple, unpaired): not "
                f"{sorted(AROMATIC_TYPES[label])}")
    return None


def _element_list(elements) -> str:
    elements = list(elements)
    return elements[0] if len(elements) == 1 else \
        ", ".join(elements[:-1]) + " and " + elements[-1]


def universal_check(molecule, elements=("C", "H"), *, states=None,
                    hydrogen_partner: str = "carbon") -> Optional[DomainVerdict]:
    """The species rules every domain shares, in ``check``'s order: only
    ``elements``; no charge; lone pairs (``states``: element -> the allowed
    ``valence_state`` set, ``NEUTRAL_STATES``; None: no lone pair at all);
    no u2 atom; the multiplicity convention on the input
    as given; hydrogen-only species are H or H2 (an accepted verdict); every
    H bonded to exactly one ``hydrogen_partner`` (``"carbon"`` or
    ``"heavy"``).  None: the rest is the caller's."""
    atoms = molecule.atoms
    if any(a.element.symbol not in elements for a in atoms):
        return refused(f"an element other than {_element_list(elements)}")
    if any(a.charge != 0 for a in atoms):
        return refused("a charged atom")
    if states is None:
        if any(a.lone_pairs != 0 for a in atoms):
            return refused("an atom with a lone pair")
    else:
        for a in atoms:
            if any(bond.is_benzene() for bond in a.bonds.values()):
                why = aromatic_refusal(a)
                if why:
                    return refused(why)
                continue
            state = valence_state(a)
            if state not in states.get(a.element.symbol, ()):
                return refused(f"an atom ({a.element.symbol}) with {state[0]} lone pairs and "
                               f"valence {state[1]:g}: not a neutral valence state "
                               f"{sorted(states.get(a.element.symbol, ()))}")
    if any(a.radical_electrons >= 2 for a in atoms):
        return refused("an atom with two unpaired electrons (u2)")
    # the declared convention is checked on the input as given, BEFORE any
    # normalization (canonical_form, molecule_from_rows) rebuilds it at
    # maximum multiplicity and so could change what RMG reads
    unpaired = sum(int(a.radical_electrons) for a in atoms)
    if int(molecule.multiplicity) != unpaired + 1:
        return refused(f"multiplicity convention violated: multiplicity "
                       f"{molecule.multiplicity} with {unpaired} unpaired electrons "
                       f"(declared: {MULTIPLICITY_CONVENTION.split(':')[0]}, "
                       f"{unpaired + 1})", multiplicity=int(molecule.multiplicity))
    if all(a.is_hydrogen() for a in atoms):
        if len(atoms) == 1 and atoms[0].radical_electrons == 1:
            return _accepted("small:H", radical_sites=[(_key(atoms[0]), "h_atom")])
        if len(atoms) == 2 and not any(a.radical_electrons for a in atoms):
            return _accepted("small:H2")
        return refused("a hydrogen-only species other than H and H2")
    for a in atoms:
        if not a.is_hydrogen():
            continue
        partner_ok = len(a.bonds) == 1 and (
            next(iter(a.bonds)).is_carbon() if hydrogen_partner == "carbon"
            else not next(iter(a.bonds)).is_hydrogen())
        if not partner_ok:
            return refused(f"a hydrogen not bonded to exactly one {hydrogen_partner}"
                           if hydrogen_partner != "heavy"
                           else "a hydrogen not bonded to exactly one heavy atom")
    return None


def check(molecule, separation: int = RADICAL_SEPARATION,
          max_sites: int = MAX_RADICAL_SITES) -> DomainVerdict:
    """The v1 PS domain verdict of one molecule (see the module docstring)."""
    prefix = universal_check(molecule, ("C", "H"))
    if prefix is not None:
        return prefix
    first = _check_form(molecule, separation, max_sites)
    if first.accepted:
        return first
    return _check_form(canonical_form(molecule), separation, max_sites)


def declared_interim_check(molecule, elements=("C", "H")) -> DomainVerdict:
    """INTERIM domain verdict for declared chains (``windows.DECLARED_INTERIM``,
    policy ``declared_interim_v1``).  A named placeholder: a later I-165
    slice (the declared chain domain) replaces it.

    Accepts a molecule iff every element is in ``elements`` (C, H and the
    elements of the witness's own declarations: ``windows.policy_for`` binds
    them; never read from the provider registry); no atom is charged or u2;
    the multiplicity convention holds on the input as given; every atom is
    in a neutral valence state (``NEUTRAL_STATES``: lone pairs and
    bond-order sum, hypervalent S and P included); every H is bonded to exactly one heavy atom and a
    hydrogen-only species is H or H2; at most ``MAX_RADICAL_SITES`` unpaired
    electrons; and ``spin_verdict`` is supported when there are two or
    more.  The verdict kind is ``declared_interim``; ``details`` gives the
    radical sites (atom ids).

    Deliberately NOT applied: ``_check_form``'s PS shape rules (simple-path
    backbone, phenyl pendants, declared end groups and radical kinds) and the
    radical separation beyond ``spin_verdict``.  That widening is the point:
    a branched PP or an ester-bearing PMMA chain is in this domain."""
    prefix = universal_check(molecule, tuple(sorted(set(elements) | {"C", "H"})),
                             states=NEUTRAL_STATES, hydrogen_partner="heavy")
    if prefix is not None and not prefix.accepted:
        return prefix
    sites = [_key(a) for a in molecule.atoms if a.radical_electrons]
    unpaired = sum(int(a.radical_electrons) for a in molecule.atoms)
    if unpaired > MAX_RADICAL_SITES:
        return refused(f"{unpaired} unpaired electrons > MAX_RADICAL_SITES "
                       f"({MAX_RADICAL_SITES})")
    if unpaired >= 2:
        spin = spin_verdict(molecule)
        if not spin.supported:
            return refused(f"spin: {spin.reason}")
    return _accepted("declared_interim", radical_sites=sites)


def _check_form(form, separation, max_sites) -> DomainVerdict:
    """The chain verdict read on one resonance form."""
    if any(a.radical_electrons >= 2 for a in form.atoms):
        return refused("an atom with two unpaired electrons (u2) in this resonance form")
    rings, why = _phenyl_rings(form)
    if rings is None:
        return refused(why)
    ring_atoms, ipso_of = rings
    backbone = [a for a in form.atoms if a.is_carbon() and a not in ring_atoms]
    if not backbone:
        return refused("no backbone carbon (an arene, not a chain)")
    bset = set(backbone)
    bb = {a: [b for b in a.bonds if b in bset] for a in backbone}
    edges = sum(len(v) for v in bb.values()) // 2
    if any(len(v) > 2 for v in bb.values()):
        return refused("a branched backbone (a carbon with three backbone neighbours)")
    if edges != len(backbone) - 1:
        return refused("the backbone is not one simple path")
    pendant = {a: [i for i in ipso_of if ipso_of[i] is a] for a in backbone}
    if any(len(v) > 1 for v in pendant.values()):
        return refused("two phenyls on one backbone carbon")
    termini = [a for a in backbone if len(bb[a]) <= 1]
    # walk the path from one terminus
    order = [termini[0]]
    while len(order) < len(backbone):
        nxt = [b for b in bb[order[-1]] if b not in order]
        order.append(nxt[0])
    for a, b in zip(order, order[1:]):
        o = a.bonds[b].order
        if o == 1.0:
            continue
        if o != 2.0:
            return refused("a backbone bond that is neither single nor double")
        if not (a is order[0] or b is order[-1]):
            return refused("backbone unsaturation away from a chain end")
    for a in backbone:
        if sum(1 for b, bond in a.bonds.items() if b in bset and bond.order == 2.0) > 1:
            return refused("a cumulated double bond")
    ends = []
    for t in ([order[0], order[-1]] if len(order) > 1 else [order[0]]):
        end = _end_group(t, bb, pendant)
        if end is None:
            return refused(f"an undeclared chain-end group at a backbone terminus")
        ends.append(end)
    if len(order) == 1:
        t = order[0]
        h = sum(1 for b in t.bonds if b.is_hydrogen())
        if pendant[t]:
            kind = "benzylic_radical" if t.radical_electrons else "benzyl"
        elif h == 4:
            kind = "methane"
        elif h == 3 and t.radical_electrons == 1:
            kind = "methyl_radical"
        else:
            return refused("an undeclared one-carbon species")
        ends = [kind]
    sites = []
    for a in form.atoms:
        if not a.radical_electrons:
            continue
        if a not in bset:
            return refused("a radical site off the backbone")
        if any(bond.order > 1.0 for b, bond in a.bonds.items() if b in bset):
            return refused("a vinylic radical (radical on a double-bonded carbon)")
        if any(bond.order == 2.0 for b in bb[a] for _c, bond in b.bonds.items()):
            return refused("an allylic radical (next to a chain-end double bond)")
        if len(order) == 1:
            kind = "benzylic" if pendant[a] else "methyl"
        elif pendant[a]:
            kind = "benzylic"
        elif len(bb[a]) == 1:
            kind = "primary_end"
        else:
            kind = "secondary_backbone"
        sites.append((a, kind))
    if len(sites) > max_sites:
        return refused(f"{len(sites)} radical sites (at most {max_sites} per chain)")
    if len(sites) > 1:
        spin = spin_verdict(form, separation)
        if not spin.supported:
            return refused(f"separated-radical policy: {spin.reason}")
    return _accepted("ps_chain", backbone_carbons=len(backbone), phenyls=len(ipso_of),
                     ends=sorted(ends),
                     radical_sites=sorted(((_key(a), k) for a, k in sites),
                                          key=lambda x: (x[0] is None, x[0] or 0, x[1])))


def _end_group(t, bb, pendant) -> Optional[str]:
    """The declared end group at backbone terminus ``t`` (None if undeclared)."""
    h = sum(1 for b in t.bonds if b.is_hydrogen())
    rad = t.radical_electrons
    if len(bb[t]) == 0:
        return "single"
    (n,) = bb[t]
    order = t.bonds[n].order
    if order == 2.0:
        # both alkene carbons are checked: the neighbour's substituents, H
        # count (by valence: one double bond, its other backbone bond, a
        # phenyl if declared) and that it is no radical
        n_h = sum(1 for b in n.bonds if b.is_hydrogen())
        n_other = len(bb[n]) - 1          # backbone neighbours of n other than t
        if n.radical_electrons or rad:
            return None
        if pendant[t] and h == 1 and not pendant[n] and n_h == 2 - n_other:
            return "styryl"               # PhCH=CH-
        if not pendant[t] and h == 2 and pendant[n] and n_h == 1 - n_other:
            return "vinylidene"           # CH2=C(Ph)-
        return None
    if pendant[t]:
        if h == 2 and not rad:
            return "benzyl"
        if h == 1 and rad == 1:
            return "benzylic_radical"
        return None
    if h == 3 and not rad:
        return "methyl"
    if h == 2 and rad == 1:
        return "primary_radical"
    return None


def check_rows(rows) -> DomainVerdict:
    """``check`` of an atom-id -> row component (rmg_moves row shape)."""
    from rmgpy.kmc import rmg_moves as rm
    return check(rm.molecule_from_rows(rows))
