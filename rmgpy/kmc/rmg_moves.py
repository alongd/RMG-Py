"""RMG-backed mapped-move oracle: directed, atom-mapped elementary moves taken
from RMG's own family machinery, before any degeneracy collapse.

``generate_moves(family, reactants, direction)`` takes RMG ``Molecule``
reactants whose atoms carry stable integer ids (``atom.id``; the same ids as a
``MeltState`` built from them, see ``molecule_from_rows``) and returns every
template match RMG's family finds as one ``Move``.  It drives RMG's own
primitives (``_match_reactant_to_template``, ``is_molecule_forbidden``,
``apply_recipe``, ``_create_reaction``, ``get_reaction_template_labels``) in
the loop structure of ``KineticsFamily._generate_reactions``; it never calls
``_generate_reactions`` itself (which drops filtered mappings silently), and it
never calls ``find_degenerate_reactions`` or ``calculate_degeneracy``.  It
shares no code with ``window_matcher`` or the test-side transition oracle.

Pinned policies
---------------
* **One move per mapping.**  A mapping is (template slot order, the resonance
  form of every reactant, one subgraph isomorphism per template reactant as
  RMG returns it).  Nothing is merged: not by product isomorphism, not by
  template, not by degeneracy, not across resonance forms.  ``apply_recipe``
  is called with ``unique=False``; RMG threads ``unique`` into
  ``ReactionRecipe._apply`` but never reads it, so no mapped duplicate is
  dropped on the recipe path either way.
* **Reactant resonance** (``resonance=``): ``"all"`` (default) enumerates every
  REACTIVE form ``generate_resonance_structures(keep_isomorphic=True,
  filter_structures=True)`` returns (``keep_isomorphic=True`` is required: two
  isomorphic forms place electrons on different atom ids).  RMG also returns
  the input itself, marked not reactive, when the input is not one of its
  reactive forms (e.g. a Kekule benzyl); that extra form depends on how the
  reactant was written, so it is dropped and the enumerated set is the same
  for every starting form.  ``"given"`` uses the form passed in only.  Forms
  are ordered by their integer-keyed signature, so ``Move.forms`` indices are
  deterministic.
* **Product resonance**: a product is compared modulo resonance by
  ``product_key`` (canonical = least integer-keyed signature over the
  *reactive* forms RMG generates for each product molecule, atom ids kept) and
  ``same_product`` (RMG ``is_isomorphic`` with the identity initial map,
  ``strict=True``, across those forms).  Both use RMG's resonance code only.
* **Filters are statuses, never silent drops**: ``kept``, ``forbidden``
  (family forbidden structure on a labelled reactant or on a product),
  ``recipe_failed`` (``InvalidActionError``/``KekulizationError``/
  ``ActionError``, or ``apply_recipe`` returned no products: product count,
  net charge or product-multiplicity check), ``unchanged_species`` (products
  are the reactants), ``charged_species``, ``unbalanced``, ``excluded`` (the
  caller's ``select`` refused the mapping; the recipe was not applied).
  ``constraint`` is reserved: species / model-growth constraints are OFF
  (``fails_species_constraints`` is never called).
* **Two reactants** must be two different ``Molecule`` objects with disjoint
  atom ids, also for two copies of one species.  As in RMG, both template slot
  orders are tried, so a symmetric family (``R_Recombination``) yields every
  pair twice, as ``order=(0, 1)`` and ``order=(1, 0)``, with equal products;
  ``Disproportionation`` yields (acceptor A, donor B) and (acceptor B,
  donor A) as different moves with different products.
* **Pathway vs product identity**: ``Move.pathway`` = (family, template node
  labels); ``product_key(move)`` is the product identity.  They are kept apart.
* **Determinism**: moves are sorted by integer keys and carry no object
  identity; the output is identical across processes and ``PYTHONHASHSEED``.

``rmg_degeneracy`` reports RMG's final reaction degeneracy (the full
``generate_reactions_from_families`` pipeline) for comparison only.
"""

from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass
from typing import Callable, Dict, Iterable, Optional, Sequence, Tuple

from rmgpy.exceptions import ActionError, InvalidActionError, KekulizationError
from rmgpy.molecule.element import get_element
from rmgpy.molecule.molecule import Atom, Bond, Molecule
from rmgpy.molecule.resonance import generate_resonance_structures
from rmgpy.reaction import same_species_lists

STATUSES = ("kept", "forbidden", "recipe_failed", "unchanged_species", "charged_species",
            "unbalanced", "constraint", "excluded")
RESONANCE_POLICIES = ("all", "given")
DIRECTIONS = ("forward", "reverse")


class RmgMovesError(ValueError):
    """Bad input, or an RMG step lost the atom ids."""


class AtomIdError(RmgMovesError):
    """Atom ids are missing, duplicated, or were not preserved by RMG."""


# ---------------------------------------------------------------------------
# rows <-> molecules
# ---------------------------------------------------------------------------
# A row is (element, isotope or None, charge, radical electrons, lone pairs,
# ((neighbour id, bond order), ...)) -- the same row shape a MeltState node
# has, so a MeltState component converts directly.

def _isotope(atom):
    iso = atom.element.isotope
    return None if iso == -1 else int(iso)


def molecule_rows(molecule) -> Dict[int, Tuple]:
    """atom id -> row for an RMG Molecule whose atoms carry ids.  Refused
    (ValueError) when two atoms share an id: an unnumbered molecule (every
    id -1) would otherwise collapse onto one row."""
    ids = [atom.id for atom in molecule.atoms]
    if len(set(ids)) != len(ids):
        repeated = sorted({i for i in ids if ids.count(i) > 1})
        raise ValueError(f"molecule_rows needs distinct atom ids; atom id {repeated[:4]} "
                         f"is repeated (an unnumbered molecule carries id -1 on every atom)")
    rows = {}
    for atom in molecule.atoms:
        edges = tuple(sorted((other.id, float(bond.order)) for other, bond in atom.bonds.items()))
        rows[atom.id] = (atom.element.symbol, _isotope(atom), int(atom.charge),
                         int(atom.radical_electrons), int(atom.lone_pairs), edges)
    return rows


def molecule_from_rows(rows: Dict[int, Tuple]) -> Molecule:
    """An RMG Molecule (atom ids = row keys) from one connected set of rows."""
    mol = Molecule()
    atoms = {}
    for i in sorted(rows):
        element, isotope, charge, radical, lone_pairs, _edges = rows[i]
        el = get_element(element) if isotope is None else get_element(element, isotope)
        atom = Atom(element=el, radical_electrons=radical, charge=charge, lone_pairs=lone_pairs)
        atom.id = i
        atoms[i] = atom
        mol.add_atom(atom)
    for i in sorted(rows):
        for j, order in rows[i][5]:
            if j not in rows:
                raise RmgMovesError(f"atom {i} is bonded to {j}, which is not in the rows")
            if i < j:
                mol.add_bond(Bond(atoms[i], atoms[j], order=float(order)))
    mol.update_multiplicity()
    mol.update(sort_atoms=False)
    return mol


def split_rows(rows: Dict[int, Tuple]) -> Tuple[Dict[int, Tuple], ...]:
    """Connected components of a row mapping, ordered by least atom id."""
    seen, parts = set(), []
    for start in sorted(rows):
        if start in seen:
            continue
        stack, members = [start], []
        seen.add(start)
        while stack:
            a = stack.pop()
            members.append(a)
            for b, _o in rows[a][5]:
                if b not in seen:
                    seen.add(b)
                    stack.append(b)
        parts.append({a: rows[a] for a in sorted(members)})
    return tuple(parts)


def molecules_from_state(state, components=None) -> Tuple[Molecule, ...]:
    """RMG Molecules (atom ids = MeltState atom ids) of the given (default:
    all) components of a MeltState, ordered by component id."""
    cids = state.components.ids() if components is None else components
    out = []
    for cid in cids:
        rows = {}
        for a in state.components.atoms(cid):
            n = state.atoms[a]
            rows[a] = (n.element, n.isotope, n.charge, n.radical, n.lone_pairs,
                       tuple(sorted((b, float(o)) for b, o in n.edges.items())))
        out.append(molecule_from_rows(rows))
    return tuple(out)


def _signature(rows: Dict[int, Tuple]) -> Tuple:
    return tuple((i,) + tuple(rows[i][:5]) + (rows[i][5],) for i in sorted(rows))


def _ids(molecule) -> Tuple[int, ...]:
    return tuple(sorted(a.id for a in molecule.atoms))


def _check_ids(reactants):
    seen = set()
    for k, mol in enumerate(reactants):
        if not isinstance(mol, Molecule):
            raise RmgMovesError(f"reactant {k} is a {type(mol).__name__}, not a Molecule")
        for atom in mol.atoms:
            if isinstance(atom.id, bool) or not isinstance(atom.id, int) or atom.id < 0:
                raise AtomIdError(f"reactant {k} has an atom without a non-negative integer id "
                                  f"({atom.id!r})")
            if atom.id in seen:
                raise AtomIdError(f"atom id {atom.id} occurs twice across the reactants")
            seen.add(atom.id)
    if len({id(m) for m in reactants}) != len(reactants):
        raise RmgMovesError("reactants must be different Molecule objects")


# ---------------------------------------------------------------------------
# reactant resonance forms
# ---------------------------------------------------------------------------
def reactant_forms(molecule, resonance: str = "all") -> Tuple[Molecule, ...]:
    """The reactant forms the generator enumerates, ordered by signature.

    ``"all"``: the reactive forms of ``generate_resonance_structures(
    keep_isomorphic=True)`` (all of them if none is reactive); ``"given"``: a
    copy of the input only.  Each form keeps the input atom ids (checked)."""
    if resonance not in RESONANCE_POLICIES:
        raise RmgMovesError(f"resonance must be one of {RESONANCE_POLICIES}, got {resonance!r}")
    base = molecule.copy(deep=True)
    for atom, original in zip(base.atoms, molecule.atoms):
        atom.id = original.id
    if resonance == "given" or len(base.atoms) == 1:
        forms = [base]
    else:
        forms = generate_resonance_structures(base, keep_isomorphic=True, filter_structures=True)
        forms = [f for f in forms if f.reactive] or list(forms)
    ids = _ids(molecule)
    keyed = []
    for form in forms:
        if _ids(form) != ids:
            raise AtomIdError("generate_resonance_structures did not preserve the atom ids")
        keyed.append((_signature(molecule_rows(form)), form))
    keyed.sort(key=lambda kv: kv[0])
    for (s1, _), (s2, _) in zip(keyed, keyed[1:]):
        if s1 == s2:
            raise RmgMovesError("RMG returned the same mapped resonance form twice")
    return tuple(form for _sig, form in keyed)


# ---------------------------------------------------------------------------
# moves
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Move:
    """One RMG template match (one mapping), with its outcome.

    ``order[s]`` is the input reactant matched to template slot ``s``;
    ``forms[k]`` the resonance-form index of input reactant ``k``;
    ``labels`` the reacting atoms ``(slot, label, atom id)``; ``mapped`` every
    mapped atom ``(slot, template structure index, template atom index, atom
    id)``; ``template`` the matched tree node labels (kept and
    ``unchanged_species`` moves; RMG never reads them for the latter);
    ``products`` the mapped product graph ``((atom id, row), ...)`` in the
    reactant atom-id namespace (None unless the recipe ran; a move forbidden
    on its reactant still carries the product the recipe makes, for
    comparison only);
    ``product_labels`` the product labels ``(label, atom id)`` after RMG's
    relabelling."""

    family: str
    direction: str
    order: Tuple[int, ...]
    forms: Tuple[int, ...]
    form_reactive: Tuple[bool, ...]
    labels: Tuple[Tuple[int, str, int], ...]
    mapped: Tuple[Tuple[int, int, int, int], ...]
    status: str
    reason: str
    template: Optional[Tuple[str, ...]]
    products: Optional[Tuple[Tuple[int, Tuple], ...]]
    product_labels: Optional[Tuple[Tuple[str, int], ...]]

    @property
    def pathway(self):
        """Pathway identity: (family, template node labels)."""
        return (self.family, self.template)

    def label_map(self) -> Dict[Tuple[int, str], int]:
        """(slot, label) -> atom id of the reacting atoms."""
        return {(s, label): a for s, label, a in self.labels}

    def atom(self, label: str, slot: Optional[int] = None) -> int:
        hits = [a for s, lab, a in self.labels if lab == label and (slot is None or s == slot)]
        if len(hits) != 1:
            raise KeyError(f"label {label!r} (slot {slot}) maps to {hits}")
        return hits[0]

    def product_rows(self) -> Dict[int, Tuple]:
        if self.products is None:
            raise RmgMovesError(f"a {self.status} move has no product")
        return dict(self.products)

    def sort_key(self):
        return (self.order, self.forms, self.labels, self.mapped)


@dataclass(frozen=True)
class Candidate:
    """What a caller's ``select`` sees before the recipe runs."""

    family: str
    direction: str
    order: Tuple[int, ...]
    forms: Tuple[int, ...]
    labels: Tuple[Tuple[int, str, int], ...]
    molecules: Tuple[Molecule, ...]  # the forms used, by input reactant index

    def atom(self, atom_id):
        for mol in self.molecules:
            for atom in mol.atoms:
                if atom.id == atom_id:
                    return atom
        raise KeyError(atom_id)


def _template_reactants(family, forward, n_reactants):
    """Template reactant groups exactly as ``_generate_reactions`` selects them."""
    if forward:
        template, reactant_num = family.forward_template, family.reactant_num
    else:
        template, reactant_num = family.reverse_template, family.product_num
    if family.auto_generated and reactant_num != n_reactants:
        return None
    if n_reactants > len(template.reactants):
        try:
            return list(template.reactants[0].item.split())
        except AttributeError:
            return [x.item for x in template.reactants]
    return [x.item for x in template.reactants]


def _structures(family, template_reactant):
    """The concrete groups ``_match_reactant_to_template`` matches against."""
    from rmgpy.data.base import Entry, LogicNode
    struct = template_reactant[0] if isinstance(template_reactant, list) else template_reactant
    if isinstance(struct, Entry):
        struct = struct.item
    if isinstance(struct, LogicNode):
        return list(struct.get_possible_structures(family.groups.entries))
    return [struct]


def _mapping_record(slot, mapping, structures):
    out = []
    for atom, tatom in mapping.items():
        for si, group in enumerate(structures):
            if any(tatom is g for g in group.atoms):
                ti = next(k for k, g in enumerate(group.atoms) if g is tatom)
                out.append((slot, si, ti, atom.id, tatom.label))
                break
        else:
            raise RmgMovesError("a template atom is not in any template structure")
    return out


def _forbidden_reason(family, molecule, where):
    """Which forbidden entry (diagnostic only; the verdict is RMG's)."""
    names = []
    forbidden = family.forbidden
    if forbidden is not None:
        from rmgpy.molecule.group import Group
        labelled = molecule.get_all_labeled_atoms()
        for name, entry in sorted(forbidden.entries.items()):
            item = entry.item
            try:
                if isinstance(item, Group):
                    if all(lab in labelled for lab in item.get_all_labeled_atoms()) and \
                            molecule.is_subgraph_isomorphic(item, generate_initial_map=True):
                        names.append(name)
                elif item.is_isomorphic(molecule):
                    names.append(name)
            except Exception:  # diagnostics never change the verdict
                continue
    return f"forbidden {where} ({', '.join(names) or 'family rule'})"


def generate_moves(family, reactants: Sequence[Molecule], direction: str = "forward", *,
                   resonance: str = "all",
                   select: Optional[Callable[[Candidate], bool]] = None) -> Tuple[Move, ...]:
    """Every directed, atom-mapped move of ``family`` on ``reactants``.

    ``family`` is an RMG ``KineticsFamily``; ``direction`` is ``"forward"``
    (the family's template) or ``"reverse"`` (its reverse template; refused
    for own-reverse families, whose forward template already covers both
    directions).  ``select``, if given, is consulted per mapping before the
    recipe runs; refused mappings are returned with status ``excluded``."""
    if direction not in DIRECTIONS:
        raise RmgMovesError(f"direction must be one of {DIRECTIONS}, got {direction!r}")
    forward = direction == "forward"
    if not forward and family.reverse_template is None:
        raise RmgMovesError(f"{family.label} is its own reverse: use direction='forward'")
    reactants = tuple(reactants)
    _check_ids(reactants)
    template_reactants = _template_reactants(family, forward, len(reactants))
    if template_reactants is None or len(template_reactants) != len(reactants) \
            or len(reactants) not in (1, 2):
        return ()
    structures = [_structures(family, t) for t in template_reactants]
    forms = [reactant_forms(m, resonance) for m in reactants]
    moves = []
    if len(reactants) == 1:
        for fa, mol in enumerate(forms[0]):
            for mapping in family._match_reactant_to_template(mol, template_reactants[0]):
                moves.append(_one_move(family, forward, (0,), (fa,), (mol,), [mol], [mapping],
                                       structures, select))
    else:
        for fa, mol_a in enumerate(forms[0]):
            for fb, mol_b in enumerate(forms[1]):
                # A in slot 0, B in slot 1 (RMG passes the structures as [B, A])
                maps_a = family._match_reactant_to_template(mol_a, template_reactants[0])
                maps_b = family._match_reactant_to_template(mol_b, template_reactants[1])
                for map_a in maps_a:
                    for map_b in maps_b:
                        moves.append(_one_move(family, forward, (0, 1), (fa, fb), (mol_a, mol_b),
                                               [mol_b, mol_a], [map_b, map_a], structures, select,
                                               slots=(1, 0)))
                # B in slot 0, A in slot 1 (the reactants are always different objects)
                maps_a = family._match_reactant_to_template(mol_a, template_reactants[1])
                maps_b = family._match_reactant_to_template(mol_b, template_reactants[0])
                for map_a in maps_a:
                    for map_b in maps_b:
                        moves.append(_one_move(family, forward, (1, 0), (fa, fb), (mol_a, mol_b),
                                               [mol_a, mol_b], [map_a, map_b], structures, select,
                                               slots=(1, 0)))
    moves.sort(key=Move.sort_key)
    for a, b in zip(moves, moves[1:]):
        if a.sort_key() == b.sort_key():
            raise RmgMovesError("two mappings share one key")
    return tuple(moves)


def _one_move(family, forward, order, form_idx, molecules, structs, maps, structures, select,
              slots=(0,)):
    """Run one mapping through RMG's per-mapping steps and classify it.

    ``structs``/``maps`` are in the order RMG passes them to the recipe;
    ``slots[k]`` is the template slot of ``structs[k]``."""
    label = family.label
    direction = "forward" if forward else "reverse"
    records = []
    for k, mapping in enumerate(maps):
        records.extend(_mapping_record(slots[k] if len(maps) > 1 else 0, mapping,
                                       structures[slots[k] if len(maps) > 1 else 0]))
    labels = tuple(sorted((s, lab, a) for s, _si, _ti, a, lab in records if lab))
    mapped = tuple(sorted((s, si, ti, a) for s, si, ti, a, _lab in records))
    reactive = tuple(bool(m.reactive) for m in molecules)

    def done(status, reason, template=None, products=None, product_labels=None):
        return Move(label, direction, order, form_idx, reactive, labels, mapped, status, reason,
                    template, products, product_labels)

    if select is not None:
        candidate = Candidate(label, direction, order, form_idx, labels, molecules)
        if not select(candidate):
            return done("excluded", "select() refused the mapping; recipe not applied")

    # -- as _generate_product_structures: label, reactant forbidden check
    for struct in structs:
        struct.clear_labeled_atoms()
    for mapping in maps:
        for atom, tatom in mapping.items():
            atom.label = tatom.label
    reactant_forbidden = None
    for struct in structs:
        if family.is_molecule_forbidden(struct):
            reactant_forbidden = _forbidden_reason(family, struct, "reactant")
            break
    # -- recipe (unique=False; RMG never reads it).  For a forbidden reactant
    # RMG stops here; the recipe still runs so the move carries the product it
    # would have made (for comparison only; the status stays ``forbidden``).
    failure = None
    try:
        products = family.apply_recipe(structs, forward=forward, unique=False, relabel_atoms=True)
    except (InvalidActionError, KekulizationError) as exc:
        failure = f"{type(exc).__name__}: {exc}"
    except ActionError as exc:
        failure = f"ActionError: {exc}"
    else:
        if not products:
            failure = ("apply_recipe returned no products (product count, net charge or "
                       "product multiplicity check)")
    if failure is not None:
        if reactant_forbidden:
            return done("forbidden", f"{reactant_forbidden}; recipe also failed: {failure}")
        return done("recipe_failed", failure)
    rows = {}
    for p in products:
        rows.update(molecule_rows(p))
    reactant_ids = sorted(a.id for m in molecules for a in m.atoms)
    if sorted(rows) != reactant_ids or len(rows) != sum(len(p.atoms) for p in products):
        raise AtomIdError(f"{label}: apply_recipe did not preserve the reactant atom ids")
    product_rows = tuple(sorted(rows.items()))
    product_labels = tuple(sorted((a.label, a.id) for p in products for a in p.atoms if a.label))
    if reactant_forbidden:
        return done("forbidden", reactant_forbidden, products=product_rows,
                    product_labels=product_labels)
    for p in products:
        if family.is_molecule_forbidden(p):
            return done("forbidden", _forbidden_reason(family, p, "product"),
                        products=product_rows, product_labels=product_labels)
    # species / model-growth constraints are OFF (never consulted)
    # -- as _create_reaction
    if same_species_lists(list(structs), products, save_order=family.save_order):
        # RMG drops the move here (its _create_reaction returns None).  The
        # template labels are read first, on the reaction _create_reaction
        # would have built, so a caller can rate the move by its template.
        template, why = _unchanged_species_template(family, structs, products, forward)
        return done("unchanged_species", "products are the reactants" + why, template=template,
                    products=product_rows, product_labels=product_labels)
    rxn = family._create_reaction(list(structs), products, forward)
    if rxn is None:
        if not family.allow_charged_species and any(
                s.get_net_charge() != 0 for s in list(structs) + list(products)):
            return done("charged_species", "a charged species in a family that forbids them",
                        products=product_rows, product_labels=product_labels)
        return done("unbalanced", "_create_reaction refused the reaction (not balanced)",
                    products=product_rows, product_labels=product_labels)
    # -- as the end of _generate_reactions: restore labels, read the template
    for reactant in rxn.reactants:
        reactant.clear_labeled_atoms()
    for lab, atom in rxn.labeled_atoms["reactants"].items():
        if isinstance(atom, list):
            for a in atom:
                a.label = lab
        else:
            atom.label = lab
    template = tuple(family.get_reaction_template_labels(rxn))
    for struct in structs:
        struct.clear_labeled_atoms()
    return done("kept", "", template=template, products=product_rows,
                product_labels=product_labels)


def _unchanged_species_template(family, structs, products, forward):
    """Template labels of a move whose products are its reactants: the
    reaction ``_create_reaction`` builds, minus its same-species refusal, and
    the label restore + ``get_reaction_template_labels`` of the kept path.
    Returns (labels or None, reason suffix)."""
    from rmgpy.data.kinetics.family import TemplateReaction
    from rmgpy.exceptions import UndeterminableKineticsError
    rxn = TemplateReaction(reactants=list(structs) if forward else list(products),
                           products=list(products) if forward else list(structs),
                           degeneracy=1, reversible=family.reversible, family=family.label,
                           is_forward=forward, electrons=family.electrons)
    for key, species_list in zip(["reactants", "products"], [rxn.reactants, rxn.products]):
        for species in species_list:
            rxn.labeled_atoms[key] = dict(rxn.labeled_atoms[key], **species.get_all_labeled_atoms())
    for reactant in rxn.reactants:
        reactant.clear_labeled_atoms()
    for lab, atom in rxn.labeled_atoms["reactants"].items():
        if isinstance(atom, list):
            for a in atom:
                a.label = lab
        else:
            atom.label = lab
    try:
        template = tuple(family.get_reaction_template_labels(rxn))
    except UndeterminableKineticsError as exc:
        return None, f"; template undeterminable: {exc}"
    finally:
        for struct in structs:
            struct.clear_labeled_atoms()
    return template, ""


def status_counts(moves: Iterable[Move]) -> Dict[str, int]:
    counts = Counter(m.status for m in moves)
    return {s: counts.get(s, 0) for s in STATUSES}


# ---------------------------------------------------------------------------
# bridge: product identity modulo resonance (RMG resonance + isomorphism)
# ---------------------------------------------------------------------------
DEFAULT_KEY_CACHE_ENTRIES = 4096


def _aromatic_bonds(molecule) -> int:
    return sum(1 for bond in molecule.get_all_edges() if bond.is_benzene())


class ProductKeys:
    """Product identity modulo resonance, with a bounded cache.

    A product component's key is the least integer-keyed signature over the
    reactive forms ``generate_resonance_structures(keep_isomorphic=True)``
    returns for it (all forms if none is reactive), generated from its most
    aromatic form (see ``component``).  Every form generated is cached under
    the same key; ``cache=False`` recomputes every time.  The cache holds at most
    ``max_entries`` signatures, least recently used evicted first (a miss only
    recomputes, so eviction never changes a key)."""

    def __init__(self, cache: bool = True, max_entries: int = DEFAULT_KEY_CACHE_ENTRIES):
        if isinstance(max_entries, bool) or not isinstance(max_entries, int) or max_entries < 1:
            raise RmgMovesError(f"max_entries must be a positive integer, got {max_entries!r}")
        self.cache = cache
        self.max_entries = max_entries
        self._known = OrderedDict()
        self.computed = 0

    def _remember(self, sig, key):
        known = self._known
        if sig in known:
            known.move_to_end(sig)
        known[sig] = key
        while len(known) > self.max_entries:
            known.popitem(last=False)

    def forms(self, rows: Dict[int, Tuple]) -> Tuple[Molecule, ...]:
        mol = molecule_from_rows(rows)
        if len(mol.atoms) == 1:
            return (mol,)
        ids = _ids(mol)
        forms = tuple(generate_resonance_structures(mol, keep_isomorphic=True,
                                                    filter_structures=True))
        for f in forms:
            if _ids(f) != ids:
                raise AtomIdError("generate_resonance_structures did not preserve the atom ids")
        return forms

    def component(self, rows: Dict[int, Tuple]) -> Tuple:
        sig = _signature(rows)
        if self.cache and sig in self._known:
            self._known.move_to_end(sig)
            return self._known[sig]
        self.computed += 1
        forms = self.forms(rows)
        sigs = [(_signature(molecule_rows(f)), bool(f.reactive)) for f in forms]
        # Re-root at the most aromatic form.  RMG can return a Kekule input as
        # a REACTIVE form beside the aromatic one (seen on PS chains: the ring
        # of a recombined benzylic radical, written in either Kekule pattern),
        # so the set generated from the input depends on the input; the set
        # generated from the most aromatic form (least signature among those
        # with the most 1.5-order bonds) does not.
        aromatic = max(_aromatic_bonds(f) for f in forms)
        root = min(s for (s, _r), f in zip(sigs, forms) if _aromatic_bonds(f) == aromatic)
        seen = list(sigs)
        if root != sig:
            rooted = self.forms(dict((r[0], r[1:]) for r in root))
            sigs = [(_signature(molecule_rows(f)), bool(f.reactive)) for f in rooted]
            seen += sigs
        reactive = [s for s, r in sigs if r] or [s for s, _r in sigs]
        key = min(reactive)
        if self.cache:
            for s, _r in seen:
                if s != sig:
                    self._remember(s, key)
            self._remember(sig, key)
        return key

    def __call__(self, item) -> Tuple:
        """Key of a Move (its product) or of an atom-id -> row mapping."""
        rows = item.product_rows() if isinstance(item, Move) else item
        return tuple(sorted(self.component(part) for part in split_rows(rows)))


_DEFAULT_KEYS = ProductKeys()


def product_key(item) -> Tuple:
    """The matcher's comparison form of a move (or of product rows): the full
    mapped product of the touched components modulo resonance."""
    return _DEFAULT_KEYS(item)


def same_product(rows_a: Dict[int, Tuple], rows_b: Dict[int, Tuple]) -> bool:
    """Atom-id-preserving equality modulo resonance by RMG isomorphism: per
    component, some RMG form of ``a`` is isomorphic to some RMG form of ``b``
    under the identity atom-id map (``strict=True``: electrons, bond orders
    and charges must agree).  With a complete initial map RMG also validates
    the map strictly (``is_mapping_valid``), so electrons are compared even
    before the isomorphism search."""
    parts_a, parts_b = split_rows(rows_a), split_rows(rows_b)
    if [sorted(p) for p in parts_a] != [sorted(p) for p in parts_b]:
        return False
    for pa, pb in zip(parts_a, parts_b):
        forms_b = _DEFAULT_KEYS.forms(pb)
        hit = False
        for form_a in _DEFAULT_KEYS.forms(pa):
            by_id = {a.id: a for a in form_a.atoms}
            for form_b in forms_b:
                initial = {by_id[b.id]: b for b in form_b.atoms}
                if form_a.is_isomorphic(form_b, initial_map=initial, strict=True):
                    hit = True
                    break
            if hit:
                break
        if not hit:
            return False
    return True


# ---------------------------------------------------------------------------
# comparison only: RMG's final degeneracy
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RmgReaction:
    """One reaction of RMG's full pipeline, with its final degeneracy."""

    family: str
    template: Tuple[str, ...]
    is_forward: bool
    degeneracy: float
    products: Tuple[str, ...]  # SMILES, sorted
    product_molecules: Tuple[Molecule, ...]


def rmg_degeneracy(kinetics_database, reactants: Sequence[Molecule], family: str,
                   direction: str = "forward") -> Tuple[RmgReaction, ...]:
    """RMG's final reactions and degeneracies for ``reactants`` in ``family``
    (``generate_reactions_from_families``: resonance, both template orders,
    ``find_degenerate_reactions``, the A + A one-half, the forward-direction
    recompute for reactions found backwards).  For comparison only: moves are
    never derived from this.  ``direction`` keeps the reactions whose
    ``is_forward`` matches."""
    from rmgpy.species import Species
    species = []
    for mol in reactants:
        spc = Species(molecule=[mol.copy(deep=True)])
        spc.generate_resonance_structures(keep_isomorphic=False)
        species.append(spc)
    rxns = kinetics_database.generate_reactions_from_families(species, only_families=[family])
    out = []
    want_forward = direction == "forward"
    for rxn in rxns:
        if bool(rxn.is_forward) != want_forward:
            continue
        prods = rxn.products if rxn.is_forward else rxn.reactants
        mols = tuple(p.molecule[0] if hasattr(p, "molecule") else p for p in prods)
        out.append(RmgReaction(rxn.family, tuple(rxn.template or ()), bool(rxn.is_forward),
                               float(rxn.degeneracy), tuple(sorted(m.to_smiles() for m in mols)),
                               mols))
    out.sort(key=lambda r: (r.template, r.products, r.degeneracy))
    return tuple(out)


def degeneracy_table(moves: Iterable[Move], rmg_reactions: Sequence[RmgReaction]):
    """Group kept moves by (family, template labels, product species) and set
    RMG's final degeneracy beside the directed move count and the count of
    distinct mapped products.  Rows: dict per group (reported, not asserted)."""
    from rmgpy.species import Species
    groups = defaultdict(list)
    species_cache = []

    def species_index(rows):
        mol = molecule_from_rows(rows)
        for k, spc in enumerate(species_cache):
            if spc.is_isomorphic(mol):
                return k
        spc = Species(molecule=[mol.copy(deep=True)])
        spc.generate_resonance_structures(keep_isomorphic=False)
        species_cache.append(spc)
        return len(species_cache) - 1

    for m in moves:
        if m.status != "kept":
            continue
        prods = tuple(sorted(species_index(p) for p in split_rows(m.product_rows())))
        groups[(m.family, tuple(sorted(m.template)), prods)].append(m)
    table = []
    for (fam, template, prods), members in sorted(groups.items()):
        smiles = tuple(sorted(species_cache[k].molecule[0].to_smiles() for k in prods))
        match = [r for r in rmg_reactions
                 if tuple(sorted(r.template)) == template and len(r.products) == len(prods)
                 and sorted(species_index(molecule_rows(m)) for m in r.product_molecules)
                 == sorted(prods)]
        table.append({
            "family": fam, "template": template, "products": smiles,
            "directed_moves": len(members),
            "distinct_mapped": len({product_key(m) for m in members}),
            "rmg_degeneracy": [r.degeneracy for r in match],
        })
    return table
