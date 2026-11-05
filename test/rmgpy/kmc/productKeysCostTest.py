"""Cost of ``rmg_moves.ProductKeys`` on multi-ring aromatic products (I-305).

A product key is the least reactive signature over RMG's resonance forms of
the product's most aromatic form.  Its cost is RMG's resonance generation.
These tests pin what is saved without changing a single key:

- a root's forms are generated once per ProductKeys, however many inputs
  (Kekule variants of spectator rings, patterns of the reacting ring) share it;
- a cache entry is stored compactly, so the default bound holds the keys of
  many large channel sets in less memory than the old 4096 full signatures;
- every key equals the uncached key of a frozen copy of the pre-I-305
  algorithm (``reference_key``) on multi-ring products, radicals and
  diradicals, Kekule variants and atom-id permutations, fed in random order
  through one shared cache.
"""

import itertools
import random
import sys

import pytest

from rmgpy.kmc import rmg_moves as rm
from rmgpy.molecule import Molecule
from rmgpy.molecule.resonance import generate_resonance_structures


# --------------------------------------------------------------------------
# the pre-I-305 key, frozen (ProductKeys.component, cache=False, at c7741afef)
# --------------------------------------------------------------------------
def _reference_forms(rows):
    mol = rm.molecule_from_rows(rows)
    if len(mol.atoms) == 1:
        return (mol,)
    return tuple(generate_resonance_structures(mol, keep_isomorphic=True, filter_structures=True))


def _reference_component(rows):
    sig = rm._signature(rows)
    forms = _reference_forms(rows)
    sigs = [(rm._signature(rm.molecule_rows(f)), bool(f.reactive)) for f in forms]
    aromatic = max(rm._aromatic_bonds(f) for f in forms)
    root = min(s for (s, _r), f in zip(sigs, forms) if rm._aromatic_bonds(f) == aromatic)
    if root != sig:
        rooted = _reference_forms(dict((r[0], r[1:]) for r in root))
        sigs = [(rm._signature(rm.molecule_rows(f)), bool(f.reactive)) for f in rooted]
    reactive = [s for s, r in sigs if r] or [s for s, _r in sigs]
    return min(reactive)


def reference_key(rows):
    return tuple(sorted(_reference_component(part) for part in rm.split_rows(rows)))


# --------------------------------------------------------------------------
# products
# --------------------------------------------------------------------------
UNIT = "CC(c1ccccc1)"
ENDS = {"radical": "C[CH]c1ccccc1", "alkene": "C=Cc1ccccc1", "ring": "CC1C=CC=C[CH]1"}


def ps_smiles(units, end="radical", head="C"):
    """A polystyrene-like segment with ``units`` phenyl rings (``head``
    "[CH2]" with end "radical" makes a diradical)."""
    return head + "C(c1ccccc1)" + UNIT * max(units - 2, 0) + ENDS[end]


def rows_of(smiles, ids=None):
    m = Molecule().from_smiles(smiles)
    ids = list(range(len(m.atoms))) if ids is None else list(ids)
    for atom, i in zip(m.atoms, ids):
        atom.id = i
    return rm.molecule_rows(m)


def aromatic_rows(rows):
    """The rows of the RMG form with the most benzene bonds."""
    return rm.molecule_rows(max(_reference_forms(rows), key=rm._aromatic_bonds))


def benzene_rings(rows):
    """The 6-cycles of 1.5-order bonds, each in ring order."""
    adj = {i: sorted(j for j, o in r[5] if o == 1.5) for i, r in rows.items()}
    rings, used = [], set()
    for start in sorted(adj):
        if start in used or len(adj[start]) != 2:
            continue
        cycle, prev, cur, closed = [start], None, start, False
        for _step in range(6):
            nxt = [j for j in adj[cur] if j != prev]
            if len(nxt) != 1 and cur != start:
                break                      # a fused-ring junction: not walked
            prev, cur = cur, nxt[0]
            if cur == start:
                closed = True
                break
            cycle.append(cur)
        if closed and len(cycle) == 6:
            rings.append(cycle)
            used.update(cycle)
    return rings


def kekule(rows, ring, phase):
    """``rows`` with ``ring`` written in one of its two Kekule patterns."""
    orders = {frozenset((ring[k], ring[(k + 1) % 6])): (2.0 if (k + phase) % 2 == 0 else 1.0)
              for k in range(6)}
    return {i: r[:5] + (tuple(sorted((j, orders.get(frozenset((i, j)), o)) for j, o in r[5])),)
            for i, r in rows.items()}


def kekule_variants(rows, limit=None):
    """Every combination of (aromatic, pattern 0, pattern 1) over the rings."""
    rings = benzene_rings(rows)
    out = []
    for choice in itertools.product((None, 0, 1), repeat=len(rings)):
        r = rows
        for ring, phase in zip(rings, choice):
            if phase is not None:
                r = kekule(r, ring, phase)
        out.append(r)
    return out if limit is None else out[:limit]


def relabel(rows, mapping):
    return {mapping[i]: r[:5] + (tuple(sorted((mapping[j], o) for j, o in r[5])),)
            for i, r in rows.items()}


def _deep_size(obj, seen):
    if id(obj) in seen:
        return 0
    seen.add(id(obj))
    size = sys.getsizeof(obj)
    if isinstance(obj, dict):
        size += sum(_deep_size(k, seen) + _deep_size(v, seen) for k, v in obj.items())
    elif isinstance(obj, (tuple, list)):
        size += sum(_deep_size(x, seen) for x in obj)
    return size


@pytest.fixture
def counted(monkeypatch):
    """Count the RMG resonance generations rmg_moves makes."""
    calls = []
    real = rm.generate_resonance_structures

    def counting(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)

    monkeypatch.setattr(rm, "generate_resonance_structures", counting)
    return calls


# --------------------------------------------------------------------------
# cost
# --------------------------------------------------------------------------
def test_a_shared_root_is_generated_once(counted):
    """The nine Kekule variants of a two-ring benzylic radical share one
    aromatic root.  A variant that misses the cache costs one generation (its
    own forms, to find the root); the root's forms are generated once."""
    base = aromatic_rows(rows_of("CC(c1ccccc1)C[CH]c1ccccc1"))
    variants = kekule_variants(base)
    assert len({rm._signature(v) for v in variants}) == 9
    want = {reference_key(v) for v in variants}
    assert len(want) == 1
    del counted[:]
    keys = rm.ProductKeys()
    got = {keys(v) for v in variants}
    assert got == want
    assert keys.computed >= 6          # the variants really miss the cache
    assert len(counted) == keys.computed


def test_a_cache_entry_is_compact():
    """A 161-atom product's cache entries cost under 4 KB each (a full
    integer-keyed signature costs ~58 KB), so the default bound of such keys
    stays below the memory of the old 4096 full signatures."""
    rows = rows_of(ps_smiles(10))
    assert len(rows) == 161
    keys = rm.ProductKeys()
    keys(rows)
    entries = len(keys._known)
    assert entries >= 2
    per_entry = _deep_size(keys._known, set()) / entries
    assert per_entry < 4096
    assert rm.DEFAULT_KEY_CACHE_ENTRIES >= 65536
    assert per_entry * rm.DEFAULT_KEY_CACHE_ENTRIES < 4096 * 58000


def test_a_hit_returns_the_key_as_computed():
    rows = rows_of(ps_smiles(3))
    keys = rm.ProductKeys()
    first = keys(rows)
    second = keys(rows)
    assert keys.computed == 1
    assert first == second == reference_key(rows)
    assert repr(first) == repr(second)


def test_int_and_float_rows_share_an_entry():
    """Rows whose bond orders arrive as ints hit the entry their float twin
    made, as the tuple-keyed cache did (1 == 1.0)."""
    rows = rows_of("CC[CH]c1ccccc1")
    as_int = {i: r[:5] + (tuple((j, int(o)) if o == int(o) else (j, o) for j, o in r[5]),)
              for i, r in rows.items()}
    assert any(type(o) is int for r in as_int.values() for _j, o in r[5])
    keys = rm.ProductKeys()
    a = keys(rows)
    b = keys(as_int)
    assert keys.computed == 1 and a == b


def test_the_root_memo_is_bounded():
    keys = rm.ProductKeys(max_entries=3)
    for n in range(1, 8):
        keys(rows_of("C" * n + "[CH]c1ccccc1"))
        assert len(keys._known) <= 3 and len(keys._rooted) <= 3


# --------------------------------------------------------------------------
# identity: randomized differential against the frozen pre-I-305 key
# --------------------------------------------------------------------------
def _differential_inputs(rng):
    inputs = []
    for smiles in [ps_smiles(2), ps_smiles(3), ps_smiles(3, "alkene"), ps_smiles(3, "ring"),
                   ps_smiles(3, head="[CH2]"), "C[CH]c1ccccc1", "[CH2]c1ccc(C[CH]c2ccccc2)cc1",
                   "C=CC(c1ccccc1)C[C](C)c1ccccc1", "c1ccc2ccccc2c1C[CH]C", "CCCC"]:
        base = rows_of(smiles)
        variants = kekule_variants(aromatic_rows(base), limit=6) + [base]
        n = len(base)
        for v in variants:
            inputs.append(v)
            # an order-preserving shift, and a random permutation of the atom ids
            inputs.append(relabel(v, {i: i + 1000 for i in range(n)}))
            perm = list(range(n))
            rng.shuffle(perm)
            inputs.append(relabel(v, dict(zip(range(n), perm))))
    # two components in one rows mapping
    a, b = rows_of("C[CH]c1ccccc1"), rows_of("CC(c1ccccc1)C=C", ids=range(500, 600))
    inputs.append({**a, **b})
    return inputs


def test_keys_equal_the_frozen_uncached_key_in_random_order():
    rng = random.Random(305)
    inputs = _differential_inputs(rng)
    assert len(inputs) > 150
    want = [reference_key(r) for r in inputs]
    assert len(set(want)) > 30
    for trial, bound in enumerate([8, 64, rm.DEFAULT_KEY_CACHE_ENTRIES]):
        order = list(range(len(inputs)))
        rng.shuffle(order)
        keys = rm.ProductKeys(max_entries=bound)
        for k in order:
            assert keys(inputs[k]) == want[k], (trial, k)
        if bound == rm.DEFAULT_KEY_CACHE_ENTRIES:
            assert keys.computed < len(inputs)   # the shared cache really hits
    fresh = rm.ProductKeys(cache=False)
    for k in order[:30]:
        assert fresh(inputs[k]) == want[k]
