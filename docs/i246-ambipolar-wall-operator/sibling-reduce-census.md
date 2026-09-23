# Sibling `__reduce__` census — proposed ticket **I-247**

Commissioned while fixing `PlasmaReactor.__reduce__` on `i246-ambipolar-wall-operator`. The
question was whether the sibling reactors carry the same hand-enumeration gap for their own
later-added parameters.

**They do, and the gap is pre-existing and wider than the plasma one.** Nothing in this census
was fixed — the branch was not widened. This file exists to name the ticket.

Every line below was measured on this build, not inferred from reading.

## Result

| reactor | `__init__` params | `__reduce__` carries | round trip | severity |
|---|---|---|---|---|
| `PlasmaReactor` | 18 | **18** (fixed on this branch) | preserved by value | — |
| `SimpleReactor` | 9 | 5 | **reconstructs SILENTLY, wrong** | **HIGH** |
| `MBSampledReactor` | 9 | 4 | `TypeError` | MEDIUM |
| `LiquidReactor` | 8 | inherits base's 1 | `TypeError` | MEDIUM |
| `SurfaceReactor` | 13 | inherits base's 1 | `TypeError` | MEDIUM |
| `ReactionSystem` (base) | 3 | 1 | — | LOW |

## The split that matters

The reactors fail in two different ways, and only one of them is dangerous.

**Silent, and therefore the real finding — `SimpleReactor`.** It defines its own `__reduce__`
carrying `(T, P, initial_mole_fractions, n_sims, termination)` and omitting four parameters. The
arity still satisfies `__init__`, so a round trip *succeeds* and quietly substitutes defaults:

```
sensitivity_threshold  1e-09      ->  0.001      (the default)
const_spc_names        ['X']      ->  None
```

This is exactly the plasma defect — a valid object that has silently lost what was omitted — and
`SimpleReactor` is the most-used reactor in the tree, so it is ranked above the others despite the
plasma fix being the one that prompted the census.

**Loud — `LiquidReactor`, `SurfaceReactor`, `MBSampledReactor`.** These cannot reconstruct at all:

- `LiquidReactor` and `SurfaceReactor` define no `__reduce__` and inherit
  `ReactionSystem.__reduce__` (`base.pyx:224`), which returns `(self.__class__, (self.termination,))`.
  Unpickling a `LiquidReactor` therefore calls `LiquidReactor([])` and raises
  `TypeError: __init__() takes at least 2 positional arguments (1 given)`. Note the *shape* of the
  bug is worse than omission: `termination` is being passed into the `T` slot. It is saved from
  silent corruption only by the arity check.
- `MBSampledReactor` carries `(T, P, initial_mole_fractions, termination)` against an `__init__` of
  `(T, P, initial_mole_fractions, k_sampling, constantSpeciesList, termination, ...)`, so
  `termination` is aimed at the **`k_sampling`** slot. It raises
  `TypeError: __init__() takes at least 6 positional arguments (4 given)`. Again a positional
  misalignment that arity happens to catch — adding two more parameters to that call would convert
  it from loud to silent.

## Why this recurs

`__reduce__` here is a hand-written enumeration of a constructor signature with no mechanism tying
the two together. Nothing fails when it falls behind, because the omission is only visible as a
*default value where a set value should be*. The plasma instance went unnoticed through the whole
of I-246 — implementation, verification, a 27-check verifier and a 588-point sweep — and was found
only by an outside reviewer.

The durable fix is not "update five more enumerations"; it is to make falling behind fail. The
guard used on this branch generalises: parse `__init__`'s parameter list and assert each name
appears in `__reduce__`'s body. That check is in
`test/rmgpy/solver/plasmaWallTest.py::test_reduce_enumerates_every_constructor_parameter` and is
written against `plasma.pyx` only; lifting it to a parametrised check over all five solver modules
is the substance of the proposed ticket.

Note for whoever takes it: `inspect.signature` does **not** work here. `PlasmaReactor` and its
siblings are Cython `cdef class`es whose `__init__` introspects as `(*args, **kwargs)`; a first
draft of the guard compared argument counts via `inspect` and was red for that reason rather than
for the defect. The check has to read the `.pyx` source.

## Scope note

`is_pickleable`/round-trip behaviour of a whole reaction system has no production caller in this
tree — nothing pickles or deep-copies one, and the only nearby `deepcopy` copies
`initial_mole_fractions`. So this is a latent break of the public API contract rather than a live
bug, which is why it is a separate ticket and not an I-246 blocker.
