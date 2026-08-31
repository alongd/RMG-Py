# I-198 — close the "found-but-unset `rmgrc` silently loads the default database" hole

`Settings.require_database_directory()` (`rmgpy/__init__.py`) refused two ways for a run to proceed
against a database the user never named — no `rmgrc` found at all, and a `database.directory` naming
a non-existent directory — but left a third, likelier door open. An `rmgrc` that **was** found yet
produced no `database.directory` substring match (the key absent, an INI `[section]` spelling, or a
mis-spelling not containing the substring) left the value at the compiled-in default, which is a real
existing directory on most checkouts, so both prior checks passed and the run proceeded silently
against the wrong database. This ticket adds a third refusal on that exact condition (see the rule
section for what it does and does not catch), in the same voice as the existing two.

Measured at branch `i198-rmgrc-guard` @ `b12ef4c43` (base, zero commits behind `plasma`), `rmg_env`,
this worktree's own `rmgpy` forced ahead of the shared `easy-install.pth` pin (`/home/alon/Code/RMG-Py`)
via `PYTHONPATH`. `/proc/loadavg` was ~5 throughout (moderate; sub-second unit tests unaffected).

## Questions 1–4 (probed, not assumed)

**Reproduction of the three ticket rows** (each run as `settings.require_database_directory()` from a
working dir holding the given `rmgrc`):

| `rmgrc` contents | `filename` | `sources["database.directory"]` | outcome |
|---|---|---|---|
| `database.directory = /this/path/does/not/exist/input` | `rmgrc` | `from rmgrc` | **SettingsError raised** (door 2) |
| `[database]` / `directory = /nonexistent/path/input` | `rmgrc` | `Default, relative to RMG-Py source code` | **returned default** (the hole) |
| `totally.bogus = 1` | `rmgrc` | `Default, relative to RMG-Py source code` | **returned default** (the hole) |

All three reproduced exactly as the ticket stated.

**Q1 — what `reset()` sets, and is `sources` a reliable discriminator?**
`reset()` (`rmgpy/__init__.py:196`–`202`) sets `self.filename = None`, `self["database.directory"]` to
`os.path.realpath(os.path.join(rmgpy_dir, "..", "..", "RMG-database", "input"))`
(= `/home/alon/Code/RMG-database/input` here, which **exists** → the `isdir` check passes), and
`self.sources["database.directory"] = "Default, relative to RMG-Py source code"`. The loader
(`load()`, line ~134) overwrites the source to `"from {filename}"` **only** when a `database.directory`
line matches. So `sources["database.directory"].startswith("from ")` is a true, reliable discriminator
between "set from the file" and "still the compiled-in default". The run banner already prints
`... (Default, relative to RMG-Py source code)` in the hole cases — the evidence was on screen and
nothing acted on it. **Discriminator confirmed reliable.**

**Q2 — is `test_data.directory` exposed to the same hole?** Yes: same parse loop, same substring
match. Nothing guards it — there is no `require_test_data_directory()`, and it is not read on the run
path this ticket concerns (its default exists and feeds tests). **Out of scope**, recorded.

**Q3 — how loose is the substring match?** `line.find("database.directory") != -1` matches any line
*containing* the substring. Measured mis-fires:
- `database.directory_old = /should/not/win` → wrongly sets `database.directory` (source `from rmgrc`, value `/should/not/win`).
- `test_data.directory = /foo/database.directory/bar` (substring in the **value path**) → wrongly sets `database.directory` to `/foo/database.directory/bar`.
- Correct spelling wins normally; a `# database.directory = ...` comment (even leading-space) is correctly ignored (comment stripping runs first).

The two branches are `elif`, so both cannot fire for one line, but a single line can fire the *wrong*
branch relative to intent. **This is a separate weakness and is left out of scope**: both mis-fires set
source to `from rmgrc`, so they are the "wrong value from a real-looking line" class — which this
guard neither catches nor should — not the typo/section hole. Tightening the match would change what
currently-working `rmgrc` files parse to (a guard firing *differently*), which the Verifier does not
require and which the ticket explicitly permits declining.

**Q4 — callers of `require_database_directory()` and direct readers.** Exactly one caller:
`rmgpy/rmg/input.py:99`, the `database()` directive, invoked as the input file is read. This fix is
**complete for the RMG run path and out of scope by design for Arkane and every direct
`settings["database.directory"]` reader** — `arkane/input.py:90`, `arkane/modelchem.py`,
`arkane/encorr/*`, `rmgpy/rmg/model.py:1720/1848`, `rmgpy/data/kinetics/family.py:857/4078`,
`scripts/*`, etc. — which read the setting without consulting the guard and are unchanged (a known,
separate gap). The ticket does not widen the guard to new call sites; it closes the third door on the
one path that already consults it.

## The rule the guard enforces (and what it misses)

The refusal fires **iff no line in the found `rmgrc` matched the loader's `database.directory`
substring scan**, leaving `sources["database.directory"] == Settings.DEFAULT_SOURCE`. That covers an
absent key, an INI `[section]` spelling, and a mis-spelling that does **not** contain the substring
(e.g. `databse.directory`). It does **not** cover `database.directory_old = /x`: that line *matches*
the substring scan, so the loader sets the value and rewrites the source to `from <file>`, and it
sails past this refusal — caught only by the existing directory check, and only if the mis-parsed
value is not an existing directory. That loose-substring-match shape is the weakness left out of
scope (Q3). The earlier draft of this report and the commit message described the guard as covering
typos generally; it does not — the rule above is what it enforces.

## The fix

`rmgpy/__init__.py` — one refusal added inside `require_database_directory()`, between the existing
"no file" (door 1) and "non-existent directory" (door 2) checks. **The discriminator is
`sources["database.directory"] == Settings.DEFAULT_SOURCE` ("the value IS the compiled-in default"),
not `startswith("from ")` ("the value came from a file").** The distinction is load-bearing: a value
assigned programmatically (`settings["database.directory"] = ...`) carries source `"-"`, is not the
default, and must not be refused — keying on provenance would false-refuse a correctly-set value, a
false positive worse than the hole. `DEFAULT_SOURCE` is now a named class constant shared by `reset()`
and the guard so the two cannot drift apart. The message names the file read, states no
`database.directory` line took effect, shows the default it refuses to fall back to, gives the exact
one-line spelling, and points at both templates — matching the shape and tone of the two existing
refusals. `rmgpy/__init__.py` is not cythonized; no rebuild of it is required. (The worktree's
extensions were built once via `make build` so the shared `test/conftest.py` DB-singleton fixture
could import the full stack.)

## Deliberate compatibility break

An `rmgrc` that sets only `test_data.directory` and relied on the compiled-in `database.directory`
default now hard-fails in `require_database_directory()`. **This is intended, not incidental**:
relying on a silently-chosen database default is exactly the failure this guard exists to stop. It is
recorded here and in the method docstring so it is a documented break rather than a discovered one.

## Verification

| # | Check | Result |
|---|---|---|
| 1 | Q1–4 answered with measurements | above |
| 2 | New test RED before / GREEN after, INI-section + typo'd key covered | RED **3 failed, 3 passed** (holes fail) → GREEN **7 passed** (6 original + programmatic-set) |
| 3 | Two existing refusals still refuse; correct spelling resolves; programmatic set NOT refused | `test_missing_file_still_refused`, `test_nonexistent_directory_still_refused`, `test_correctly_spelled_resolves`, `test_programmatic_set_not_refused` (all GREEN) |
| 4 | No regression in the caller's suite | `test/rmgpy/rmg/inputTest.py` **101 passed, 1 skipped** both before (HEAD) and after (patched), one process each — *with an rmgrc present, see trap below* |
| 5 | 5 torr argon deck end-to-end | `PYTHON_EXIT=0`, one `MODEL GENERATION COMPLETED`, core **3 species / 1 reaction** |
| 6 | One suite per process | every suite above run alone |

**Item-1 behavioural proof.** Against the pre-review guard (`0b941b27d`, `startswith("from ")`), a
programmatically-set value `settings["database.directory"] = <existing dir>` **raised
`SettingsError`** — the false positive. Against the reviewed guard (`== DEFAULT_SOURCE`) the same call
**returns the directory**. Confirmed both directions directly.

**Trap for the next person (no code change, report note).** `test/rmgpy/rmg/inputTest.py` run in a
fresh worktree with **no `rmgrc` present** yields **22 failures** — all from the PRE-EXISTING
no-`rmgrc` refusal (door 1, shipped in i149's `fccda5fc0`), triggered because every `database()` call
resolves `require_database_directory()` and finds no file. They are **not** from this change. Place a
valid `rmgrc` in cwd first; then baseline (`b12ef4c43`) and branch are **both 101 passed / 1 skipped**
(measured both). This trap cost a review diagnosis; it will cost the next person one.

Commands and raw output are recorded in the contract note `docs/contracts/i198-rmgrc-guard.md` and
this ticket's conversation. The new test file is `test/rmgpy/i198RmgrcGuardTest.py`.

## What was deliberately not done

- The loose substring match (Q3) — separate weakness, out of scope, argued above.
- `test_data.directory` (Q2) — same hole, no guard, off the run path, out of scope.
- Widening the guard to Arkane / direct readers (Q4) — the guard stays on the one path that consults it.
- No fallback, no default-guessing, no format/`configparser` rewrite (forbidden by the ticket).
