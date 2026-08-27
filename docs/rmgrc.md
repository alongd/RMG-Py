# The `rmgrc` run-configuration file (plasma lineage)

RMG resolves `database.directory` from an `rmgrc` file. On this lineage that file is **untracked
on purpose** and must stay that way. This note is the thing to read before you wonder where your
`rmgrc` went or why a run stopped.

## Why it is untracked

A tracked `rmgrc` pins `database.directory` at a path that is correct for exactly one checkout.
Every branch cut from a branch that tracks it inherits that pointer verbatim — at a database the
new branch does not own. That produced four consecutive integration passes whose green test
suites measured the wrong (or a since-moved) database, and it took a third pass to notice. Each
pass "fixed" its own pointer and committed the fix, which is precisely what made the next branch
inherit a wrong one. See the repointing commits: `fa3b58f05` (origin pin) then `ae52de38a`,
`9874fc706`, `34de41c76`, `cd6cde425`.

The file is now git-ignored (`.gitignore`: `/rmgrc`) and a copy-me template ships beside it
(`rmgrc.template`).

## Fresh-worktree setup (the whole thing)

```bash
cp rmgrc.template rmgrc          # from the repo root of your worktree
# edit ./rmgrc so database.directory points at YOUR RMG-database checkout
# (plasma worktrees sit beside RMG-database-plasma, so the template's default is usually right)
```

If you skip this, the run stops in its first second naming the file it wanted — by design. There
is no fallback to a default database, because a silent default is exactly what caused the
meaningless greens above.

## What fails, and how, when it is wrong

Enforced in `Settings.require_database_directory()` (`rmgpy/__init__.py`), called from the
`database()` directive while the input file is read (`rmgpy/rmg/input.py`):

- **No `rmgrc` found anywhere** (`./rmgrc`, `~/.rmg/rmgrc`, `rmgpy/rmgrc`) → `SettingsError`,
  non-zero exit, naming the file and telling you to copy the template.
- **`database.directory` names a non-existent directory** → `SettingsError`, non-zero exit,
  naming the resolved path and the file it came from. No fall-back, no sibling search, no guess.

A correctly-configured run logs the resolved path (`Using RMG database at: <path> (<source>)`)
and otherwise behaves exactly as before.

Note: because the guard fires inside `database()`, unit tests that call that directive now require
a configured `rmgrc` too (the worktree's `./rmgrc` satisfies this). On a checkout with no `rmgrc`
at all, those tests fail loudly rather than silently reading a default database — the same policy,
applied to the suite.

## Merge and checkout consequences

1. **Merging a branch that still carries pin commits into a branch where the file is untracked
   and ignored.** Git sees the file *modified* on the incoming side and *deleted* on ours — a
   modify/delete conflict, not a silent overwrite. Resolve it by keeping the file **removed**
   (`git rm rmgrc` during the merge); do not accept the incoming tracked copy. `.gitignore` does
   not auto-resolve this: ignore rules never override an explicit tree entry coming in through a
   merge.

2. **Checking out one of the ~31 existing branches that still track `rmgrc`, after this lands.**
   You get that branch's **stale tracked copy** on disk — git restores tracked files from the
   checked-out commit regardless of `.gitignore`, because ignore rules only affect *untracked*
   files. So an old branch still hands you its old pin. Untracking on this branch does **not**
   retroactively clean the others; a branch stops carrying the file only once it merges or rebases
   past this change and drops the file itself. Until then, re-point or delete the restored `rmgrc`
   by hand after checking such a branch out.

3. **What a worker must do on a fresh worktree.** The two lines under *Fresh-worktree setup*
   above: copy the template to `./rmgrc`, set `database.directory` to your own checkout.
