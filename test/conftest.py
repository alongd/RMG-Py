import multiprocessing

multiprocessing.set_start_method('fork')

import os

from openbabel import pybel

pybel.ob.obErrorLog.SetOutputLevel(0)


################################################################################
# Guard: the test suite must not modify the RMG-database checkout it is run
# against. Tests that need to exercise a save/load round trip must do it in a
# temporary directory; writing into settings["database.directory"] leaves the
# user's database checkout dirty and makes unrelated work look broken.
#
# Implemented as a session-level snapshot rather than an autouse per-test
# fixture for two reasons: stat-ing ~1800 files once per session is free
# whereas once per test is not, and tests that legitimately create and then
# remove files under the database directory (see TestMain in
# test/rmgpy/rmg/mainTest.py, which writes seed libraries and rmtree's them in
# teardown) are judged on the state they leave behind, which is what actually
# harms the checkout.
################################################################################

_DATABASE_SNAPSHOT = None
_DATABASE_DIRECTORY = None


def _snapshot_database_directory(directory):
    """Map every file under `directory` to a cheap (size, mtime) fingerprint.

    `__pycache__` is skipped: several database files (kinetics/families/
    recommended.py, quantum_corrections/data.py) are imported as modules, so
    CPython writes bytecode caches next to them. Those are git-ignored by
    RMG-database and regenerated on demand, so they are not pollution.
    """
    snapshot = {}
    for dirpath, dirnames, filenames in os.walk(directory):
        for ignored in ('.git', '__pycache__'):
            if ignored in dirnames:
                dirnames.remove(ignored)
        for filename in filenames:
            path = os.path.join(dirpath, filename)
            try:
                stat = os.stat(path)
            except OSError:
                continue
            snapshot[path] = (stat.st_size, stat.st_mtime_ns)
    return snapshot


def pytest_sessionstart(session):
    global _DATABASE_SNAPSHOT, _DATABASE_DIRECTORY

    from rmgpy import settings

    directory = settings["database.directory"]
    if not directory or not os.path.isdir(directory):
        return
    _DATABASE_DIRECTORY = directory
    _DATABASE_SNAPSHOT = _snapshot_database_directory(directory)


def pytest_sessionfinish(session, exitstatus):
    if _DATABASE_SNAPSHOT is None:
        return

    after = _snapshot_database_directory(_DATABASE_DIRECTORY)

    added = sorted(set(after) - set(_DATABASE_SNAPSHOT))
    removed = sorted(set(_DATABASE_SNAPSHOT) - set(after))
    modified = sorted(p for p in set(after) & set(_DATABASE_SNAPSHOT) if after[p] != _DATABASE_SNAPSHOT[p])

    if not (added or removed or modified):
        return

    lines = [
        "",
        "=" * 78,
        "DATABASE POLLUTION: the test suite modified the RMG-database checkout at",
        "    {0}".format(_DATABASE_DIRECTORY),
        "",
        "Tests must not write under settings['database.directory']. Redirect the",
        "write to a temporary directory (see test_write_entry_to_database in",
        "test/rmgpy/data/surfaceTest.py for the pattern).",
        "",
    ]
    for label, paths in (("modified", modified), ("added", added), ("removed", removed)):
        for path in paths[:20]:
            lines.append("  {0:8s} {1}".format(label, path))
        if len(paths) > 20:
            lines.append("  {0:8s} ... and {1} more".format(label, len(paths) - 20))
    lines.append("=" * 78)
    lines.append("")

    reporter = session.config.pluginmanager.get_plugin("terminalreporter")
    message = "\n".join(lines)
    if reporter is not None:
        reporter.write_line(message, red=True, bold=True)
    else:
        print(message)

    session.exitstatus = 1
