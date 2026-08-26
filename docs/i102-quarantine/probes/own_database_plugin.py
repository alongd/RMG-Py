"""
A pytest plugin that points ``database.directory`` at this ticket's data worktree.

The committed ``rmgrc`` at the code repo root pins the database relatively, to the
shared plasma checkout, and that file is not the place to record a per-ticket
override. Running the suite under this plugin exercises the same tests against the
database that actually carries the quarantine manifest::

    pytest test/rmgpy/ -p own_database_plugin

Without it the manifest-dependent tests skip, loudly, naming the database they looked
in -- which is the correct behaviour on any checkout that predates the marker.
"""

DATABASE = "/home/alon/Code/RMG-database-i102-quarantine/input"


def pytest_configure(config):
    from rmgpy import settings

    settings["database.directory"] = DATABASE
