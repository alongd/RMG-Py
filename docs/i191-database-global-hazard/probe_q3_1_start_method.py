#!/usr/bin/env python3
"""
Q3-1: what multiprocessing start method is actually in force for a plain
(non-pytest) interpreter under this environment, in this worktree.

This is an in-process measurement (no Pool is created; get_start_method()
only reports what would be used).
"""
import multiprocessing

import rmgpy

print('database.directory =', rmgpy.settings['database.directory'])
print('multiprocessing.get_start_method() =', multiprocessing.get_start_method())
print('multiprocessing.get_start_method(allow_none=True) =',
      multiprocessing.get_start_method(allow_none=True))
print('platform default via get_all_start_methods() =', multiprocessing.get_all_start_methods())
