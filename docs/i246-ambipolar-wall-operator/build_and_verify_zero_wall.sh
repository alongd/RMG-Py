#!/usr/bin/env bash
# Build a second extension module from the PRE-CHANGE plasma.pyx and compare it,
# bit for bit, against this branch's build on a wall-less reactor.
#
# The temporary source, its generated C and its .so are removed on exit, so the
# worktree is left exactly as it was found.
set -euo pipefail

BASE_SHA=311818121
ROOT=/home/alon/Code/RMG-Py-i246-ambipolar-wall-operator
HERE="$ROOT/docs/i246-ambipolar-wall-operator"
TMPNAME=plasma_base_i246

export PATH=/home/alon/anaconda3/envs/rmg_env/bin:$PATH
cd "$ROOT"

cleanup() {
    rm -f "rmgpy/solver/${TMPNAME}.pyx" "rmgpy/solver/${TMPNAME}.c" \
          "rmgpy/solver/${TMPNAME}".*.so "$HERE/_build_base_ext.py"
    rm -rf "build/temp.linux-x86_64-cpython-39/rmgpy/solver/${TMPNAME}.o"
}
trap cleanup EXIT

git show "${BASE_SHA}:rmgpy/solver/plasma.pyx" > "rmgpy/solver/${TMPNAME}.pyx"
echo "extracted pre-change source: $(wc -l < rmgpy/solver/${TMPNAME}.pyx) lines from ${BASE_SHA}"

cat > "$HERE/_build_base_ext.py" <<PYEOF
import numpy
from setuptools import setup, Extension
from Cython.Build import cythonize

setup(
    ext_modules=cythonize(
        [Extension("rmgpy.solver.${TMPNAME}",
                   ["rmgpy/solver/${TMPNAME}.pyx"],
                   include_dirs=[".", numpy.get_include()])],
        compiler_directives={"language_level": 3},
    ),
    script_args=["build_ext", "--inplace"],
)
PYEOF

python "$HERE/_build_base_ext.py"
echo "built rmgpy/solver/${TMPNAME}"

cd "$HERE"
PYTHONPATH="$ROOT:." python verify_zero_wall_bitwise.py
