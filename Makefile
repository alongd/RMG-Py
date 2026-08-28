################################################################################
#
#   Makefile for RMG Py
#
################################################################################

CC=gcc
CXX=g++

# Sentinel file written after a successful pip editable install.
# Lives in the source tree; deleted by `make clean`.
INSTALL_SENTINEL = .installed

################################################################################
#
#   Shared-environment guard
#
#   `pip install -e .` does not install into this checkout; it rewrites the
#   *environment's* editable-install record so that `import rmgpy` resolves to
#   whichever checkout ran it. Several worktrees of this repository routinely
#   share one conda environment, so an install run from one worktree silently
#   repoints every other worktree's imports at the wrong source tree.
#
#   The default goal therefore refuses instead of installing. The supported
#   route for day-to-day work is `make build`, which compiles extensions in
#   place and never touches the environment. The two targets that do mutate the
#   environment are named `unsafe-*` and require an explicit opt-in:
#
#       make unsafe-install-shared-env   CONFIRM_SHARED_ENV_MUTATION=yes
#       make unsafe-uninstall-shared-env CONFIRM_SHARED_ENV_MUTATION=yes
#
#   The guard is deliberately unconditional: it inspects neither the
#   environment's editable-install record nor any path, so it cannot be
#   defeated by that record being absent, stale, or pointed anywhere in
#   particular. `CONFIRM_SHARED_ENV_MUTATION` is intentionally left empty here.
#
################################################################################

CONFIRM_SHARED_ENV_MUTATION =

# Printed by both the default-goal guard and the opt-in check, so the refusal
# reads the same wherever it comes from. Single quotes keep the backticks
# literal rather than letting the shell run them.
define SHARED_ENV_REFUSAL
	echo 'Refusing to modify the shared RMG environment.' >&2; \
	echo 'Use `make build` for an in-place worktree build.' >&2; \
	echo 'Editable installation requires an explicit maintenance procedure:' >&2; \
	echo '    make unsafe-install-shared-env CONFIRM_SHARED_ENV_MUTATION=yes' >&2; \
	echo '' >&2; \
	echo 'Targets that do not touch the environment: build, check, clean,' >&2; \
	echo 'clean-solver, decython, documentation, q2dtor, test, test-all,' >&2; \
	echo 'test-unittests, test-functional, test-database.' >&2
endef

define REQUIRE_SHARED_ENV_CONFIRMATION
	@ if [ "$(CONFIRM_SHARED_ENV_MUTATION)" != "yes" ]; then \
		$(SHARED_ENV_REFUSAL); \
		exit 1; \
	fi
endef

.PHONY : guard all build check check-pydas clean clean-solver install \
         unsafe-install-shared-env unsafe-uninstall-shared-env \
         decython documentation test test-all test-unittests test-functional \
         test-database q2dtor scoop

# Default goal: refuse, before pip is invoked and before anything is mutated.
.DEFAULT_GOAL := guard

guard:
	@ $(SHARED_ENV_REFUSAL)
	@ exit 1

# `all` and `install` used to perform the editable install. They are kept as
# aliases of the guard so that muscle memory and old scripts fail loudly rather
# than silently repointing the environment.
all: guard

install: guard

check:
	@ python utilities.py check-dependencies
	@ python utilities.py check-pydas

# Worktree-scoped: writes only rmgpy/solver/settings.pxi in this checkout.
check-pydas:
	@ python utilities.py check-pydas

# Incremental in-place build; skips pip entirely and never touches the
# environment's site-packages or editable-install record.
build: check-pydas
	python setup.py build_ext --inplace

# Maintenance only. Rewrites the shared environment's editable-install record.
unsafe-install-shared-env:
	$(REQUIRE_SHARED_ENV_CONFIRMATION)
	@ python utilities.py check-pydas
	python -m pip install --no-build-isolation -vv -e .
	@ touch $(INSTALL_SENTINEL)

# Maintenance only. Removes the package from the shared environment, which
# breaks every other worktree relying on that editable install.
unsafe-uninstall-shared-env:
	$(REQUIRE_SHARED_ENV_CONFIRMATION)
	python -m pip uninstall --yes reactionmechanismgenerator
	@ rm -f $(INSTALL_SENTINEL)

documentation:
	$(MAKE) -C documentation html
	@ echo "Start at: documentation/build/html/index.html"

# Worktree-scoped: utilities.py clean only removes build artefacts under this
# checkout. Removing the package from the shared environment is a separate,
# opt-in target (unsafe-uninstall-shared-env).
clean:
	@ python utilities.py clean
	@ rm -f $(INSTALL_SENTINEL)

clean-solver:
	@ python utilities.py clean-solver

q2dtor:
	@ echo -e "\nInstalling Q2DTor...\n"
	@ echo -e "Q2DTor is a software for calculating the partition functions and themodynamic properties\
	of molecular systems with two or more torsional modes developed by David Ferro Costas (david.ferro@usc.es)\
	 and Antonio Fernandez Ramos (qf.ramos@usc.es) at the Universidade de Santiago de Compostela. Arkane can\
	  integrate Q2DTor to compute the quantum mechanical partition function of 2D rotors.  \n\nFor use of Q2DTor\
 and HinderedRotor2D within Arkane please cite:  \n\nD. Ferro-Costas, M. N. D. S.Cordeiro, D. G. Truhlar, A.\
		  Fernández-Ramos, Comput. Phys. Commun. 232, 190-205, 2018.\n"
	@ read -p "Press ENTER to continue" dummy
	@ mkdir -p external
	@ git clone https://github.com/cathedralpkg/Q2DTor external/Q2DTor

decython:
	# de-cythonize all but the 'minimal'. Helpful for debugging in "pure python" mode.
	find . -name *.so ! \( -name _statmech.so -o -name quantity.so -o -regex '.*rmgpy/solver/.*' \) -exec rm -f '{}' \;
	find . -name *.pyc -exec rm -f '{}' \;

test-all:
	python -m pytest

test test-unittests:
	python -m pytest -m "not functional and not database"

test-functional:
	python -m pytest -m "functional"

test-database:
	python -m pytest -m "database"

eg0: build
	mkdir -p testing/eg0
	rm -rf testing/eg0/*
	cp examples/rmg/superminimal/input.py testing/eg0/input.py
	@ echo "Running eg0: superminimal (H2 oxidation) example"
	python rmg.py testing/eg0/input.py

eg1: build
	mkdir -p testing/eg1
	rm -rf testing/eg1/*
	cp examples/rmg/minimal/input.py testing/eg1/input.py
	coverage erase
	@ echo "Running eg1: minimal (ethane pyrolysis) example with coverage tracking AND profiling"
	coverage run rmg.py -p testing/eg1/input.py
	coverage report
	coverage html

eg2: build
	mkdir -p testing/eg2
	rm -rf testing/eg2/*
	cp examples/rmg/1,3-hexadiene/input.py testing/eg2/input.py
	coverage erase
	@ echo "Running eg2: 1,3-hexadiene example with profiling"
	python rmg.py -p testing/eg2/input.py

eg3: build
	mkdir -p testing/eg3
	rm -rf testing/eg3/*
	cp examples/rmg/liquid_phase/input.py testing/eg3/input.py
	coverage erase
	@ echo "Running eg3: liquid_phase example with profiling"
	python rmg.py -p testing/eg3/input.py

eg5: build
	mkdir -p testing/eg5
	rm -rf testing/eg5/*
	cp examples/rmg/heptane-eg5/input.py testing/eg5/input.py
	@ echo "Running eg5: heptane example"
	python rmg.py testing/eg5/input.py

eg6: build
	mkdir -p testing/eg6
	rm -rf testing/eg6/*
	cp examples/rmg/ethane-oxidation/input.py testing/eg6/input.py
	@ echo "Running eg6: ethane-oxidation example"
	python rmg.py testing/eg6/input.py

eg7: build
	mkdir -p testing/eg7
	rm -rf testing/eg7/*
	cp examples/rmg/gri_mech_rxn_lib/input.py testing/eg7/input.py
	@ echo "Running eg7: gri_mech_rxn_lib example"
	python rmg.py testing/eg7/input.py

scoop: build
	mkdir -p testing/scoop
	rm -rf testing/scoop/*
	cp examples/rmg/minimal/input.py testing/scoop/input.py
	coverage erase
	@ echo "Running minimal example with SCOOP"
	python -m scoop -n 2 rmg.py -v testing/scoop/input.py

eg4: build
	mkdir -p testing/eg4
	rm -rf testing/eg4/*
	cp examples/thermoEstimator/input.py testing/eg4/input.py
	@ echo "Running thermo data estimator example. This tests QM."
	python scripts/thermoEstimator.py testing/eg4/input.py

# RMS reactor examples (require Julia)
eg8: build
	mkdir -p testing/eg8
	rm -rf testing/eg8/*
	cp examples/rmg/rms_constant_V/input.py testing/eg8/input.py
	@ echo "Running RMS constantVIdealGasReactor example (requires Julia)"
	python rmg.py testing/eg8/input.py

eg9: build
	mkdir -p testing/eg9
	rm -rf testing/eg9/*
	cp examples/rmg/nox_transitory_edge/input.py testing/eg9/input.py
	@ echo "Running RMS constantTPIdealGasReactor example (requires Julia)"
	python rmg.py testing/eg9/input.py

eg10: build
	mkdir -p testing/eg10
	rm -rf testing/eg10/*
	cp examples/rmg/liquid_cat/input.py testing/eg10/input.py
	@ echo "Running RMS liquidSurfaceReactor example (requires Julia)"
	python rmg.py testing/eg10/input.py
