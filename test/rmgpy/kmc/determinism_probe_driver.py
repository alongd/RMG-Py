"""Cold-compile a reduced fixture under one hidden-state perturbation (I-122).

Run in a fresh interpreter per perturbation; the artifact lands in OUTPUT.

    python determinism_probe_driver.py DATABASE OUTPUT PERTURBATION \
        --families F1,F2 --proxies P1,P2      # or --proxies all (full PS set)

Perturbations
-------------
baseline                 nothing
counter_start:N          RMG's atom-id counter at N before the compile
counter_at_proxies:N     the counter at N right before the proxies are built,
                         so N near 32767 makes the wrap land inside them
counter_before_each:N    the counter advanced by N before every generation
permute_atoms            every participant's atom list shuffled per call
reverse_atoms            every participant's atom list reversed per call
reverse_reactants        the reactant list reversed per call
shuffle_output           the generator's returned reaction list shuffled
reverse_output           the generator's returned reaction list reversed
shuffle_degenerate       the list entering RMG's degeneracy collapse shuffled
heap_noise               objects allocated and freed around every call
proxy_flag_on            every participant flagged as a polymer proxy
species_proxy_flag_only  the Species flag set without its molecules'
clear_labels             every label cleared on every returned reaction's participants
scatter_labels           each returned participant form keeps or loses its reaction
                         labels by a coin flip (RMG leaves a run-dependent subset)
preset_sorting / preset_connectivity / preset_ignore_h
                         scratch vertex fields pre-set on every participant
tolerance:N[:M]          missing-reverse tolerance N (M dropped reactions
                         already counted)
"""

import argparse
import hashlib
import random
import sys
import time
from pathlib import Path

from rmgpy.kmc import compile as compile_module
from rmgpy.kmc import compiler as compiler_module
from rmgpy.kmc.compile import CompileRequest, compile_artifact
from rmgpy.molecule.molecule import Atom, Molecule


def counter_now() -> int:
    """Read RMG's atom-id counter (consumes one id)."""
    molecule = Molecule()
    molecule.add_atom(Atom(element="C"))
    molecule.assign_atom_ids()
    return molecule.atoms[0].id


def advance_counter(count: int) -> None:
    count %= 65536
    block = Molecule()
    for _ in range(1000):
        block.add_atom(Atom(element="C"))
    while count >= 1000:
        block.assign_atom_ids()
        count -= 1000
    if count:
        rest = Molecule()
        for _ in range(count):
            rest.add_atom(Atom(element="C"))
        rest.assign_atom_ids()


def set_counter(target: int) -> None:
    advance_counter(target - (counter_now() + 1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("database")
    parser.add_argument("output")
    parser.add_argument("perturbation")
    parser.add_argument("--families", required=True)
    parser.add_argument("--proxies", required=True)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    perturbation = args.perturbation
    families = tuple(args.families.split(","))
    wanted = None if args.proxies == "all" else set(args.proxies.split(","))
    rng = random.Random(args.seed)
    heap = []

    original_proxy_set = compile_module.ps_proxy_set

    def proxy_set(*options, **keywords):
        if perturbation.startswith("counter_at_proxies:"):
            set_counter(int(perturbation.split(":")[1]))
        return tuple(
            proxy for proxy in original_proxy_set(*options, **keywords)
            if wanted is None or proxy.site_type in wanted
        )

    compile_module.ps_proxy_set = proxy_set

    if perturbation.startswith("tolerance:"):
        from rmgpy.data.kinetics import family as family_module

        parts = perturbation.split(":")
        family_module.TOLERATE_MISSING_REVERSE = int(parts[1])
        if len(parts) > 2:
            family_module._dropped_reverse_counts["preseeded"] = int(parts[2])

    if perturbation == "shuffle_degenerate":
        import rmgpy.data.kinetics.database as database_module

        collapse = database_module.find_degenerate_reactions

        def shuffled_collapse(reactions, *options, **keywords):
            reactions = list(reactions)
            rng.shuffle(reactions)
            return collapse(reactions, *options, **keywords)

        database_module.find_degenerate_reactions = shuffled_collapse

    def factory(default):
        def provider(reactants, products=None, only_families=None, resonance=True):
            if perturbation.startswith("counter_before_each:"):
                advance_counter(int(perturbation.split(":")[1]))
            if perturbation == "heap_noise":
                for _ in range(rng.randrange(1, 4000)):
                    heap.append(rng.choice([Atom(element="C"), {}, [], Molecule(), object()]))
                if rng.random() < 0.3:
                    del heap[: rng.randrange(0, len(heap) + 1)]
            if perturbation == "reverse_reactants":
                reactants = list(reversed(reactants))
            for species in reactants:
                if perturbation == "proxy_flag_on":
                    species.is_polymer_proxy = True
                if perturbation == "species_proxy_flag_only":
                    species._is_polymer_proxy = True
                for molecule in species.molecule:
                    if perturbation == "permute_atoms":
                        atoms = list(molecule.atoms)
                        rng.shuffle(atoms)
                        molecule.atoms = atoms
                    if perturbation == "reverse_atoms":
                        molecule.atoms = list(reversed(molecule.atoms))
                    for index, atom in enumerate(molecule.atoms):
                        if perturbation == "preset_sorting":
                            atom.sorting_label = rng.randrange(-5, 500)
                        elif perturbation == "preset_connectivity":
                            atom.connectivity1 = rng.randrange(0, 50)
                            atom.connectivity2 = rng.randrange(0, 50)
                            atom.connectivity3 = rng.randrange(0, 50)
                        elif perturbation == "preset_ignore_h":
                            atom.ignore = atom.is_hydrogen() and index % 3 == 0
            reactions = default(reactants, products, only_families, resonance)
            if perturbation in ("clear_labels", "scatter_labels"):
                for reaction in reactions:
                    for participant in list(reaction.reactants) + list(reaction.products):
                        for molecule in getattr(participant, "molecule", [participant]):
                            if perturbation == "clear_labels" or rng.random() < 0.5:
                                molecule.clear_labeled_atoms()
            if perturbation == "shuffle_output":
                reactions = list(reactions)
                rng.shuffle(reactions)
            if perturbation == "reverse_output":
                reactions = list(reversed(reactions))
            return reactions

        return provider

    if perturbation.startswith("counter_start:"):
        set_counter(int(perturbation.split(":")[1]))
    before = counter_now()
    start = time.time()
    path, _ = compile_artifact(
        CompileRequest(
            database_path=args.database,
            output_dir=args.output,
            families=families,
            family_universe=families,
            database_sha="determinism-probe",
            legacy_provider="polystyrene",
        ),
        generation_provider_factory=factory,
    )
    digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    print(
        f"RESULT {perturbation} {digest} {time.time() - start:.1f}s "
        f"counter {before} -> {counter_now()}",
        flush=True,
    )


if __name__ == "__main__":
    sys.exit(main())
