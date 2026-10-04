"""Full-molecule C3 oracle, independent of compiler structures and scheduling."""

import hashlib
import json
import os
import subprocess
from collections import defaultdict
from pathlib import Path

from rmgpy.molecule.molecule import Molecule
from rmgpy.molecule.resonance import generate_aromatic_resonance_structure
from rmgpy.species import Species


CANDIDATE_FAMILIES = (
    "Disproportionation",
    "H_Abstraction",
    "R_Addition_MultipleBond",
    "R_Recombination",
    "intra_H_migration",
)


def full_molecule_inputs(radius):
    units = 2 * radius + 3
    structures = {
        "pristine": "CC(c1ccccc1)" * units,
        "interior_radical": "".join(
            "C[C](c1ccccc1)" if index == units // 2 else "CC(c1ccccc1)"
            for index in range(units)
        ),
        "doubly_featured": "C[C](c1ccccc1)" * 2 + "CC(c1ccccc1)" * (units - 2),
        "end_radical": "[CH2]C(c1ccccc1)" + "CC(c1ccccc1)" * (units - 1),
        "benzylic_end_radical": "CC(c1ccccc1)" * (units - 1) + "C[CH](c1ccccc1)",
        "benzylic_end_radical_short": "CC(c1ccccc1)" * (units - 2) + "C[CH](c1ccccc1)",
        "junction_radical": "CC(c1cc[c]cc1)" + "CC(c1ccccc1)" * (units - 1),
        "styrene": "C=Cc1ccccc1",
    }
    species = {
        name: Species(molecule=[Molecule(smiles=smiles)])
        for name, smiles in structures.items()
    }
    declarations = [(name, (name,)) for name in structures
                    if name not in ("styrene", "benzylic_end_radical_short")]
    declarations.extend(
        [
            ("end_radical+end_radical", ("end_radical", "end_radical")),
            ("junction_radical+end_radical", ("junction_radical", "end_radical")),
            ("end_radical+styrene", ("end_radical", "styrene")),
            ("benzylic_end_radical+styrene", ("benzylic_end_radical_short", "styrene")),
            ("benzylic_end_radical+end_radical", ("benzylic_end_radical", "end_radical")),
            ("benzylic_end_radical+benzylic_end_radical", ("benzylic_end_radical", "benzylic_end_radical")),
            ("junction_radical+benzylic_end_radical", ("junction_radical", "benzylic_end_radical")),
        ]
    )
    declarations.extend(
        (f"{radical}+pristine", (radical, "pristine"))
        for radical in (
            "interior_radical",
            "end_radical",
            "benzylic_end_radical",
            "junction_radical",
        )
    )
    return [
        (site, tuple(species[name] for name in participants))
        for site, participants in declarations
    ]


def oracle_cache_key(repository, database, radius):
    commits = [subprocess.check_output(
        ["git", "-C", str(repository), "rev-parse", "HEAD"], text=True
    ).strip()]
    commits.append(os.environ.get("RMG_DATABASE_SHA") or subprocess.check_output(
        ["git", "-C", str(database), "rev-parse", "HEAD"], text=True
    ).strip())
    files = [
        Path(__file__),
        Path(repository) / "rmgpy/kmc/compiler.py",
        Path(repository) / "rmgpy/kmc/reference_thermo.py",
    ]
    digest = hashlib.sha256(b"".join(path.read_bytes() for path in files)).hexdigest()
    return "-".join([*commits, digest, str(radius)])


def species_graph_key(species):
    molecules = getattr(species, "molecule", [species])
    forms = []
    for molecule in molecules:
        aromatic = generate_aromatic_resonance_structure(
            molecule, copy=True, save_order=True
        )
        forms.append((aromatic[0] if aromatic else molecule).to_smiles())
    return min(forms)


def reaction_graph_key(family, reactants, products):
    sides = sorted(
        (
            tuple(sorted(species_graph_key(species) for species in reactants)),
            tuple(sorted(species_graph_key(species) for species in products)),
        )
    )
    return family, tuple(sides)


def independent_c3_oracle(kinetics_database, repository, database, radius, cache_root):
    path = Path(cache_root) / (oracle_cache_key(repository, database, radius) + ".json")
    if path.is_file():
        return json.loads(path.read_text())
    degeneracies = defaultdict(float)
    reactions = set()
    counts = {}
    for site, participants in full_molecule_inputs(radius):
        generated = []
        for family in CANDIDATE_FAMILIES:
            generated.extend(
                kinetics_database.generate_reactions_from_families(
                    [participant.copy(deep=True) for participant in participants],
                    only_families=[family],
                    resonance=True,
                )
            )
        counts[site] = len(generated)
        for reaction in generated:
            template = ";".join(
                getattr(item, "label", str(item)) for item in (reaction.template or [])
            )
            degeneracies[(site, reaction.family, template)] += float(
                reaction.degeneracy
            )
            graph_key = reaction_graph_key(
                reaction.family, reaction.reactants, reaction.products
            )
            reactions.add((site, reaction.family, template, graph_key[1]))
    oracle = {
        "units": 2 * radius + 3,
        "counts": counts,
        "keys": [list(key) for key in sorted(degeneracies)],
        "reactions": sorted(reactions),
        "degeneracies": [
            {"key": list(key), "value": value}
            for key, value in sorted(degeneracies.items())
        ],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(oracle, sort_keys=True))
    return oracle
