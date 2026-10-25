"""Fail-closed support boundary for shared-spec method-of-moments pools.

Only operations with an exact end-inventory update and a production residual
mass-closure test belong in the allow-list below.  Adding a kernel or route is
therefore a two-part change: implement its E law, then name its proof test here.
Everything absent from the allow-list is refused at configuration time.
"""

import math
from typing import Iterable, Mapping, Optional, Sequence


SHARED_MASS_REL_TOL = 1.0e-8
SHARED_MASS_ABS_TOL_G_MOL = 1.0e-8


SHARED_KERNEL_FIELDS = (
    "k_scission",
    "k_unzip",
    "radical_qssa_unzip",
    "k_homolysis",
    "k_depropagation",
    "side_group_homolysis",
    "tail_kinetics",
)

SHARED_KERNEL_ALLOWLIST = {
    "k_unzip": "test_shared_unzip_residual_closes_dp_above_one_and_dp1_mass",
    "side_group_homolysis": (
        "test_real_sgh_z_and_gas_accounting_matches_independent_composition"
    ),
}

SHARED_REFUSED_KERNELS = frozenset(SHARED_KERNEL_FIELDS) - frozenset(
    SHARED_KERNEL_ALLOWLIST
)

# Mirrors the solver's public FLUX_* constants.  The census test cross-checks
# this table against the compiled extension so a newly added route cannot
# silently escape validation.
SHARED_ROUTE_FLAGS = {
    0: "NONE",
    1: "SAME_POOL",
    2: "MIGRATION",
    3: "SCISSION_FRAGMENT",
    4: "UNRESOLVED",
    5: "DISCRETE_CHIP",
    6: "VOLATILE_EJECTION",
    7: "MOMENT_CREDIT_CONDUIT",
}

# Step 2 needs no generated-reaction moment route.  The exact supported
# transfer is the solver-owned unzip handshake, catalogued as a kernel above.
SHARED_ROUTE_ALLOWLIST = {}
SHARED_REFUSED_ROUTES = frozenset(SHARED_ROUTE_FLAGS)


def _kernel_is_active(pool, name: str) -> bool:
    value = getattr(pool, name, None)
    if name in ("k_scission", "k_unzip"):
        return float(value or 0.0) > 0.0
    if name == "side_group_homolysis":
        return bool(value)
    return value is not None


def _is_shared(pool) -> bool:
    return bool(
        getattr(pool, "shared_mass_accounting", False)
        or getattr(pool, "linear_homopolymer_spec", None) is not None
    )


def _expected_graphs(pool) -> Mapping[int, Sequence[object]]:
    return getattr(pool, "shared_expected_explicit_dp_molecules", None) or {}


def _same_graph(candidate, reference) -> bool:
    return (
        candidate.is_isomorphic(reference)
        and int(getattr(candidate, "multiplicity", 1))
        == int(getattr(reference, "multiplicity", 1))
    )


def _is_form_of(candidate, reference) -> bool:
    """True if ``candidate`` is ``reference`` or one of its resonance forms."""
    if _same_graph(candidate, reference):
        return True
    try:
        forms = reference.copy(deep=True).generate_resonance_structures()
    except Exception:
        return False
    return any(_same_graph(candidate, form) for form in forms)


def _graphs_match(
    actual: Sequence[object], expected: Sequence[object]
) -> bool:
    """Graph identity including isotope and electronic multiplicity.

    Every member of ``actual`` must be a resonance form of some expected
    graph: one matching member must not launder the rest of the list.
    """
    actual = tuple(actual)
    expected = tuple(expected)
    if not actual or not expected:
        return False
    return all(
        any(_is_form_of(candidate, reference) for reference in expected)
        for candidate in actual
    )


def _authoritative_spec(pool):
    spec = getattr(pool, "linear_homopolymer_spec", None)
    if spec is None:
        raise ValueError(
            f"Pool {pool.label}: shared direct config requires its "
            "authoritative LinearHomopolymerSpec for graph and inventory "
            "validation"
        )
    from rmgpy.polymer_spec import BranchedMomRefusal, is_branched, refuse_schema_v2
    if is_branched(spec):
        raise BranchedMomRefusal("polymer_shared_support._authoritative_spec")
    return refuse_schema_v2(spec, "polymer_shared_support._authoritative_spec")


def _verify_core_destination(
    pool,
    core_species: Sequence[object],
    index: int,
    expected: Sequence[object],
    description: str,
) -> None:
    if not (0 <= index < len(core_species)):
        raise ValueError(
            f"Pool {pool.label}: {description} index {index} is unavailable"
        )
    actual = tuple(getattr(core_species[index], "molecule", ()) or ())
    if not actual or not expected:
        raise ValueError(
            f"Pool {pool.label}: {description} must carry an actual "
            "molecular graph"
        )
    if not _graphs_match(actual, expected):
        raise ValueError(
            f"Pool {pool.label}: {description} is not graph-isomorphic to "
            "the authoritative declared pristine graph (including isotopes)"
        )


def _authoritative_capped_graph(spec, dp: int):
    from rmgpy.polymer import stitch_molecules_by_labeled_atoms

    repeat, ends = spec.legacy_molecules()
    molecule = ends[0].copy(deep=True)
    for _ in range(int(dp)):
        molecule = stitch_molecules_by_labeled_atoms(
            molecule, repeat.copy(deep=True)
        )
    molecule = stitch_molecules_by_labeled_atoms(
        molecule, ends[1].copy(deep=True)
    )
    if molecule is None:
        raise ValueError(
            f"Could not materialize authoritative capped DP={dp} graph"
        )
    molecule.update()
    molecule.clear_labeled_atoms()
    return molecule


def _mass_matches(supplied: float, derived: float) -> bool:
    return math.isfinite(float(supplied)) and math.isclose(
        float(supplied), float(derived),
        rel_tol=SHARED_MASS_REL_TOL,
        abs_tol=SHARED_MASS_ABS_TOL_G_MOL,
    )


def _validate_authoritative_masses_and_graphs(
    pool, core_species: Sequence[object]
) -> None:
    spec = _authoritative_spec(pool)
    repeat_graph, _ends = spec.legacy_molecules()
    repeat_mass = repeat_graph.get_molecular_weight() * 1000.0
    if not _mass_matches(pool.monomer_mw_g_mol, repeat_mass):
        raise ValueError(
            f"Pool {pool.label}: shared repeat mass "
            f"{pool.monomer_mw_g_mol!r} disagrees with authoritative graph "
            f"mass {repeat_mass:.12g} g/mol (rtol={SHARED_MASS_REL_TOL:g}, "
            f"atol={SHARED_MASS_ABS_TOL_G_MOL:g} g/mol)"
        )
    # Within tolerance the graph value, not the supplied approximation, is
    # what the live mass law must use.
    object.__setattr__(pool, "monomer_mw_g_mol", repeat_mass)
    end_mass = spec.end_composition.molar_mass_g_mol()
    if not _mass_matches(pool.end_group_mw_g_mol, end_mass):
        raise ValueError(
            f"Pool {pool.label}: shared end mass "
            f"{pool.end_group_mw_g_mol!r} disagrees with authoritative graph "
            f"mass {end_mass:.12g} g/mol (rtol={SHARED_MASS_REL_TOL:g}, "
            f"atol={SHARED_MASS_ABS_TOL_G_MOL:g} g/mol)"
        )
    object.__setattr__(pool, "end_group_mw_g_mol", end_mass)

    explicit = getattr(pool, "explicit_dp_to_species_index", None) or {}
    supplied_masses = getattr(pool, "explicit_dp_mw_g_mol", None) or {}
    supplied_expected = _expected_graphs(pool)
    for dp, index in explicit.items():
        authoritative = _authoritative_capped_graph(spec, int(dp))
        declared = tuple(supplied_expected.get(dp, ()) or ())
        if not declared or not _graphs_match(declared, (authoritative,)):
            raise ValueError(
                f"Pool {pool.label}: expected explicit DP={dp} graph "
                "disagrees with the authoritative declaration"
            )
        _verify_core_destination(
            pool, core_species, int(index), (authoritative,),
            f"shared explicit DP={dp} destination",
        )
        graph_mass = authoritative.get_molecular_weight() * 1000.0
        supplied = supplied_masses.get(dp)
        if supplied is None or not _mass_matches(supplied, graph_mass):
            raise ValueError(
                f"Pool {pool.label}: explicit DP={dp} supplied mass "
                f"{supplied!r} disagrees with actual graph mass "
                f"{graph_mass:.12g} g/mol (rtol={SHARED_MASS_REL_TOL:g}, "
                f"atol={SHARED_MASS_ABS_TOL_G_MOL:g} g/mol)"
            )
        supplied_masses[dp] = graph_mass


def _verify_explicit_graph(pool, dp: int, core_species: Sequence[object]) -> None:
    explicit = getattr(pool, "explicit_dp_to_species_index", None) or {}
    expected_by_dp = _expected_graphs(pool)
    if dp not in explicit or dp not in expected_by_dp:
        raise ValueError(
            f"Pool {pool.label}: shared k_unzip requires a declared pristine "
            f"capped DP={dp} destination graph"
        )
    index = int(explicit[dp])
    if not (0 <= index < len(core_species)):
        raise ValueError(
            f"Pool {pool.label}: capped DP={dp} destination index {index} "
            "is unavailable for shared graph validation"
        )
    actual = tuple(getattr(core_species[index], "molecule", ()) or ())
    expected = tuple(expected_by_dp[dp] or ())
    if not actual or not expected:
        raise ValueError(
            f"Pool {pool.label}: capped DP={dp} destination must carry an "
            "actual molecular graph"
        )
    if not _graphs_match(actual, expected):
        raise ValueError(
            f"Pool {pool.label}: HANDSHAKE/DP1_UNZIP destination DP={dp} "
            "is not graph-isomorphic to the declared pristine capped chain"
        )


def validate_shared_spec_allowlist(
    pool,
    *,
    context: str,
    core_species: Optional[Sequence[object]] = None,
    route_records: Iterable[Mapping[str, object]] = (),
    parent_shared_spec=None,
) -> None:
    """Validate one pool against the single fail-closed shared support table.

    ``context`` is included in errors to identify whether construction,
    system initialization, or a spawned-daughter boundary found the defect.
    ``route_records`` contains only rows that touch this pool; no Step-2
    generated-reaction route is currently admitted.
    """
    child_spec = getattr(pool, "linear_homopolymer_spec", None)
    if (child_spec is not None
            and hasattr(pool, "shared_mass_accounting")
            and not pool.shared_mass_accounting):
        raise ValueError(
            f"{context}: pool '{getattr(pool, 'label', '?')}' carries a "
            "shared spec but shared_mass_accounting is false"
        )
    if parent_shared_spec is not None:
        if child_spec is None:
            raise ValueError(
                f"{context}: shared-spec spawn would create legacy daughter "
                f"'{getattr(pool, 'label', '?')}'"
            )
        if getattr(child_spec, "normalized_json", None) \
                != getattr(parent_shared_spec, "normalized_json", None):
            raise ValueError(
                f"{context}: spawned daughter '{getattr(pool, 'label', '?')}' "
                "does not carry its parent's authoritative shared spec"
            )

    if not _is_shared(pool):
        return

    active = [name for name in SHARED_KERNEL_FIELDS
              if _kernel_is_active(pool, name)]
    refused = [name for name in active if name not in SHARED_KERNEL_ALLOWLIST]
    if refused:
        raise ValueError(
            f"Pool {pool.label}: shared-spec {context} refuses unsupported "
            f"kernel(s): {', '.join(refused)}; no exact E update has a named "
            "production residual proof"
        )
    if len(active) > 1:
        raise ValueError(
            f"Pool {pool.label}: shared-spec {context} refuses kernel "
            f"combination {', '.join(active)}; combinations require their own "
            "exact E proof"
        )

    for record in route_records:
        flag = int(record.get("flag", -1))
        name = str(record.get("name", SHARED_ROUTE_FLAGS.get(flag, flag)))
        row = record.get("row", "?")
        if flag not in SHARED_ROUTE_ALLOWLIST:
            raise ValueError(
                f"Pool {pool.label}: shared-spec {context} refuses route "
                f"{name} (row {row}); generated, unresolved, demoted, and "
                "refused routes have no proven E update"
            )

    if core_species is not None:
        _validate_authoritative_masses_and_graphs(pool, core_species)

    if active == ["k_unzip"]:
        defect = float(getattr(pool, "chain_mass_defect_g_mol", 0.0) or 0.0)
        if defect != 0.0:
            raise ValueError(
                f"Pool {pool.label}: shared k_unzip/HANDSHAKE requires zero "
                "chain_mass_defect_g_mol because removing mu0 would otherwise "
                "mint the per-chain defect mass"
            )
        if getattr(pool, "monomer_poly_index", None) is None:
            raise ValueError(
                f"Pool {pool.label}: shared k_unzip requires a real monomer "
                "credit destination"
            )
        required_dp = {1}
        xs = int(getattr(pool, "xs"))
        explicit = getattr(pool, "explicit_dp_to_species_index", None) or {}
        if xs in explicit:
            required_dp.add(xs)
        expected = _expected_graphs(pool)
        missing = sorted(dp for dp in required_dp
                         if dp not in explicit or dp not in expected)
        if missing:
            raise ValueError(
                f"Pool {pool.label}: shared k_unzip requires real declared "
                f"pristine capped destination graph(s) for DP={missing}"
            )
        spec_for_product = getattr(pool, "linear_homopolymer_spec", None)
        if (spec_for_product is not None
                and spec_for_product.monomer_product is not None
                and spec_for_product.monomer_product.composition
                != spec_for_product.repeat_unit.composition):
            raise ValueError(
                f"Pool {pool.label}: declared monomer_product composition "
                f"{spec_for_product.monomer_product.composition.as_mapping()!r}"
                " differs from the repeat unit composition "
                f"{spec_for_product.repeat_unit.composition.as_mapping()!r}; "
                "unzip would not conserve elemental mass"
            )
        if core_species is not None:
            spec = _authoritative_spec(pool)
            for dp in sorted(required_dp):
                _verify_explicit_graph(pool, dp, core_species)
            if spec.monomer_product is None:
                raise ValueError(
                    f"Pool {pool.label}: shared k_unzip requires an "
                    "authoritative declared monomer_product graph"
                )
            expected_monomer = (
                spec.monomer_product.to_legacy_molecule({}),
            )
            _verify_core_destination(
                pool,
                core_species,
                int(pool.monomer_poly_index),
                expected_monomer,
                "shared k_unzip destination for the declared monomer",
            )

    if active == ["side_group_homolysis"]:
        if core_species is not None:
            spec = _authoritative_spec(pool)
            repeat, _ends = spec.legacy_molecules()
            from rmgpy.solver.polymer import (
                side_group_site_atom_indices,
                validate_side_group_homolysis,
            )
            normalized = validate_side_group_homolysis(
                pool.label, pool.side_group_homolysis, monomer=repeat
            )
            indices = getattr(pool, "side_group_gas_indices", None) or ()
            if len(indices) != len(normalized):
                raise ValueError(
                    f"Pool {pool.label}: shared side_group_homolysis "
                    "destination indices do not cover every channel"
                )
            from rmgpy.molecule import Molecule
            for channel, index in zip(normalized, indices):
                gas_product = Molecule().from_smiles(
                    channel["gas_product"]
                )
                gas_atom = gas_product.atoms[0]
                site_indices = side_group_site_atom_indices(
                    repeat,
                    gas_atom.symbol,
                    channel["site_selector"],
                )
                source_isotopes = {
                    repeat.atoms[site_index].element.isotope
                    for site_index in site_indices
                }
                if source_isotopes != {gas_atom.element.isotope}:
                    raise ValueError(
                        f"Pool {pool.label}: shared side_group_homolysis "
                        f"channel '{channel['label']}' source isotope(s) "
                        f"{sorted(source_isotopes)} disagree with "
                        "gas_product "
                        f"isotope {gas_atom.element.isotope}"
                    )
                expected = (gas_product,)
                _verify_core_destination(
                    pool,
                    core_species,
                    int(index),
                    expected,
                    "shared side_group_homolysis destination",
                )
