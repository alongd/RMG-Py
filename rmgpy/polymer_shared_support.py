"""Fail-closed support boundary for shared-spec method-of-moments pools.

Only operations with an exact end-inventory update and a production residual
mass-closure test belong in the allow-list below.  Adding a kernel or route is
therefore a two-part change: implement its E law, then name its proof test here.
Everything absent from the allow-list is refused at configuration time.
"""

from typing import Iterable, Mapping, Optional, Sequence


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
    matches = any(
        candidate.is_isomorphic(reference)
        and int(getattr(candidate, "multiplicity", 1))
        == int(getattr(reference, "multiplicity", 1))
        for candidate in actual
        for reference in expected
    )
    if not matches:
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
        if core_species is not None:
            for dp in sorted(required_dp):
                _verify_explicit_graph(pool, dp, core_species)
