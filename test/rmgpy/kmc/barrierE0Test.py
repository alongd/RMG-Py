"""Contract tests for the compiler-local zero-K barrier provider."""

import math

import pytest

from rmgpy.kmc.barrier_e0 import FixedBBarrierE0Provider
from rmgpy.kmc.compiler import EventSetCompiler
from rmgpy.kinetics.arrhenius import ArrheniusEP
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import ThermoData


def _thermo(*, E0=None):
    return ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
        Cpdata=([70, 95, 115, 130, 150, 160, 175], "J/(mol*K)"),
        H298=(50, "kJ/mol"),
        S298=(250, "J/(mol*K)"),
        Cp0=(30, "J/(mol*K)"),
        CpInf=(200, "J/(mol*K)"),
        E0=E0,
    )


def _reaction(thermo):
    reactant = Species(label="reactant")
    reactant.thermo = thermo
    product = Species(label="product")
    product.thermo = _thermo()
    return Reaction(reactants=[reactant], products=[product])


@pytest.mark.parametrize("value", [0.0, -1.0, math.inf, math.nan, True])
def test_fixed_b_provider_rejects_invalid_b(value):
    with pytest.raises((TypeError, ValueError), match="finite positive number"):
        FixedBBarrierE0Provider(value)


def test_provider_fits_missing_e0_with_the_declared_fixed_b():
    shared = _thermo()
    reaction = _reaction(shared)
    expected = shared.to_wilhoit(B=900.0).E0.value_si
    optimized = shared.to_wilhoit().E0.value_si
    assert expected != pytest.approx(optimized, abs=1.0)

    prepared, assignments = FixedBBarrierE0Provider(900.0).prepare_reaction(
        reaction
    )

    assert prepared.reactants[0].thermo.E0.value_si == pytest.approx(expected)
    assert assignments[0] == {
        "role": "reactant",
        "index": 0,
        "origin": "provider",
        "E0_J_per_mol": pytest.approx(expected),
        "temperature_grid_K": [300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0],
        "weights": [1.0] * 7,
    }


def test_provider_never_mutates_shared_thermo():
    shared = _thermo()
    reaction = _reaction(shared)

    prepared, _ = FixedBBarrierE0Provider(900.0).prepare_reaction(reaction)

    assert reaction.reactants[0].thermo is shared
    assert reaction.reactants[0].thermo.E0 is None
    assert reaction.products[0].thermo.E0 is None
    assert prepared.reactants[0].thermo is not shared
    assert prepared.reactants[0].thermo.E0 is not None


def test_provider_never_overwrites_supplied_e0():
    supplied = _thermo(E0=(12.5, "kJ/mol"))
    reaction = _reaction(supplied)

    prepared, assignments = FixedBBarrierE0Provider(900.0).prepare_reaction(
        reaction
    )

    assert prepared.reactants[0].thermo.E0.value_si == pytest.approx(12500.0)
    assert assignments[0]["origin"] == "supplied"
    assert assignments[0]["temperature_grid_K"] is None
    assert assignments[0]["weights"] is None


def test_provider_on_rejects_unfittable_thermo_without_fallback():
    malformed = _thermo()
    malformed.Cp0 = None

    with pytest.raises(ValueError, match="fixed-B E0 fit failed"):
        FixedBBarrierE0Provider(900.0).prepare_reaction(_reaction(malformed))


def test_compiler_uses_provider_only_for_the_barrier_copy():
    reactant_thermo = _thermo()
    product_thermo = ThermoData(
        Tdata=([300, 400, 500, 600, 800, 1000, 1500], "K"),
        Cpdata=([95, 122, 143, 158, 177, 188, 202], "J/(mol*K)"),
        H298=(65, "kJ/mol"),
        S298=(270, "J/(mol*K)"),
        Cp0=(30, "J/(mol*K)"),
        CpInf=(230, "J/(mol*K)"),
    )
    reaction = _reaction(reactant_thermo)
    reaction.products[0].thermo = product_thermo
    reaction.kinetics = ArrheniusEP(
        A=(1.0e6, "s^-1"), n=0.0, alpha=0.0, E0=(0.0, "kJ/mol")
    )

    class Assignment:
        def assign_reaction(self, assigned):
            assert assigned is reaction
            return [{"role": "reactant"}, {"role": "product"}]

    compiler = EventSetCompiler.__new__(EventSetCompiler)
    compiler.temperature_grid = (300.0, 600.0)
    compiler.thermo_assignment = Assignment()
    compiler.thermo_database = object()
    compiler.barrier_e0_provider = FixedBBarrierE0Provider(900.0)

    table, source = compiler._rate_table(reaction)

    assert table is not None
    conversion = source["kinetics_conversion"]
    assert conversion["barrier_e0_provider"] == compiler.barrier_e0_provider.provenance
    assert [item["origin"] for item in conversion["barrier_e0_assignments"]] == [
        "provider",
        "provider",
    ]
    assert reaction.reactants[0].thermo is reactant_thermo
    assert reaction.products[0].thermo is product_thermo
    assert reactant_thermo.E0 is None
    assert product_thermo.E0 is None
