from contextlib import contextmanager
from types import SimpleNamespace

import numpy as np

from rmgpy.kinetics import Arrhenius
from rmgpy.data import rmg as rmg_data
from rmgpy.molecule import Molecule
from rmgpy.species import Species

from examples.rmg.polymer_rate_estimation import (
    ProductionRateEstimator,
    ProductionReaction,
    TEMPERATURES,
    _requested_orientation,
    _requested_view,
    fit_channel,
    write_artifacts,
)


def make_species(smiles):
    return Species(molecule=[Molecule().from_smiles(smiles)])


@contextmanager
def restored_rmg_database():
    previous_database = rmg_data.database
    try:
        yield
    finally:
        rmg_data.database = previous_database


def test_requested_orientation_matches_unordered_sides():
    ethane = make_species("CC")
    methyl = make_species("[CH3]")
    reaction = SimpleNamespace(
        reactants=[ethane],
        products=[methyl, methyl],
    )

    forward, products = _requested_orientation(reaction, [Molecule().from_smiles("CC")])
    assert forward
    assert products == [methyl, methyl]

    forward, products = _requested_orientation(
        reaction,
        [Molecule().from_smiles("[CH3]"), Molecule().from_smiles("[CH3]")],
    )
    assert not forward
    assert products == [ethane]


def test_requested_view_uses_reverse_template_and_degeneracy():
    ethane = make_species("CC")
    methyl = make_species("[CH3]")
    reaction = SimpleNamespace(
        reactants=[ethane],
        products=[methyl, methyl],
        degeneracy=1.0,
        template=[SimpleNamespace(label="forward")],
        reverse=SimpleNamespace(
            degeneracy=2.0,
            template=[SimpleNamespace(label="reverse")],
        ),
    )

    view = _requested_view(
        reaction,
        [Molecule().from_smiles("[CH3]"), Molecule().from_smiles("[CH3]")],
    )

    assert view.products == [ethane]
    assert view.degeneracy == 2.0
    assert view.template_labels == ["reverse"]


class ConstantKinetics:
    def get_rate_coefficient(self, temperature):
        return 8.0


class ConstantEquilibriumReaction:
    kinetics = ConstantKinetics()

    def get_equilibrium_constant(self, temperature, type):
        assert type == "Kc"
        return 4.0


def test_production_reaction_preserves_requested_physical_direction():
    reaction = ConstantEquilibriumReaction()
    forward = ProductionReaction(
        reaction=reaction,
        requested_forward_in_final=True,
        selected_kinetics={"requested_degeneracy": 3.0},
    )
    reverse = ProductionReaction(
        reaction=reaction,
        requested_forward_in_final=False,
        selected_kinetics={"requested_degeneracy": 2.0},
    )

    assert np.array_equal(forward.requested_forward_rates([500.0]), [8.0])
    assert np.array_equal(forward.requested_reverse_rates([500.0]), [2.0])
    assert np.array_equal(reverse.requested_forward_rates([500.0]), [2.0])
    assert np.array_equal(reverse.requested_reverse_rates([500.0]), [8.0])
    assert reverse.requested_degeneracy == 2.0


def test_real_pe_own_reverse_intra_h_migration_fixture():
    """Pin the production-path own-reverse PE transfer row.

    Expected values were derived from the rmg_env production run against
    /home/alon/runs/phase2b/database/input on 2026-10-07. The tight tolerance
    catches changes in production selection or barrier processing while
    allowing floating-point noise.
    """
    with restored_rmg_database():
        estimator = ProductionRateEstimator(
            "/home/alon/runs/phase2b/database/input", ["intra_H_migration"]
        )
        reactions = estimator.generate(["[CH2]CCCCCCC"], "intra_H_migration")
    own_reverse = [
        reaction for reaction in reactions if not reaction.requested_forward_in_final
    ]

    assert len(reactions) == 6
    assert len(own_reverse) == 1
    reaction = own_reverse[0]
    assert reaction.requested_degeneracy == 2.0
    assert reaction.reaction.degeneracy == 3.0
    assert reaction.selected_kinetics["requested_template_labels"] == [
        "R6H_SSSSS",
        "C_rad_out_2H",
        "Cs_H_out_H/(NonDeC/Cs)",
    ]
    assert reaction.selected_kinetics["final_template_labels"] == [
        "R3H_SS_Cs",
        "C_rad_out_H/NonDeC",
        "Cs_H_out_2H",
    ]
    np.testing.assert_allclose(
        reaction.requested_forward_rates([300.0, 600.0, 1000.0]),
        [9.076264062166428e-15, 0.24188767379152315, 62944.04667768389],
        rtol=1e-8,
        atol=0.0,
    )


def test_fit_and_artifacts_report_rmg_and_solver_laws(tmp_path):
    source = Arrhenius(
        A=(2.5e7, "s^-1"),
        n=0.75,
        Ea=(45000.0, "J/mol"),
        T0=(1.0, "K"),
    )
    rates = np.array(
        [source.get_rate_coefficient(temperature) for temperature in TEMPERATURES]
    )
    channel = fit_channel(
        "example",
        rates,
        "s^-1",
        {"family": "test", "surrogate": "A -> B"},
    )

    assert np.isclose(channel["A"], 2.5e7)
    assert np.isclose(channel["n"], 0.75)
    assert np.isclose(channel["Ea_J_mol"], 45000.0)
    assert channel["max_relative_fit_error"] < 1e-12
    assert channel["max_relative_qssa_law_error"] > 0.0

    write_artifacts(tmp_path, tmp_path / "database", [channel])
    header = (tmp_path / "rate_points.csv").read_text().splitlines()[0]
    assert header == "channel,T_K,k_RMG,k_fit,relative_error"
    table = (tmp_path / "rate_table.md").read_text()
    assert "max RMG fit error" in table
    assert "max solver-law error" in table
