"""Output-only bounded Arrhenius fits for tabulated kinetics."""

import pytest

from rmgpy.kinetics import KineticsData
from rmgpy.kinetics.export import fit_kinetics_data_for_export


def test_fit_uses_only_source_grid_sets_bounds_and_reports_error_without_mutating_table():
    table = KineticsData(
        Tdata=([300, 400, 500, 700, 1000], "K"),
        kdata=([1e-18, 1.5e-18, 2.0e-18, 3.2e-18, 5e-18], "m^3/(mol*s)"),
        Tmin=(250, "K"), Tmax=(1500, "K"),
        comment="Source: fixture table; SHA-256 abc",
    )

    fit = fit_kinetics_data_for_export(table, n_gas=2)

    assert (fit.Tmin.value_si, fit.Tmax.value_si) == (300, 1000)
    assert "Source: fixture table" in fit.comment
    assert "maximum relative fit error" in fit.comment
    assert "The source table is the record" in fit.comment
    assert isinstance(table, KineticsData)
    assert (table.Tmin.value_si, table.Tmax.value_si) == (250, 1500)


def test_non_arrhenius_table_is_refused_when_fit_error_exceeds_25_percent():
    table = KineticsData(
        Tdata=([300, 400, 500, 600, 700, 800], "K"),
        kdata=([1e-18, 1e-14, 1e-18, 1e-14, 1e-18, 1e-14], "m^3/(mol*s)"),
        Tmin=(300, "K"), Tmax=(800, "K"),
    )

    with pytest.raises(ValueError, match="exceeds 25%"):
        fit_kinetics_data_for_export(table, n_gas=2)
