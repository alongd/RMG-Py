"""Bounded, explicitly lossy kinetics conversions used only by mechanism writers."""

import numpy as np

from rmgpy.kinetics.arrhenius import Arrhenius
from rmgpy.kinetics.kineticsdata import KineticsData
from rmgpy.kinetics.model import get_rate_coefficient_units_from_reaction_order


MAX_KINETICS_DATA_EXPORT_RELATIVE_ERROR = 0.25


def fit_kinetics_data_for_export(kinetics, n_gas, n_surf=0):
    """Fit table data over its own grid for output; never mutate solver kinetics."""
    if not isinstance(kinetics, KineticsData):
        raise TypeError("Expected KineticsData for output-only Arrhenius fitting")

    temperatures = kinetics.Tdata.value_si
    rates = kinetics.kdata.value_si
    if len(temperatures) < 4:
        raise ValueError("At least four KineticsData points are required for a three-parameter Arrhenius fit")

    fit = Arrhenius()
    units = get_rate_coefficient_units_from_reaction_order(n_gas, n_surf)
    fit.fit_to_data(temperatures, rates, units)
    fit.Tmin = (float(temperatures[0]), "K")
    fit.Tmax = (float(temperatures[-1]), "K")
    relative_errors = np.abs(
        np.array([fit.get_rate_coefficient(float(t)) for t in temperatures]) - rates
    ) / np.abs(rates)
    max_index = int(np.argmax(relative_errors))
    max_error = float(relative_errors[max_index])
    if max_error > MAX_KINETICS_DATA_EXPORT_RELATIVE_ERROR:
        raise ValueError(
            "KineticsData output fit maximum relative error "
            f"{100 * max_error:.2f}% exceeds 25% at {temperatures[max_index]:g} K"
        )

    record = kinetics.comment.strip()
    if record:
        record += "\n"
    fit.comment = (
        f"{record}Output-only modified Arrhenius fit over the source table range "
        f"{temperatures[0]:g}-{temperatures[-1]:g} K; maximum relative fit error "
        f"{100 * max_error:.2f}% at {temperatures[max_index]:g} K. The source table is the record."
    )
    return fit
