"""Bounded, explicitly lossy kinetics conversions used only by mechanism writers."""

import numpy as np
from scipy.optimize import least_squares
from scipy.special import logsumexp

from rmgpy import constants
from rmgpy.kinetics.arrhenius import Arrhenius, MultiArrhenius
from rmgpy.kinetics.kineticsdata import KineticsData
from rmgpy.kinetics.model import get_rate_coefficient_units_from_reaction_order

MAX_KINETICS_DATA_EXPORT_RELATIVE_ERROR = 0.25
TARGET_KINETICS_DATA_EXPORT_RELATIVE_ERROR = 0.10


def fit_kinetics_data_for_export(kinetics, n_gas, n_surf=0):
    """Fit positive Arrhenius sums on every source node without mutating the table.

    A fixed scaled initialization makes the bounded log-rate optimization
    deterministic. The original single-term least-squares fit is preferred when
    it meets the target; additional terms are output representations only.
    """
    if not isinstance(kinetics, KineticsData):
        raise TypeError("Expected KineticsData for output-only Arrhenius fitting")
    temperatures = np.array(kinetics.Tdata.value_si, copy=True)
    rates = np.array(kinetics.kdata.value_si, copy=True)
    if (len(temperatures) < 4 or len(temperatures) != len(rates)
            or not np.all(np.isfinite(temperatures)) or not np.all(np.isfinite(rates))
            or np.any(temperatures <= 0) or np.any(rates <= 0)
            or np.any(np.diff(temperatures) <= 0)):
        raise ValueError("At least four finite positive ordered KineticsData points are required")
    units = get_rate_coefficient_units_from_reaction_order(n_gas, n_surf)
    bounds = dict(Tmin=(float(temperatures[0]), "K"), Tmax=(float(temperatures[-1]), "K"))
    single = Arrhenius().fit_to_data(temperatures, rates, units)
    single.Tmin, single.Tmax = bounds["Tmin"], bounds["Tmax"]

    def errors(fit):
        return np.abs(np.array([fit.get_rate_coefficient(float(t)) for t in temperatures]) / rates - 1)

    candidates = [(float(max(errors(single))), single)]
    if candidates[0][0] > TARGET_KINETICS_DATA_EXPORT_RELATIVE_ERROR:
        log_t = np.log(temperatures / 1000.)
        inverse_t = 1000. / temperatures
        scale = float(np.exp(np.mean(np.log(rates))))
        target = np.log(rates / scale)
        for count in (2, 3):
            # Fixed slopes/curvatures separate terms without random restarts.
            initial = np.array([[-np.log(count), 1.5, e] for e in np.linspace(0., 5., count)])
            def residual(parameters):
                a, n, e = parameters.reshape(count, 3).T
                return logsumexp(a[:, None] + n[:, None] * log_t - e[:, None] * inverse_t, axis=0) - target
            result = least_squares(residual, initial.ravel(),
                                   bounds=(np.tile([-200., -100., -100.], count),
                                           np.tile([200., 100., 100.], count)),
                                   max_nfev=4000, ftol=1e-11, xtol=1e-11, gtol=1e-11)
            # Refine the largest relative residuals rather than allowing a few
            # cold nodes to be hidden by a mean log-error objective.
            def peak_residual(parameters):
                error = np.expm1(np.clip(residual(parameters), -150., 150.))
                return error * np.abs(error)
            refined = least_squares(peak_residual, result.x,
                                    bounds=(np.tile([-200., -100., -100.], count),
                                            np.tile([200., 100., 100.], count)),
                                    max_nfev=4000, ftol=1e-11, xtol=1e-11, gtol=1e-11)
            from rmgpy.quantity import RateCoefficient
            unit_scale = RateCoefficient(1., units).value_si
            terms = [Arrhenius(A=(float(np.exp(a) * scale / unit_scale), units),
                               n=float(n), Ea=(float(e * constants.R * 1000.), "J/mol"),
                               T0=(1000., "K"), **bounds)
                     for a, n, e in refined.x.reshape(count, 3)]
            fit = MultiArrhenius(arrhenius=terms, **bounds)
            error = float(max(errors(fit)))
            candidates.append((error, fit))
            if error <= TARGET_KINETICS_DATA_EXPORT_RELATIVE_ERROR:
                break
    max_error, fit = min(candidates, key=lambda candidate: candidate[0])
    relative_errors = errors(fit)
    max_index = int(np.argmax(relative_errors))
    if not np.isfinite(max_error) or max_error > MAX_KINETICS_DATA_EXPORT_RELATIVE_ERROR:
        raise ValueError("KineticsData output fit maximum relative error "
                         f"{100 * max_error:.2f}% exceeds 25% at {temperatures[max_index]:g} K")
    record = kinetics.comment.strip()
    fit.comment = (f"Output-only modified Arrhenius fit ({len(fit.arrhenius) if isinstance(fit, MultiArrhenius) else 1} terms) "
                   + f"over the source table range {temperatures[0]:g}-{temperatures[-1]:g} K; "
                   + f"maximum relative fit error {100 * max_error:.2f}% at {temperatures[max_index]:g} K "
                   + "over every Tdata point. The source table is the record."
                   + ("\n" + record if record else ""))
    return fit
