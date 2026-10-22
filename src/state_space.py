"""Pooled Gompertz state-space model for regional census counts.

For one species, with ``x`` the true log abundance in region ``r`` and year
``t`` and ``y = log(count + 1)`` the observed value:

.. code-block:: text

    x[r,t] - m[r] = b * (x[r,t-1] - m[r]) + c * rain[t] + w,   w ~ N(0, q^2)
    y[r,t]        = x[r,t] + v,                                v ~ N(0, o[era]^2)

``rain[t]`` is the covariate acting on the transition into year ``t``;
callers decide its lag. ``m`` is region-specific; ``b``, ``c`` and ``q`` are
shared across regions, which is how the regional series pool information.
There is one observation sd per survey era (optionally shared). With
series this short, process and observation noise are not separately
identifiable: one of the two sds often sits on its lower bound. One-step
forecasts depend on their combined variance, which is identified, so
forecasts are scored but the fitted parameters are not interpreted. The Kalman filter skips the
update step in unobserved years, so missing counts need no imputation.
Fitting is by maximum likelihood; this is an engineering-style state
estimator, not a Bayesian hierarchical model.
"""

import warnings
from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import minimize

B_BOUNDS = (-0.95, 0.99)
LOG_SD_BOUNDS = (-6.0, 1.5)


@dataclass
class SpeciesSeries:
    """Aligned regional observations for one species.

    Attributes:
        years: Consecutive years covered, shape ``(T,)``.
        regions: Region names, length ``R``.
        y: Log counts, shape ``(T, R)``, NaN where unobserved.
        era: Survey-era index (0 or 1) for each year, shape ``(T,)``.
    """

    years: np.ndarray
    regions: Sequence[str]
    y: np.ndarray
    era: np.ndarray


@dataclass
class Params:
    """Model parameters for one species.

    Attributes:
        m: Region mean log levels, shape ``(R,)``.
        b: Density-dependence coefficient (1 is a random walk).
        c: Rainfall coefficient on the log scale.
        q: Process-noise sd.
        o: Observation-noise sd per survey era, shape ``(2,)``.
    """

    m: np.ndarray
    b: float
    c: float
    q: float
    o: np.ndarray


def prepare_series(
    regional: pd.DataFrame,
    species: str,
    first_year: int = 1965,
    last_year: int = 1997,
    last_corrected_year: int = 1976,
    min_obs: int = 8,
) -> SpeciesSeries:
    """Build the year-by-region log-count matrix for one species.

    Args:
        regional: Long table from ``data/raw/regional_counts.csv``.
        species: Species name as in the ``name`` column.
        first_year: First year of the matrix.
        last_year: Last year of the matrix.
        last_corrected_year: Last year of the first survey era.
        min_obs: Regions with fewer observations are dropped.

    Returns:
        The aligned series.

    Raises:
        ValueError: If no region has enough observations.
    """
    rows = regional[regional["name"] == species]
    wide = rows.pivot_table(index="year", columns="region", values="count")
    years = np.arange(first_year, last_year + 1)
    wide = wide.reindex(years)
    keep = [region for region in wide.columns if wide[region].notna().sum() >= min_obs]
    if not keep:
        raise ValueError(f"No region has {min_obs} observations for '{species}'.")
    return SpeciesSeries(
        years=years,
        regions=keep,
        y=np.log(wide[keep].to_numpy(dtype=float) + 1.0),
        era=(years > last_corrected_year).astype(int),
    )


def kalman_filter(params: Params, y: np.ndarray, rain: np.ndarray, era: np.ndarray):
    """Run the filter and return the log-likelihood and final state.

    The first year starts from the stationary distribution of the centred
    state, which exists because ``b`` is bounded below 1.

    Args:
        params: Model parameters.
        y: Log counts, shape ``(T, R)``, NaN where unobserved.
        rain: Covariate for the transition into each year, shape ``(T,)``.
        era: Survey-era index per year, shape ``(T,)``.

    Returns:
        ``(loglik, d, p)``: total log-likelihood, and the filtered mean and
        variance of the centred state ``x - m`` in the last year, each of
        shape ``(R,)``.
    """
    n_regions = y.shape[1]
    d = np.zeros(n_regions)
    p = np.full(n_regions, params.q**2 / (1.0 - params.b**2))
    loglik = 0.0
    for t in range(y.shape[0]):
        if t > 0:
            d = params.b * d + params.c * rain[t]
            p = params.b**2 * p + params.q**2
        seen = ~np.isnan(y[t])
        if not seen.any():
            continue
        s = p[seen] + params.o[era[t]] ** 2
        innovation = y[t, seen] - params.m[seen] - d[seen]
        loglik -= 0.5 * float(np.sum(np.log(2.0 * np.pi * s) + innovation**2 / s))
        gain = p[seen] / s
        d[seen] = d[seen] + gain * innovation
        p[seen] = (1.0 - gain) * p[seen]
    return loglik, d, p


def _unpack(
    theta: np.ndarray, n_regions: int, use_rain: bool, shared_o: bool
) -> Params:
    """Map an unconstrained optimiser vector to model parameters.

    Args:
        theta: Optimiser vector: region means, ``b`` on a logit scale,
            ``c`` (if used), then log ``q`` and one or two log ``o``.
        n_regions: Number of regions.
        use_rain: Whether ``theta`` includes ``c``.
        shared_o: Whether one log ``o`` serves both eras.

    Returns:
        The model parameters.
    """
    low, high = B_BOUNDS
    b = low + (high - low) / (1.0 + np.exp(-theta[n_regions]))
    rest = theta[n_regions + 1 :]
    c = rest[0] if use_rain else 0.0
    rest = rest[1:] if use_rain else rest
    return Params(
        m=theta[:n_regions].copy(),
        b=float(b),
        c=float(c),
        q=float(np.exp(rest[0])),
        o=np.exp(np.repeat(rest[1], 2) if shared_o else rest[1:3]),
    )


def _logit_b(b: float) -> float:
    """Inverse of the bounded transform for ``b``.

    Args:
        b: Value strictly inside ``B_BOUNDS``.

    Returns:
        The unconstrained value.
    """
    low, high = B_BOUNDS
    share = (b - low) / (high - low)
    return float(np.log(share / (1.0 - share)))


def fit(
    series: SpeciesSeries,
    rain: np.ndarray,
    use_rain: bool = True,
    fixed_o: Optional[float] = None,
    shared_o: bool = False,
) -> Params:
    """Fit the model by maximum likelihood from several starting points.

    Args:
        series: Observations, already masked to the fitting window.
        rain: Covariate for the transition into each year, shape ``(T,)``.
        use_rain: If false, ``c`` is fixed at zero (the ablation model).
        fixed_o: If given, both observation sds are fixed at this value
            (sensitivity runs only).
        shared_o: If true, both survey eras share one observation sd
            (tried during development, not used in the frozen design).

    Returns:
        The best parameters found.
    """
    n_regions = len(series.regions)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        m0 = np.nanmean(series.y, axis=0)
    # A region with no observations in the window starts at the overall mean.
    m0 = np.where(np.isnan(m0), np.nanmean(series.y), m0)
    bounds = [(None, None)] * n_regions + [(-8.0, 8.0)]
    bounds += [(-3.0, 3.0)] if use_rain else []
    bounds += [LOG_SD_BOUNDS] * (2 if shared_o else 3)

    def negative_loglik(theta):
        """Objective for the optimiser."""
        params = _unpack(theta, n_regions, use_rain, shared_o)
        if fixed_o is not None:
            params.o = np.array([fixed_o, fixed_o])
        return -kalman_filter(params, series.y, rain, series.era)[0]

    best = None
    for b0 in (0.3, 0.7, 0.95):
        for q0, o0 in ((0.2, 0.1), (0.1, 0.3)):
            sds = [q0, o0] if shared_o else [q0, o0, o0]
            theta = np.concatenate(
                [m0, [_logit_b(b0)], [0.0] if use_rain else [], np.log(sds)]
            )
            result = minimize(negative_loglik, theta, method="L-BFGS-B", bounds=bounds)
            if best is None or result.fun < best.fun:
                best = result
    params = _unpack(best.x, n_regions, use_rain, shared_o)
    if fixed_o is not None:
        params.o = np.array([fixed_o, fixed_o])
    return params


def forecast_next(
    params: Params, series: SpeciesSeries, rain: np.ndarray, origin_index: int
):
    """One-step predictive distribution for the year after the origin.

    Only observations up to and including ``origin_index`` are used.

    Args:
        params: Parameters fitted on data up to the origin.
        series: Full observations; later years are ignored.
        rain: Covariate for the transition into each year, shape ``(T,)``.
        origin_index: Index of the origin year in ``series.years``.

    Returns:
        ``(mean, sd)`` of the predicted log count for each region in year
        ``origin_index + 1``, each of shape ``(R,)``.
    """
    window = slice(0, origin_index + 1)
    _, d, p = kalman_filter(params, series.y[window], rain[window], series.era[window])
    t_next = origin_index + 1
    d_next = params.b * d + params.c * rain[t_next]
    p_next = params.b**2 * p + params.q**2
    sd = np.sqrt(p_next + params.o[series.era[t_next]] ** 2)
    return params.m + d_next, sd
