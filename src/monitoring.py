"""Simulated survey designs and threshold alerts for version 3.

Version 3 asks how well a decline alert built on the version 2 Kalman
filter works, not whether counts can be forecast. This module simulates
regional log counts from a no-rain ``Params`` operating model, thins them
to a survey interval, filters every simulated path at once and turns the
filtered state into the probability that a core-region index has fallen
below a threshold. The pre-registered design is summarised in the README.

The vectorised filter repeats the update equations of
``state_space.kalman_filter`` for the no-rain model (``c = 0``) and keeps
the filtered state for every year rather than only the last; a test checks
the two agree.
"""

from dataclasses import replace
from typing import Sequence

import numpy as np
import pandas as pd
from scipy.stats import norm

from state_space import Params, SpeciesSeries

BURN_IN_YEARS = 10
SCENARIO_YEARS = 15
MISS_WINDOW_YEARS = 5
FALSE_ALARM_WINDOW_YEARS = 10


def core_regions(
    series: SpeciesSeries, baseline_end: int, min_mean_count: float = 50.0
) -> list:
    """Indices of regions holding enough animals for a percentage threshold.

    Args:
        series: Regional observations for one species.
        baseline_end: Last year of the baseline window.
        min_mean_count: Smallest mean count over the baseline window.

    Returns:
        Column indices into ``series.regions``.

    Raises:
        ValueError: If no region reaches ``min_mean_count``.
    """
    counts = np.exp(series.y[series.years <= baseline_end]) - 1.0
    means = np.nanmean(counts, axis=0)
    core = [int(i) for i in np.flatnonzero(means >= min_mean_count)]
    if not core:
        raise ValueError(f"No region averages {min_mean_count} animals.")
    return core


def with_noise_split(params: Params, share: float, obs_sd: float) -> Params:
    """Operating model with a chosen process share of the fitted variance.

    Process and observation noise are not separately identified, so the
    fitted total ``q^2 + o^2`` (survey era 2) is kept and divided between
    them. The survey being designed then gets its own observation sd.

    Args:
        params: Fitted no-rain parameters.
        share: Fraction of the fitted total variance that is process noise.
        obs_sd: Observation sd of the simulated survey, on the log scale.

    Returns:
        A copy of ``params`` with the new ``q`` and both ``o`` set to
        ``obs_sd``.
    """
    total = params.q**2 + params.o[1] ** 2
    return replace(
        params,
        q=float(np.sqrt(share * total)),
        o=np.array([obs_sd, obs_sd]),
    )


def simulate(
    params: Params,
    n_sims: int,
    rng: np.random.Generator,
    decline: float = 0.0,
    drift: float = 0.0,
    burn_in: int = BURN_IN_YEARS,
    n_years: int = SCENARIO_YEARS,
):
    """Simulate true and observed log counts for every region.

    The centred state starts from its stationary distribution and follows
    the Gompertz process. A true decline lowers the level by
    ``log(1 - decline)`` each scenario year; a detection drift lowers only
    the observations by ``log(1 - drift)`` each scenario year.

    Args:
        params: Operating model (``c`` is ignored).
        n_sims: Number of simulated paths.
        rng: Random generator.
        decline: Annual fractional decline of the true population.
        drift: Annual fractional decline in detection.
        burn_in: Stable years before the scenario starts.
        n_years: Scenario years.

    Returns:
        ``(true_centred, observed)``: the true log abundance minus ``m``
        and the observed log counts, each of shape ``(n_sims, T, R)`` with
        ``T = burn_in + n_years``.
    """
    n_regions = len(params.m)
    n_total = burn_in + n_years
    scenario_year = np.clip(np.arange(n_total) - burn_in + 1, 0, None)
    level = scenario_year * np.log1p(-decline)
    bias = scenario_year * np.log1p(-drift)

    d = np.empty((n_sims, n_total, n_regions))
    stationary_sd = params.q / np.sqrt(1.0 - params.b**2)
    d[:, 0] = rng.normal(0.0, stationary_sd, (n_sims, n_regions))
    shocks = rng.normal(0.0, params.q, (n_sims, n_total, n_regions))
    for t in range(1, n_total):
        d[:, t] = params.b * d[:, t - 1] + shocks[:, t]

    true_centred = d + level[None, :, None]
    noise = rng.normal(0.0, params.o[1], (n_sims, n_total, n_regions))
    observed = params.m + true_centred + bias[None, :, None] + noise
    return true_centred, observed


def survey_years(n_total: int, interval: int, onset: int = BURN_IN_YEARS):
    """Mask of years in which a survey is flown.

    Surveys are aligned so the scenario's first year is always surveyed,
    whatever the interval.

    Args:
        n_total: Number of years.
        interval: Years between surveys.
        onset: Index of the first scenario year.

    Returns:
        Boolean array of shape ``(n_total,)``.
    """
    return (np.arange(n_total) - onset) % interval == 0


def filter_paths(params: Params, observed: np.ndarray):
    """Filter every simulated path and keep each year's filtered state.

    Args:
        params: Model the manager uses; must have ``c = 0``. Observations
            use ``o[1]`` throughout.
        observed: Log counts, shape ``(n_sims, T, R)``, NaN where not
            surveyed.

    Returns:
        ``(d, p)``: filtered mean and variance of the centred state, each
        of shape ``(n_sims, T, R)``.

    Raises:
        ValueError: If ``params.c`` is not zero.
    """
    if params.c != 0.0:
        raise ValueError("filter_paths supports the no-rain model only.")
    n_sims, n_total, n_regions = observed.shape
    d = np.zeros((n_sims, n_regions))
    p = np.full((n_sims, n_regions), params.q**2 / (1.0 - params.b**2))
    obs_var = params.o[1] ** 2
    d_out = np.empty_like(observed)
    p_out = np.empty_like(observed)
    for t in range(n_total):
        if t > 0:
            d = params.b * d
            p = params.b**2 * p + params.q**2
        y_t = observed[:, t]
        seen = ~np.isnan(y_t)
        gain = p / (p + obs_var)
        innovation = np.where(seen, y_t - params.m - d, 0.0)
        d = d + np.where(seen, gain * innovation, 0.0)
        p = np.where(seen, (1.0 - gain) * p, p)
        d_out[:, t] = d
        p_out[:, t] = p
    return d_out, p_out


def crossing_probability(
    d: np.ndarray, p: np.ndarray, core: Sequence[int], log_threshold: float
) -> np.ndarray:
    """Probability that the core-region index is below a threshold.

    The index is the mean of the centred log state over the core regions.
    Regional states are independent in the model, so its variance is the
    sum of the regional variances over the number of regions squared.

    Args:
        d: Filtered means, shape ``(..., R)``.
        p: Filtered variances, same shape.
        core: Region indices forming the index.
        log_threshold: Threshold on the index, e.g. ``log(0.7)``.

    Returns:
        Probabilities with the region axis removed.
    """
    core = list(core)
    mean = d[..., core].mean(axis=-1)
    sd = np.sqrt(p[..., core].sum(axis=-1)) / len(core)
    return norm.cdf((log_threshold - mean) / sd)


def _first_true(mask: np.ndarray) -> np.ndarray:
    """Index of the first true value along the last axis.

    Args:
        mask: Boolean array.

    Returns:
        Integer indices, -1 where no value is true.
    """
    first = mask.argmax(axis=-1)
    return np.where(mask.any(axis=-1), first, -1)


def alert_metrics(
    prob: np.ndarray,
    true_index: np.ndarray,
    surveyed: np.ndarray,
    alpha: float,
    log_threshold: float,
    declining: bool,
    onset: int = BURN_IN_YEARS,
) -> dict:
    """Pre-registered alert metrics for one design cell.

    Alerts are raised only in surveyed years from the scenario onset on.
    False-alarm metrics apply to runs without a true decline (including
    detection drift); delay and miss rate apply to declining runs.

    Args:
        prob: Crossing probabilities, shape ``(n_sims, T)``.
        true_index: True core-region index, shape ``(n_sims, T)``.
        surveyed: Survey mask, shape ``(T,)``.
        alpha: Alert when ``prob > alpha``.
        log_threshold: Threshold on the index.
        declining: Whether the scenario has a true decline.
        onset: Index of the first scenario year.

    Returns:
        Dict of metric name to value; metrics that do not apply are NaN.
    """
    scenario = np.arange(prob.shape[1]) >= onset
    alerts = (prob > alpha) & surveyed & scenario
    nan = float("nan")
    if not declining:
        window = scenario & (
            np.arange(prob.shape[1]) < onset + FALSE_ALARM_WINDOW_YEARS
        )
        return {
            "false_alarm_per_survey_year": float(alerts[:, surveyed & scenario].mean()),
            "any_false_alarm_10y": float(alerts[:, window].any(axis=1).mean()),
            "median_delay_years": nan,
            "miss_rate": nan,
            "n_crossed": 0,
        }
    crossed = _first_true((true_index < log_threshold) & scenario)
    first_alert = _first_true(alerts)
    has_cross = crossed >= 0
    delay = np.where(first_alert >= 0, first_alert - crossed, np.inf)[has_cross]
    finite = delay[np.isfinite(delay)]
    return {
        "false_alarm_per_survey_year": nan,
        "any_false_alarm_10y": nan,
        "median_delay_years": float(np.median(finite)) if finite.size else nan,
        "miss_rate": float(np.mean(delay > MISS_WINDOW_YEARS)) if delay.size else nan,
        "n_crossed": int(has_cross.sum()),
    }


def reliability(prob: np.ndarray, below: np.ndarray, n_bins: int = 10):
    """Observed frequency of the event in bins of predicted probability.

    Args:
        prob: Predicted probabilities, any shape.
        below: Whether the event happened, same shape.
        n_bins: Number of equal-width bins on [0, 1].

    Returns:
        One row per non-empty bin with its mid-point, count, mean
        predicted probability and observed frequency.
    """
    prob, below = np.ravel(prob), np.ravel(below)
    bins = np.minimum((prob * n_bins).astype(int), n_bins - 1)
    rows = []
    for b in np.unique(bins):
        inside = bins == b
        rows.append(
            {
                "bin_mid": (b + 0.5) / n_bins,
                "n": int(inside.sum()),
                "mean_prob": float(prob[inside].mean()),
                "observed_freq": float(below[inside].mean()),
            }
        )
    return pd.DataFrame(rows)
