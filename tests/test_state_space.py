"""Tests for the pooled Gompertz state-space model."""

import numpy as np
import pandas as pd
import pytest

from state_space import Params, SpeciesSeries, fit, forecast_next, kalman_filter
from state_space import prepare_series


def _params(**overrides):
    """Parameters for one region with round numbers."""
    values = {
        "m": np.array([5.0]),
        "b": 0.5,
        "c": 0.0,
        "q": 0.2,
        "o": np.array([0.1, 0.1]),
    }
    values.update(overrides)
    return Params(**values)


def test_filter_matches_hand_worked_two_steps():
    """Log-likelihood and update agree with the Kalman equations by hand."""
    params = _params()
    y = np.array([[5.1], [5.3]])
    loglik, d, p = kalman_filter(params, y, np.zeros(2), np.zeros(2, dtype=int))

    # Year 0: stationary prior var q^2 / (1 - b^2).
    p0 = 0.04 / 0.75
    s0 = p0 + 0.01
    expected = -0.5 * (np.log(2 * np.pi * s0) + 0.1**2 / s0)
    d0, p0 = (p0 / s0) * 0.1, p0 * (1 - p0 / s0)
    # Year 1: predict, then update on the innovation.
    d1, p1 = 0.5 * d0, 0.25 * p0 + 0.04
    s1 = p1 + 0.01
    innovation = 0.3 - d1
    expected -= 0.5 * (np.log(2 * np.pi * s1) + innovation**2 / s1)

    assert loglik == pytest.approx(expected)
    assert d[0] == pytest.approx(d1 + (p1 / s1) * innovation)
    assert p[0] == pytest.approx(p1 * (1 - p1 / s1))


def test_filter_skips_update_in_unobserved_years():
    """A missing year only predicts, so it adds nothing to the likelihood."""
    params = _params()
    y_full = np.array([[5.1], [np.nan]])
    loglik_full, _, p = kalman_filter(
        params, y_full, np.zeros(2), np.zeros(2, dtype=int)
    )
    loglik_one, _, p_one = kalman_filter(
        params, y_full[:1], np.zeros(1), np.zeros(1, dtype=int)
    )
    assert loglik_full == pytest.approx(loglik_one)
    assert p[0] == pytest.approx(0.25 * p_one[0] + 0.04)


def _simulate(rng, b=0.6, c=0.25, q=0.15, o=0.1, n_years=33, n_regions=4):
    """Simulated pooled series with known parameters and some gaps."""
    m = np.array([7.0, 8.0, 6.5, 6.0])[:n_regions]
    rain = rng.normal(size=n_years)
    d = rng.normal(scale=q / np.sqrt(1 - b**2), size=n_regions)
    y = np.empty((n_years, n_regions))
    for t in range(n_years):
        if t > 0:
            d = b * d + c * rain[t] + rng.normal(scale=q, size=n_regions)
        y[t] = m + d + rng.normal(scale=o, size=n_regions)
    y[[1, 2, 3, 11, 13]] = np.nan
    series = SpeciesSeries(
        years=np.arange(1965, 1965 + n_years),
        regions=[f"R{i}" for i in range(n_regions)],
        y=y,
        era=(np.arange(n_years) > 11).astype(int),
    )
    return series, rain


def test_fit_recovers_simulated_parameters():
    """Pooled fitting recovers b and c from series of the real length."""
    estimates = []
    for seed in range(5):
        series, rain = _simulate(np.random.default_rng(seed))
        params = fit(series, rain)
        estimates.append((params.b, params.c))
    b_hat, c_hat = np.mean(estimates, axis=0)
    assert b_hat == pytest.approx(0.6, abs=0.15)
    assert c_hat == pytest.approx(0.25, abs=0.08)


def test_ablation_fixes_rain_coefficient_at_zero():
    """The no-rainfall model never uses the covariate."""
    series, rain = _simulate(np.random.default_rng(0))
    assert fit(series, rain, use_rain=False).c == 0.0


def test_forecast_ignores_years_after_origin():
    """Changing observations after the origin cannot change the forecast."""
    series, rain = _simulate(np.random.default_rng(1))
    params = fit(series, rain)
    altered = SpeciesSeries(series.years, series.regions, series.y.copy(), series.era)
    altered.y[21:] += 5.0
    mean, sd = forecast_next(params, series, rain, 20)
    mean_alt, sd_alt = forecast_next(params, altered, rain, 20)
    np.testing.assert_allclose(mean, mean_alt)
    np.testing.assert_allclose(sd, sd_alt)


def test_prepare_series_drops_sparse_regions_and_logs_counts():
    """Regions below the observation threshold are excluded."""
    regional = pd.DataFrame(
        {
            "name": "Kudu",
            "year": list(range(1965, 1975)) + [1970],
            "region": ["Central"] * 10 + ["South"],
            "count": [99] * 10 + [5],
        }
    )
    series = prepare_series(regional, "Kudu", 1965, 1976, min_obs=8)
    assert series.regions == ["Central"]
    assert series.y[0, 0] == pytest.approx(np.log(100.0))
    assert np.isnan(series.y[-1, 0])
    assert series.era.tolist() == [0] * 12


def test_shared_observation_sd_is_equal_across_eras():
    """The shared option fits one survey sd; the default fits one per era."""
    series, rain = _simulate(np.random.default_rng(2))
    shared = fit(series, rain, shared_o=True)
    separate = fit(series, rain)
    assert shared.o[0] == shared.o[1]
    assert separate.o.shape == (2,)
