"""Tests for the version 3 survey-design and alert module."""

import numpy as np
import pytest

from monitoring import alert_metrics, core_regions, crossing_probability
from monitoring import filter_paths, reliability, simulate, survey_years
from monitoring import with_noise_split
from state_space import Params, SpeciesSeries, fit, kalman_filter


def _params(n_regions=2, b=0.7, q=0.15, o=0.2):
    """No-rain parameters with round numbers."""
    return Params(m=np.full(n_regions, 5.0), b=b, c=0.0, q=q, o=np.array([o, o]))


def test_vectorised_filter_matches_state_space_filter_with_gaps():
    """Every path's final filtered state equals the frozen version 2 filter."""
    params = _params()
    rng = np.random.default_rng(1)
    _, observed = simulate(params, 3, rng)
    observed[:, ~survey_years(observed.shape[1], 2)] = np.nan
    d, p = filter_paths(params, observed)
    n_total = observed.shape[1]
    for path in range(3):
        _, d_ref, p_ref = kalman_filter(
            params, observed[path], np.zeros(n_total), np.ones(n_total, dtype=int)
        )
        np.testing.assert_allclose(d[path, -1], d_ref)
        np.testing.assert_allclose(p[path, -1], p_ref)


def test_survey_years_always_include_scenario_onset():
    """The first scenario year is surveyed at every interval."""
    for interval in (1, 2, 3):
        mask = survey_years(25, interval, onset=10)
        assert mask[10]
        assert mask.sum() == len(range(10 % interval, 25, interval))


def test_filter_rejects_rain_model():
    """The vectorised filter refuses parameters with a rainfall term."""
    params = _params()
    params.c = 0.1
    with pytest.raises(ValueError):
        filter_paths(params, np.zeros((1, 3, 2)))


def test_decline_lowers_true_state_but_drift_only_observations():
    """A decline moves the truth; detection drift moves only the counts."""
    params = _params(q=1e-6, o=1e-6)
    true_d, obs_d = simulate(params, 1, np.random.default_rng(0), decline=0.1)
    true_b, obs_b = simulate(params, 1, np.random.default_rng(0), drift=0.1)
    expected = 15 * np.log(0.9)
    assert true_d[0, -1, 0] == pytest.approx(expected, abs=1e-4)
    assert obs_d[0, -1, 0] - 5.0 == pytest.approx(expected, abs=1e-4)
    assert true_b[0, -1, 0] == pytest.approx(0.0, abs=1e-4)
    assert obs_b[0, -1, 0] - 5.0 == pytest.approx(expected, abs=1e-4)


def test_noise_split_keeps_fitted_total_variance():
    """Changing the share moves variance between q and o, not the total."""
    fitted = _params(q=0.3, o=0.4)
    split = with_noise_split(fitted, 0.25, obs_sd=0.1)
    assert split.q**2 == pytest.approx(0.25 * 0.25)
    np.testing.assert_allclose(split.o, [0.1, 0.1])


def test_fit_recovers_total_variance_from_simulation():
    """The identified quantity q^2 + o^2 is recovered from a long series."""
    params = _params(n_regions=4, b=0.6, q=0.15, o=0.2)
    _, observed = simulate(params, 1, np.random.default_rng(3), n_years=140)
    n_total = observed.shape[1]
    series = SpeciesSeries(
        years=np.arange(n_total),
        regions=list("ABCD"),
        y=observed[0],
        era=np.ones(n_total, dtype=int),
    )
    fitted = fit(series, np.zeros(n_total), use_rain=False)
    total = fitted.q**2 + fitted.o[1] ** 2
    assert total == pytest.approx(0.15**2 + 0.2**2, rel=0.2)


def test_crossing_probability_is_calibrated_when_model_is_true():
    """With no decline, mean probability matches the frequency below threshold."""
    params = _params(n_regions=2)
    true_centred, observed = simulate(params, 4000, np.random.default_rng(5))
    d, p = filter_paths(params, observed)
    log_threshold = np.log(0.85)
    prob = crossing_probability(d, p, [0, 1], log_threshold)
    below = true_centred.mean(axis=-1) < log_threshold
    assert prob.mean() == pytest.approx(below.mean(), abs=0.01)
    table = reliability(prob, below)
    big = table[table["n"] >= 2000]
    np.testing.assert_allclose(big["observed_freq"], big["mean_prob"], atol=0.05)


def test_crossing_probability_averages_core_regions():
    """The index ignores non-core regions and shrinks variance by region count."""
    d = np.array([[0.0, -0.2, 9.0]])
    p = np.array([[0.04, 0.04, 9.0]])
    prob = crossing_probability(d, p, [0, 1], log_threshold=-0.1)
    assert prob[0] == pytest.approx(0.5)


def test_alert_metrics_delay_and_false_alarms():
    """Delay counts from the true crossing; false alarms only after onset."""
    prob = np.array([[0.9, 0.0, 0.0, 0.6, 0.9], [0.0, 0.0, 0.0, 0.0, 0.0]])
    true_index = np.array([[0.0, 0.0, -1.0, -1.0, -1.0], [0.0, 0.0, -1.0, 0.0, 0.0]])
    surveyed = np.ones(5, dtype=bool)
    declining = alert_metrics(
        prob, true_index, surveyed, 0.5, np.log(0.7), True, onset=1
    )
    assert declining["n_crossed"] == 2
    assert declining["median_delay_years"] == 1.0
    assert declining["miss_rate"] == 0.5
    stable = alert_metrics(prob, true_index, surveyed, 0.5, np.log(0.7), False, onset=1)
    assert stable["false_alarm_per_survey_year"] == pytest.approx(2 / 8)
    assert stable["any_false_alarm_10y"] == 0.5


def test_core_regions_use_baseline_mean_counts():
    """Regions below the mean-count cut are left out of the index."""
    counts = np.array([[100.0, 10.0], [60.0, 0.0], [0.0, 500.0]])
    series = SpeciesSeries(
        years=np.array([1970, 1980, 1990]),
        regions=["A", "B"],
        y=np.log(counts + 1.0),
        era=np.zeros(3, dtype=int),
    )
    assert core_regions(series, baseline_end=1985) == [0]
    with pytest.raises(ValueError):
        core_regions(series, baseline_end=1985, min_mean_count=1000)
