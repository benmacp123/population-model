"""Tests for the version 2 rolling-origin scoring helpers."""

import numpy as np
import pandas as pd
import pytest

from evaluate_state_space import (
    CONFIRMATION_SPECIES,
    DEVELOPMENT_SPECIES,
    _naive_forecasts,
    crps_normal,
    lagged_rain,
)


def test_crps_of_normal_matches_known_value():
    """CRPS at the forecast mean is sd * (sqrt(2) - 1) / sqrt(pi)."""
    value = crps_normal(np.array(0.0), np.array(2.0), np.array(0.0))
    assert value == pytest.approx(2.0 * (np.sqrt(2) - 1) / np.sqrt(np.pi))


def test_crps_of_point_forecast_is_absolute_error():
    """A zero-sd forecast scores its absolute error, not NaN."""
    value = crps_normal(np.array([0.0]), np.array([0.0]), np.array([1.5]))
    assert value[0] == pytest.approx(1.5)


def test_crps_penalises_distance_and_overconfidence():
    """A far-off or overconfident forecast scores worse."""
    near = crps_normal(np.array(0.0), np.array(1.0), np.array(0.5))
    far = crps_normal(np.array(0.0), np.array(1.0), np.array(3.0))
    tight_and_wrong = crps_normal(np.array(0.0), np.array(0.05), np.array(0.5))
    assert far > near
    assert tight_and_wrong > near


def test_lagged_rain_uses_previous_season_and_early_standardisation():
    """Lag 1 maps year t to seasonal year t-1, scaled on 1962-1979 only."""
    rainfall = pd.DataFrame(
        {"seasonal_year": range(1962, 1998), "dry_season_mm": range(36)}
    )
    years = np.array([1980, 1981])
    rain = lagged_rain(rainfall, years, "dry_season_mm", 1)
    early = np.arange(18)
    expected = (np.array([17, 18]) - early.mean()) / early.std(ddof=1)
    np.testing.assert_allclose(rain, expected)


def test_lagged_rain_rejects_missing_seasons():
    """A model year without its lagged season is an error."""
    rainfall = pd.DataFrame({"seasonal_year": [1970], "dry_season_mm": [1.0]})
    with pytest.raises(ValueError, match="No rainfall"):
        lagged_rain(rainfall, np.array([1965]), "dry_season_mm", 1)


def test_random_walk_carries_last_observation_and_widens_over_gaps():
    """The random walk forecasts the last value with sd growing with the gap."""
    years = np.arange(2000, 2005)
    y = np.array([[1.0, 2.0], [2.0, 3.0], [1.0, np.nan], [2.0, np.nan], [9.0, 9.0]])
    forecasts = _naive_forecasts(y, years, origin_index=3)
    mean, sd = forecasts["random_walk"]
    assert mean.tolist() == [2.0, 3.0]
    assert sd[1] == pytest.approx(sd[0] * np.sqrt(3))
    level_mean, _ = forecasts["calibration_mean"]
    assert level_mean.tolist() == [1.5, 2.5]


def test_species_sets_do_not_overlap():
    """Confirmation species must be unseen during development."""
    assert not set(DEVELOPMENT_SPECIES) & set(CONFIRMATION_SPECIES)
    assert "Buffalo" not in CONFIRMATION_SPECIES
