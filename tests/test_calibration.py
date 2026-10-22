"""Tests for bounded calibration and holdout forecasting."""

import pandas as pd
import pytest

from calibration import calibrate_species, holdout_forecast, near_best
from growth_model import population_change


@pytest.fixture
def traits():
    """Minimal trait record with a one-year interbirth interval."""
    return {
        "name": "TestBuck",
        "avg_weight_kg": 100.0,
        "gestation_period_months": 12,
        "interbirth_interval_days": 365.25,
        "life_span_years": 10,
        "sex_matur_m_years": 2,
        "sex_matur_f_years": 2,
        "male_ratio": 0.5,
    }


@pytest.fixture
def counts():
    """Short count series with one unreported year and a 50/50 seed split."""
    totals = {1965: 1000.0, 1967: 1100.0, 1969: 1150.0, 1970: None, 1972: 1200.0}
    rows = [
        {
            "name": "TestBuck",
            "year": year,
            "total_pop": total,
            "adult_pop": None if total is None else total / 2,
            "juvenile_pop": None if total is None else total / 2,
        }
        for year, total in totals.items()
    ]
    return pd.DataFrame(rows)


GRID = {
    "juv_death_rate": (0.2, 0.4),
    "birth_multiplier": (0.5, 1.0),
    "capacity_multiplier": (1.0,),
}


def test_birth_multiplier_scales_births(traits):
    """Births scale linearly with the birth multiplier."""
    full = population_change(traits, adult_pop=100.0, juvenile_pop=0.0)
    half = population_change(
        dict(traits, birth_multiplier=0.5), adult_pop=100.0, juvenile_pop=0.0
    )
    # With no juveniles to mature or die, juveniles next year are births only.
    assert half["juvenile_pop"] == pytest.approx(full["juvenile_pop"] / 2)


def test_maturation_is_capped_at_surviving_juveniles(traits):
    """Maturation cannot exceed the juveniles that survive the year."""
    record = dict(traits, juv_death_rate=0.9, sex_matur_m_years=1, sex_matur_f_years=1)
    result = population_change(record, adult_pop=0.0, juvenile_pop=100.0)
    # Only the 10 survivors can mature; without the cap 100 would.
    assert result["adult_pop"] == pytest.approx(10.0)
    assert result["juvenile_pop"] == pytest.approx(0.0)


def test_calibration_scores_every_grid_point_best_first(traits, counts):
    """Every grid point is scored and results are sorted best first."""
    fits = calibrate_species(traits, counts, "TestBuck", 1965, 1969, GRID)
    assert len(fits) == 4
    assert fits["calibration_mae"].is_monotonic_increasing


def test_calibration_ignores_years_after_cutoff(traits, counts):
    """Counts after the calibration end cannot influence the fit."""
    changed = counts.copy()
    changed.loc[changed["year"] == 1972, "total_pop"] = 99999.0
    original = calibrate_species(traits, counts, "TestBuck", 1965, 1969, GRID)
    altered = calibrate_species(traits, changed, "TestBuck", 1965, 1969, GRID)
    pd.testing.assert_frame_equal(original, altered)


def test_near_best_keeps_sets_within_tolerance():
    """The scenario set keeps only fits within the error tolerance."""
    fits = pd.DataFrame({"calibration_mae": [10.0, 10.5, 11.0, 12.0]})
    assert len(near_best(fits, tolerance=0.1)) == 3


def test_holdout_forecast_starts_from_anchor_and_brackets_central(traits, counts):
    """The forecast restarts at the anchor total and low <= central <= high."""
    fits = near_best(
        calibrate_species(traits, counts, "TestBuck", 1965, 1969, GRID), 1.0
    )
    forecast = holdout_forecast(traits, counts, "TestBuck", fits, 1969, 1972, 1969)
    assert forecast.index.tolist() == [1969, 1970, 1971, 1972]
    assert forecast.loc[1969, "central"] == pytest.approx(1150.0)
    assert (forecast["low"] <= forecast["central"]).all()
    assert (forecast["central"] <= forecast["high"]).all()


def test_holdout_forecast_rejects_unreported_anchor(traits, counts):
    """Restarting from a year with no reported total is an error."""
    fits = calibrate_species(traits, counts, "TestBuck", 1965, 1969, GRID)
    with pytest.raises(ValueError):
        holdout_forecast(traits, counts, "TestBuck", fits, 1970, 1972, 1969)


def test_holdout_forecast_carries_model_stage_share(traits, counts):
    """Protects the restart from re-imposing the seeded 50/50 split."""
    fits = calibrate_species(traits, counts, "TestBuck", 1965, 1969, GRID).head(1)
    seeded = holdout_forecast(traits, counts, "TestBuck", fits, 1969, 1972, 1969)
    carried = holdout_forecast(
        traits, counts, "TestBuck", fits, 1969, 1972, 1969, stage_start_year=1965
    )
    # Both restart from the observed total; only the stage split differs.
    assert carried.loc[1969, "central"] == pytest.approx(1150.0)
    assert carried.loc[1970, "central"] != pytest.approx(seeded.loc[1970, "central"])
