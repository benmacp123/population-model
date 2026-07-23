"""Tests for the two-stage population model."""

import pandas as pd
import pytest

from growth_model import (
    adult_pop_change,
    calc_carrying_capacity,
    compare_to_counts,
    forecast_population,
    get_species_counts,
    get_species_traits,
    juv_pop_change,
    population_change,
    run_growth_model,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def traits():
    """Minimal single-species trait record with round numbers."""
    return {
        "name": "TestBuck",
        "avg_weight_kg": 100.0,
        "gestation_period_months": 12,
        "life_span_years": 10,
        "sex_matur_m_years": 2,
        "sex_matur_f_years": 2,
        "male_ratio": 0.5,
        "juv_death_rate": None,
        "carrying_capacity": None,
    }


@pytest.fixture
def traits_df(traits):
    return pd.DataFrame([traits])


@pytest.fixture
def counts_df():
    return pd.DataFrame(
        [
            {"name": "TestBuck", "year": 2018, "adult_pop": 100.0, "juvenile_pop": 50.0},
            {"name": "TestBuck", "year": 2020, "adult_pop": 130.0, "juvenile_pop": 60.0},
        ]
    ).assign(total_pop=lambda df: df["adult_pop"] + df["juvenile_pop"])


# ---------------------------------------------------------------------------
# Component arithmetic
# ---------------------------------------------------------------------------


def test_adult_pop_change():
    assert (
        adult_pop_change(adult_pop=100, matured_juveniles=10, adult_death_rate=0.05)
        == 10 - 5
    )


def test_juv_pop_change():
    assert (
        juv_pop_change(
            juvenile_pop=50, births=20, juv_death_rate=0.1, matured_juveniles=5
        )
        == 20 - 5 - 5
    )


# ---------------------------------------------------------------------------
# Carrying capacity
# ---------------------------------------------------------------------------


def test_carrying_capacity_uses_explicit_value_when_supplied(traits, traits_df, counts_df):
    traits["carrying_capacity"] = 1234.0
    assert calc_carrying_capacity(traits, traits_df, counts_df, 2018) == 1234.0


def test_carrying_capacity_falls_back_to_biomass_proxy(traits, traits_df, counts_df):
    # Single species at equal weight: proxy reduces to that year's total count.
    capacity = calc_carrying_capacity(traits, traits_df, counts_df, 2018)
    assert capacity == pytest.approx(150.0)


def test_carrying_capacity_returns_zero_for_unobserved_year(traits, traits_df, counts_df):
    assert calc_carrying_capacity(traits, traits_df, counts_df, 1999) == 0.0


# ---------------------------------------------------------------------------
# Annual step
# ---------------------------------------------------------------------------


def test_population_change_returns_all_three_totals(traits):
    result = population_change(traits, adult_pop=100.0, juvenile_pop=50.0)
    assert set(result) == {"adult_pop", "juvenile_pop", "total_pop"}
    assert result["total_pop"] == pytest.approx(
        result["adult_pop"] + result["juvenile_pop"]
    )


def test_juvenile_mortality_defaults_to_adult_rate(traits):
    """With juv_death_rate unset, the fallback must equal the adult rate."""
    fallback = population_change(traits, adult_pop=100.0, juvenile_pop=50.0)

    explicit = dict(traits, juv_death_rate=1.0 / traits["life_span_years"])
    matched = population_change(explicit, adult_pop=100.0, juvenile_pop=50.0)

    assert fallback["juvenile_pop"] == pytest.approx(matched["juvenile_pop"])


def test_higher_juvenile_mortality_reduces_juveniles(traits):
    low = population_change(
        dict(traits, juv_death_rate=0.1), adult_pop=100.0, juvenile_pop=50.0
    )
    high = population_change(
        dict(traits, juv_death_rate=0.6), adult_pop=100.0, juvenile_pop=50.0
    )
    assert high["juvenile_pop"] < low["juvenile_pop"]


def test_density_dependence_halts_maturation_at_capacity(traits):
    """At or above capacity the maturation term is zeroed, so no juveniles graduate."""
    at_capacity = population_change(
        traits, adult_pop=100.0, juvenile_pop=50.0, carrying_capacity=150.0
    )
    adult_death_rate = 1.0 / traits["life_span_years"]
    expected_adults = 100.0 - (100.0 * adult_death_rate)
    assert at_capacity["adult_pop"] == pytest.approx(expected_adults)


def test_below_capacity_allows_maturation(traits):
    constrained = population_change(
        traits, adult_pop=100.0, juvenile_pop=50.0, carrying_capacity=150.0
    )
    roomy = population_change(
        traits, adult_pop=100.0, juvenile_pop=50.0, carrying_capacity=10_000.0
    )
    assert roomy["adult_pop"] > constrained["adult_pop"]


def test_populations_are_floored_at_zero(traits):
    doomed = dict(traits, life_span_years=0.5, juv_death_rate=5.0)
    result = population_change(doomed, adult_pop=10.0, juvenile_pop=10.0)
    assert result["adult_pop"] >= 0.0
    assert result["juvenile_pop"] >= 0.0


# ---------------------------------------------------------------------------
# Forecasting
# ---------------------------------------------------------------------------


def test_forecast_starts_from_observed_values(traits_df, counts_df):
    forecast = forecast_population(traits_df, counts_df, "TestBuck", 2018, 3)
    first = forecast.iloc[0]
    assert first["year"] == 2018
    assert first["adult_pop"] == 100.0
    assert first["juvenile_pop"] == 50.0


def test_forecast_returns_horizon_plus_start_year(traits_df, counts_df):
    forecast = forecast_population(traits_df, counts_df, "TestBuck", 2018, 3)
    assert len(forecast) == 4
    assert list(forecast["year"]) == [2018, 2019, 2020, 2021]


def test_forecast_rejects_unknown_species(traits_df, counts_df):
    with pytest.raises(ValueError, match="No traits found"):
        forecast_population(traits_df, counts_df, "Griffin", 2018, 3)


def test_forecast_rejects_start_year_without_observation(traits_df, counts_df):
    with pytest.raises(ValueError, match="No observed count"):
        forecast_population(traits_df, counts_df, "TestBuck", 1999, 3)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_comparison_is_exact_in_the_start_year(traits_df, counts_df):
    """The start year is initialised from observation, so error there must be zero."""
    forecast = forecast_population(traits_df, counts_df, "TestBuck", 2018, 3)
    comparison = compare_to_counts(forecast, counts_df, "TestBuck")
    start_row = comparison[comparison["year"] == 2018].iloc[0]
    assert start_row["error"] == pytest.approx(0.0)
    assert start_row["pct_error"] == pytest.approx(0.0)


def test_comparison_only_covers_overlapping_years(traits_df, counts_df):
    forecast = forecast_population(traits_df, counts_df, "TestBuck", 2018, 3)
    comparison = compare_to_counts(forecast, counts_df, "TestBuck")
    # 2019 and 2021 are projected but never observed.
    assert list(comparison["year"]) == [2018, 2020]


def test_comparison_error_columns_are_consistent(traits_df, counts_df):
    forecast = forecast_population(traits_df, counts_df, "TestBuck", 2018, 3)
    comparison = compare_to_counts(forecast, counts_df, "TestBuck")
    for _, row in comparison.iterrows():
        assert row["error"] == pytest.approx(
            row["projected_total"] - row["observed_total"]
        )
        assert row["pct_error"] == pytest.approx(
            100.0 * row["error"] / row["observed_total"]
        )


def test_comparison_is_empty_without_overlap(traits_df, counts_df):
    forecast = forecast_population(traits_df, counts_df, "TestBuck", 2018, 0)
    unrelated = counts_df.assign(year=counts_df["year"] + 100)
    assert compare_to_counts(forecast, unrelated, "TestBuck").empty


# ---------------------------------------------------------------------------
# Placeholder data and CLI
# ---------------------------------------------------------------------------


def test_placeholder_species_appear_in_both_tables():
    """Every species with counts must also have traits, or forecasting breaks."""
    trait_names = set(get_species_traits()["name"])
    count_names = set(get_species_counts()["name"])
    assert count_names.issubset(trait_names)


def test_placeholder_counts_total_is_consistent():
    counts = get_species_counts()
    assert (
        counts["total_pop"] == counts["adult_pop"] + counts["juvenile_pop"]
    ).all()


def test_run_growth_model_prints_projection_and_validation(capfd):
    run_growth_model("Zebra", start_year=2018, forecast_horizon=6)
    out = capfd.readouterr().out
    assert "Projection for Zebra from 2018:" in out
    assert "Validation against observed counts:" in out
    assert "2024:" in out