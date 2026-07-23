import pandas as pd

from growth_model import (
    adult_pop_change,
    calc_carrying_capacity,
    get_species_df,
    juv_pop_change,
    population_change,
    run_growth_model,
)


def test_adult_pop_change():
    assert adult_pop_change(adult_pop=100, matured_juveniles=10, adult_death_rate=0.05) == 10 - 5


def test_juv_pop_change():
    assert juv_pop_change(juvenile_pop=50, births=20, juv_death_rate=0.1, matured_juveniles=5) == 20 - 5 - 5


def test_calc_carrying_capacity_default():
    species_df = get_species_df()
    capacity = calc_carrying_capacity("Zebra")
    assert capacity > 0
    assert isinstance(capacity, float)


def test_population_change_uses_population_scaling():
    species_details = {
        "adult_pop": 100.0,
        "juvenile_pop": 50.0,
        "gestation_period_months": 6,
        "life_span_years": 10,
        "sex_matur_m_years": 2,
        "sex_matur_f_years": 2,
        "sex_ratio": 0.5,
    }
    new_adult_pop, new_juvenile_pop = population_change(species_details, carrying_capacity=500.0)
    assert new_adult_pop != species_details["adult_pop"]
    assert new_juvenile_pop != species_details["juvenile_pop"]
    assert new_adult_pop > 0
    assert new_juvenile_pop > 0


def test_run_growth_model_prints_forecast(capfd):
    run_growth_model("Zebra", start_year=2024, forecast_horizon=2)
    captured = capfd.readouterr()
    assert "Population forecast for Zebra:" in captured.out
    assert "2024:" in captured.out or "2025:" in captured.out
