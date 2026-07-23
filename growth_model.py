#!/usr/bin/env python3

import click
from typing import Any, Dict, Tuple, Union

import pandas as pd


def adult_pop_change(
    adult_pop: float, matured_juveniles: float, adult_death_rate: float
) -> float:
    """
    Calculate net adult population change from maturation and adult mortality.

    Args:
        adult_pop (float): Current adult population.
        matured_juveniles (float): Number of juveniles maturing into adults.
        adult_death_rate (float): Annual adult mortality rate.

    Returns:
        float: The net adult population change.
    """
    adult_deaths = adult_pop * adult_death_rate
    return matured_juveniles - adult_deaths


def juv_pop_change(
    juvenile_pop: float,
    births: float,
    juv_death_rate: float,
    matured_juveniles: float,
) -> float:
    """
    Calculate net juvenile population change from births, juvenile mortality, and maturation.

    Args:
        juvenile_pop (float): Current juvenile population.
        births (float): Number of newly born juveniles.
        juv_death_rate (float): Annual juvenile mortality rate.
        matured_juveniles (float): Number of juveniles maturing into adults.

    Returns:
        float: The net juvenile population change.
    """
    juvenile_deaths = juvenile_pop * juv_death_rate
    return births - juvenile_deaths - matured_juveniles


def get_species_df() -> pd.DataFrame:
    """
    Placeholder function to retrieve a list of species records.

    Returns:
        pd.DataFrame: A DataFrame containing species data with columns for name, population, and average weight.
    """
    return pd.DataFrame(
        [
            {
                "name": "Zebra",
                "year": 2024,
                "avg_weight_kg": 10.0,
                "adult_pop": 1000.0,
                "juvenile_pop": 500.0,
                "total_pop": 1500.0,
                "gestation_period_months": 6,
                "life_span_years": 5,
                "sex_matur_m_years": 2,
                "sex_matur_f_years": 2,
                "sex_ratio": 0.5,
            },
            {
                "name": "Waterbuck",
                "year": 2024,
                "avg_weight_kg": 20.0,
                "adult_pop": 800.0,
                "juvenile_pop": 400.0,
                "total_pop": 1200.0,
                "gestation_period_months": 8,
                "life_span_years": 6,
                "sex_matur_m_years": 3,
                "sex_matur_f_years": 3,
                "sex_ratio": 0.6,
            },
        ]
    )


def calc_carrying_capacity(species_df: pd.DataFrame, species_name: str) -> float:
    """
    Calculate the carrying capacity for a given species based on total biomass and average weight.

    Args:
        species_df (pd.DataFrame): Species data.
        species_name (str): The name of the species for which to calculate carrying capacity.

    Returns:
        float: The carrying capacity for the specified species.
    """
    if species_df.empty:
        return 1000.0

    total_kg = (species_df["total_pop"] * species_df["avg_weight_kg"]).sum()
    species_match = species_df[species_df["name"] == species_name]
    if species_match.empty:
        species_avg_weight_kg = 10.0
    else:
        species_avg_weight_kg = float(species_match.iloc[0]["avg_weight_kg"])

    return total_kg / species_avg_weight_kg


def _normalize_species_details(
    species_details: Union[pd.Series, Dict[str, Any]]
) -> Dict[str, Any]:
    if isinstance(species_details, pd.Series):
        return species_details.to_dict()
    if isinstance(species_details, dict):
        return species_details
    raise TypeError(
        "species_details must be a pandas Series or a dict containing species fields"
    )


def population_change(
    species_details: Union[pd.Series, Dict[str, Any]], carrying_capacity: float = None
) -> pd.Series:
    """
    Calculate the annual population update for a species record.

    Args:
        species_details (Union[pd.Series, Dict[str, Any]]): Species data record.
        carrying_capacity (float, optional): Carrying capacity for the species.

    Returns:
        pd.Series: Updated adult, juvenile, and total population values.
    """
    details = _normalize_species_details(species_details)
    adult_pop = float(details["adult_pop"])
    juvenile_pop = float(details["juvenile_pop"])

    birth_rate = 12.0 / float(details["gestation_period_months"])
    adult_death_rate = 1.0 / float(details["life_span_years"])
    juv_death_rate = adult_death_rate

    female_ratio = 1.0 - float(details["sex_ratio"])
    births = adult_pop * female_ratio * birth_rate

    sex_mat_rate_m = 1.0 / float(details["sex_matur_m_years"])
    sex_mat_rate_f = 1.0 / float(details["sex_matur_f_years"])
    sex_mat_rate_years = (
        sex_mat_rate_m * float(details["sex_ratio"])
        + sex_mat_rate_f * female_ratio
    )

    carrying_factor = 1.0
    if carrying_capacity and carrying_capacity > 0:
        carrying_factor = max(
            0.0,
            1.0 - (adult_pop + juvenile_pop) / carrying_capacity,
        )

    matured_juveniles = juvenile_pop * sex_mat_rate_years * carrying_factor
    adult_change = adult_pop_change(
        adult_pop=adult_pop,
        matured_juveniles=matured_juveniles,
        adult_death_rate=adult_death_rate,
    )
    juv_change = juv_pop_change(
        juvenile_pop=juvenile_pop,
        births=births,
        juv_death_rate=juv_death_rate,
        matured_juveniles=matured_juveniles,
    )

    new_adult_pop = adult_pop + adult_change
    new_juvenile_pop = juvenile_pop + juv_change
    new_total_pop = new_adult_pop + new_juvenile_pop

    return pd.Series(
        {
            "adult_pop": new_adult_pop,
            "juvenile_pop": new_juvenile_pop,
            "total_pop": new_total_pop,
        }
    )


def forecast_population(
    species_df: pd.DataFrame,
    species_name: str,
    start_year: int = 2024,
    forecast_horizon: int = 10,
) -> pd.DataFrame:
    """
    Forecast population counts for a species as a DataFrame.

    Args:
        species_df (pd.DataFrame): Species data.
        species_name (str): Name of the species to forecast.
        start_year (int): Forecast start year.
        forecast_horizon (int): Number of years to forecast.

    Returns:
        pd.DataFrame: Year-by-year forecast of adult, juvenile, and total population.
    """
    species_match = species_df[species_df["name"] == species_name]
    if species_match.empty:
        raise ValueError(f"Species '{species_name}' not found in the data.")

    species_details = species_match.iloc[0].copy()
    carrying_capacity = calc_carrying_capacity(species_df, species_name)

    rows = []
    for year in range(start_year, start_year + forecast_horizon + 1):
        rows.append(
            {
                "year": year,
                "name": species_name,
                "adult_pop": float(species_details["adult_pop"]),
                "juvenile_pop": float(species_details["juvenile_pop"]),
                "total_pop": float(species_details["total_pop"]),
            }
        )
        updated = population_change(species_details, carrying_capacity)
        species_details["adult_pop"] = updated["adult_pop"]
        species_details["juvenile_pop"] = updated["juvenile_pop"]
        species_details["total_pop"] = updated["total_pop"]

    return pd.DataFrame(rows)


def run_growth_model(species_name: str, start_year: int, forecast_horizon: int) -> None:
    """
    Run the growth model for a specified species and time period.

    Args:
        species_name (str): The name of the species to forecast.
        start_year (int): The starting year for the forecast.
        forecast_horizon (int): The number of years to forecast.
    """
    species_df = get_species_df()
    if species_df.empty:
        print("No species data available.")
        return

    forecast_df = forecast_population(species_df, species_name, start_year, forecast_horizon)
    carrying_capacity = calc_carrying_capacity(species_df, species_name)
    print(f"Carrying capacity for {species_name}: {carrying_capacity}")
    print(f"Population forecast for {species_name}:")
    for _, row in forecast_df.iterrows():
        print(
            f"  {int(row['year'])}: Adults = {row['adult_pop']}, Juveniles = {row['juvenile_pop']}, Total = {row['total_pop']}"
        )


@click.command()
@click.option(
    "--species_name",
    default="Zebra",
    help="Name of the species to forecast.",
)
@click.option(
    "--start_year",
    default=2024,
    type=int,
    help="The starting year for the forecast.",
)
@click.option(
    "--forecast_horizon",
    default=10,
    type=int,
    help="Number of years to forecast.",
)
def growth_model(
    species_name: str, start_year: int, forecast_horizon: int
) -> None:  # pylint: disable=no-value-for-parameter
    run_growth_model(species_name, start_year, forecast_horizon)


if __name__ == "__main__":
    growth_model.main()  # pylint: disable=no-value-for-parameter
