#!/usr/bin/env python3
"""
Simple two-stage (juvenile/adult) population model with backtest validation.

USAGE PATTERN
-------------
Historical census counts are supplied for multiple years. The model is
initialised from one chosen start year, projected forward, and then compared
against the observed counts for the years that follow. The gap between
projected and observed is the validation signal.

KNOWN SIMPLIFICATIONS (documented, to be improved later)
--------------------------------------------------------
1. Juvenile mortality falls back to adult mortality unless an explicit
   juvenile rate is supplied. In reality juvenile mortality is usually much
   higher and often drives the dynamics.
2. Carrying capacity falls back to a total-biomass estimate from the start
   year rather than a habitat/resource limit. This makes capacity partly
   circular (inferred from the population already present) and sensitive to
   which species are included in the dataset.
3. Birth rate is derived as 12 / gestation_months, i.e. continuous
   back-to-back breeding with no inter-birth interval and no seasonality.
4. Adult mortality is derived as 1 / life_span_years, i.e. a flat annual
   hazard rather than an age-structured survival curve.
5. Density dependence is applied only to the maturation term, not to births
   or mortality.
"""

from pathlib import Path
from typing import Any, Dict, Optional, Union

import click
import pandas as pd


# ---------------------------------------------------------------------------
# Data
#
# Traits and counts are held separately on purpose:
#   - traits = literature-derived constants (one row per species)
#   - counts = census observations (one row per species per year)
# Keeping them apart means each has its own provenance/source citation, and it
# stops a constant being accidentally treated as if it varied year to year.
# The committed CSVs are the reproducible default inputs for the CLI.
# ---------------------------------------------------------------------------


def get_species_traits() -> pd.DataFrame:
    """
    Static, literature-derived parameters for each species.

    Returns:
        pd.DataFrame: One row per species.
    """
    return load_species_traits("data/raw/species_traits.csv")


def get_species_counts() -> pd.DataFrame:
    """
    Observed census counts by species and year.

    Returns:
        pd.DataFrame: One row per species per year, with total_pop derived.
    """
    return load_species_counts("data/raw/species_counts.csv")


# ---------------------------------------------------------------------------
# Model components
# ---------------------------------------------------------------------------


def load_species_traits(path: Union[str, Path]) -> pd.DataFrame:
    """Load normalized species traits from a committed CSV file.

    Args:
        path: Path to the normalized traits CSV.

    Returns:
        A dataframe containing one row per species.
    """
    return pd.read_csv(path)


def load_species_counts(
    path: Union[str, Path], juvenile_fraction: float = 0.5
) -> pd.DataFrame:
    """Load total counts and seed the model's two population stages.

    Census data provides total populations only. This function applies the
    supplied juvenile fraction to each reported total; it does not fill
    missing census totals.

    Args:
        path: Path to the normalized counts CSV.
        juvenile_fraction: Assumed juvenile share of each total, between zero
            and one. Defaults to 0.5 as a provisional 50/50 adult/juvenile
            split. This is a modelling assumption, not an observation.

    Returns:
        A dataframe with ``adult_pop``, ``juvenile_pop``, and ``total_pop``.

    Raises:
        ValueError: If the juvenile fraction is outside the inclusive range
            from zero to one.
    """
    if not 0.0 <= juvenile_fraction <= 1.0:
        raise ValueError("juvenile_fraction must be between 0 and 1.")

    counts = pd.read_csv(path)
    counts["total_pop"] = pd.to_numeric(counts["total_pop"], errors="coerce")
    counts["juvenile_pop"] = counts["total_pop"] * juvenile_fraction
    counts["adult_pop"] = counts["total_pop"] * (1.0 - juvenile_fraction)
    return counts


def adult_pop_change(
    adult_pop: float, matured_juveniles: float, adult_death_rate: float
) -> float:
    """
    Net adult population change from maturation and adult mortality.

    Args:
        adult_pop (float): Current adult population.
        matured_juveniles (float): Juveniles maturing into adults this step.
        adult_death_rate (float): Annual adult mortality rate.

    Returns:
        float: Net adult population change.
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
    Net juvenile population change from births, mortality and maturation.

    Args:
        juvenile_pop (float): Current juvenile population.
        births (float): Newly born juveniles this step.
        juv_death_rate (float): Annual juvenile mortality rate.
        matured_juveniles (float): Juveniles maturing into adults this step.

    Returns:
        float: Net juvenile population change.
    """
    juvenile_deaths = juvenile_pop * juv_death_rate
    return births - juvenile_deaths - matured_juveniles


def calc_carrying_capacity(
    traits: Dict[str, Any],
    traits_df: pd.DataFrame,
    counts_df: pd.DataFrame,
    start_year: int,
) -> float:
    """
    Carrying capacity for a species, expressed in individuals.

    Uses the explicit trait value where one is supplied, otherwise falls back
    to a total-biomass estimate from the start year (see simplification 2).

    Args:
        traits (Dict[str, Any]): Trait record for the target species.
        traits_df (pd.DataFrame): All species traits (for weights).
        counts_df (pd.DataFrame): Observed counts for all species.
        start_year (int): Year to estimate capacity from.

    Returns:
        float: Carrying capacity in individuals.
    """
    explicit = traits.get("carrying_capacity")
    if explicit is not None and not pd.isna(explicit):
        return float(explicit)

    year_counts = counts_df[counts_df["year"] == start_year]
    if year_counts.empty:
        return 0.0

    merged = year_counts.merge(
        traits_df[["name", "avg_weight_kg"]], on="name", how="left"
    )
    total_kg = (merged["total_pop"] * merged["avg_weight_kg"]).sum()
    return float(total_kg / float(traits["avg_weight_kg"]))


def _normalize_record(record: Union[pd.Series, Dict[str, Any]]) -> Dict[str, Any]:
    """Convert a pandas or dictionary trait record to a dictionary.

    Args:
        record: Trait record supplied as a pandas Series or dictionary.

    Returns:
        The trait values as a dictionary.

    Raises:
        TypeError: If ``record`` is neither a pandas Series nor a dictionary.
    """
    if isinstance(record, pd.Series):
        return record.to_dict()
    if isinstance(record, dict):
        return record
    raise TypeError("record must be a pandas Series or a dict")


def population_change(
    traits: Union[pd.Series, Dict[str, Any]],
    adult_pop: float,
    juvenile_pop: float,
    carrying_capacity: Optional[float] = None,
) -> Dict[str, float]:
    """
    Advance a population by one year.

    Args:
        traits (Union[pd.Series, Dict[str, Any]]): Species trait record.
        adult_pop (float): Current adult population.
        juvenile_pop (float): Current juvenile population.
        carrying_capacity (Optional[float]): Capacity in individuals.

    Returns:
        Dict[str, float]: Updated adult, juvenile and total populations.
    """
    t = _normalize_record(traits)

    required_traits = [
        "gestation_period_months",
        "life_span_years",
        "male_ratio",
        "sex_matur_m_years",
        "sex_matur_f_years",
    ]
    missing_traits = [
        field for field in required_traits
        if t.get(field) is None or pd.isna(t.get(field))
    ]
    if missing_traits:
        raise ValueError(
            f"Missing required traits for population change: {missing_traits}"
        )

    birth_rate = 12.0 / float(t["gestation_period_months"])
    adult_death_rate = 1.0 / float(t["life_span_years"])

    juv_rate = t.get("juv_death_rate")
    juv_death_rate = (
        adult_death_rate if juv_rate is None or pd.isna(juv_rate) else float(juv_rate)
    )

    male_ratio = float(t["male_ratio"])
    female_ratio = 1.0 - male_ratio
    births = adult_pop * female_ratio * birth_rate

    sex_mat_rate = (1.0 / float(t["sex_matur_m_years"])) * male_ratio + (
        1.0 / float(t["sex_matur_f_years"])
    ) * female_ratio

    carrying_factor = 1.0
    if carrying_capacity and carrying_capacity > 0:
        carrying_factor = max(
            0.0, 1.0 - (adult_pop + juvenile_pop) / carrying_capacity
        )

    matured_juveniles = juvenile_pop * sex_mat_rate * carrying_factor

    new_adult_pop = max(
        0.0, adult_pop + adult_pop_change(adult_pop, matured_juveniles, adult_death_rate)
    )
    new_juvenile_pop = max(
        0.0,
        juvenile_pop
        + juv_pop_change(juvenile_pop, births, juv_death_rate, matured_juveniles),
    )

    return {
        "adult_pop": new_adult_pop,
        "juvenile_pop": new_juvenile_pop,
        "total_pop": new_adult_pop + new_juvenile_pop,
    }


def forecast_population(
    traits_df: pd.DataFrame,
    counts_df: pd.DataFrame,
    species_name: str,
    start_year: int,
    forecast_horizon: int,
) -> pd.DataFrame:
    """
    Project a species forward from an observed starting year.

    Args:
        traits_df (pd.DataFrame): Species traits.
        counts_df (pd.DataFrame): Observed counts.
        species_name (str): Species to forecast.
        start_year (int): Year to initialise from (must exist in counts).
        forecast_horizon (int): Number of years to project.

    Returns:
        pd.DataFrame: Year-by-year projected populations.
    """
    trait_match = traits_df[traits_df["name"] == species_name]
    if trait_match.empty:
        raise ValueError(f"No traits found for species '{species_name}'.")
    traits = trait_match.iloc[0].to_dict()

    start_counts = counts_df[
        (counts_df["name"] == species_name) & (counts_df["year"] == start_year)
    ]
    if start_counts.empty:
        raise ValueError(
            f"No observed count for '{species_name}' in {start_year} to start from."
        )

    adult_pop = float(start_counts.iloc[0]["adult_pop"])
    juvenile_pop = float(start_counts.iloc[0]["juvenile_pop"])
    carrying_capacity = calc_carrying_capacity(
        traits, traits_df, counts_df, start_year
    )

    rows = []
    for year in range(start_year, start_year + forecast_horizon + 1):
        rows.append(
            {
                "year": year,
                "name": species_name,
                "adult_pop": adult_pop,
                "juvenile_pop": juvenile_pop,
                "total_pop": adult_pop + juvenile_pop,
            }
        )
        updated = population_change(traits, adult_pop, juvenile_pop, carrying_capacity)
        adult_pop = updated["adult_pop"]
        juvenile_pop = updated["juvenile_pop"]

    return pd.DataFrame(rows)


def compare_to_counts(
    forecast_df: pd.DataFrame, counts_df: pd.DataFrame, species_name: str
) -> pd.DataFrame:
    """
    Compare projected populations against observed counts where both exist.

    Args:
        forecast_df (pd.DataFrame): Output of forecast_population.
        counts_df (pd.DataFrame): Observed counts.
        species_name (str): Species being validated.

    Returns:
        pd.DataFrame: Year-by-year projected vs observed, with error columns.
    """
    observed = counts_df[counts_df["name"] == species_name][
        ["year", "total_pop"]
    ].rename(columns={"total_pop": "observed_total"})
    observed = observed.dropna(subset=["observed_total"])

    merged = (
        forecast_df[["year", "total_pop"]]
        .rename(columns={"total_pop": "projected_total"})
        .merge(observed, on="year", how="inner")
    )

    if merged.empty:
        return merged

    merged["error"] = merged["projected_total"] - merged["observed_total"]
    merged["pct_error"] = 100.0 * merged["error"] / merged["observed_total"]
    return merged


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def run_growth_model(species_name: str, start_year: int, forecast_horizon: int) -> None:
    """
    Run the forecast and validation for one species and print the results.

    Args:
        species_name (str): Species to forecast.
        start_year (int): Year to initialise from.
        forecast_horizon (int): Number of years to project.
    """
    traits_df = load_species_traits("data/raw/species_traits.csv")
    counts_df = load_species_counts("data/raw/species_counts.csv")

    forecast_df = forecast_population(
        traits_df, counts_df, species_name, start_year, forecast_horizon
    )

    print(f"\nProjection for {species_name} from {start_year}:")
    for _, row in forecast_df.iterrows():
        print(
            f"  {int(row['year'])}: adults={row['adult_pop']:.0f} "
            f"juveniles={row['juvenile_pop']:.0f} total={row['total_pop']:.0f}"
        )

    comparison = compare_to_counts(forecast_df, counts_df, species_name)
    if comparison.empty:
        print("\nNo overlapping observed years to validate against.")
        return

    print("\nValidation against observed counts:")
    for _, row in comparison.iterrows():
        print(
            f"  {int(row['year'])}: projected={row['projected_total']:.0f} "
            f"observed={row['observed_total']:.0f} "
            f"error={row['error']:+.0f} ({row['pct_error']:+.1f}%)"
        )


@click.command()
@click.option(
    "--species_name", default="Wildebeest", help="Species to forecast."
)
@click.option("--start_year", default=1965, type=int, help="Year to initialise from.")
@click.option(
    "--forecast_horizon", default=6, type=int, help="Number of years to project."
)
def growth_model(species_name: str, start_year: int, forecast_horizon: int) -> None:
    """Run the population forecast CLI command.

    Args:
        species_name: Species to forecast.
        start_year: Observed year from which to initialize the forecast.
        forecast_horizon: Number of annual steps to project.
    """
    run_growth_model(species_name, start_year, forecast_horizon)


if __name__ == "__main__":
    growth_model()  # pylint: disable=no-value-for-parameter
