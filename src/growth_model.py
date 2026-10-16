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
3. Birth rate is one young per interbirth interval (365.25 / interbirth
   days) with no seasonality. Where the interval is missing it falls back to
   12 / gestation_months, i.e. continuous back-to-back breeding.
4. Adult mortality is derived as 1 / life_span_years, i.e. a flat annual
   hazard rather than an age-structured survival curve.
5. Density dependence is applied only to the maturation term, not to births
   or mortality.
"""

from pathlib import Path
from typing import Any, Dict, Optional, Union

import click
import numpy as np
import pandas as pd

ACTIVE_SPECIES = ("Wildebeest", "Waterbuck", "Kudu")
DAYS_PER_YEAR = 365.25


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


def annual_birth_rate(traits: Dict[str, Any]) -> float:
    """Young per adult female per year (see simplification 3).

    Uses the interbirth interval where available, otherwise assumes
    back-to-back breeding from gestation length.

    Args:
        traits: Species trait record.

    Returns:
        Annual births per adult female.
    """
    interval = traits.get("interbirth_interval_days")
    if interval is not None and not pd.isna(interval) and float(interval) > 0:
        return DAYS_PER_YEAR / float(interval)
    return 12.0 / float(traits["gestation_period_months"])


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
        field
        for field in required_traits
        if t.get(field) is None or pd.isna(t.get(field))
    ]
    if missing_traits:
        raise ValueError(
            f"Missing required traits for population change: {missing_traits}"
        )

    birth_rate = annual_birth_rate(t)
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
        carrying_factor = max(0.0, 1.0 - (adult_pop + juvenile_pop) / carrying_capacity)

    matured_juveniles = juvenile_pop * sex_mat_rate * carrying_factor

    new_adult_pop = max(
        0.0,
        adult_pop + adult_pop_change(adult_pop, matured_juveniles, adult_death_rate),
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
    carrying_capacity = calc_carrying_capacity(traits, traits_df, counts_df, start_year)

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


def calibration_peak_capacity(
    counts_df: pd.DataFrame, species_name: str, calibration_end: int
) -> float:
    """Diagnostic capacity: the highest reported total up to a cutoff year.

    This is circular by construction (inferred from the species' own counts)
    and is a stand-in for a habitat-based capacity, not an estimate of one.
    Only years up to ``calibration_end`` are used, so later years remain an
    independent holdout.

    Args:
        counts_df: Observed counts.
        species_name: Species to estimate capacity for.
        calibration_end: Last year included in the calibration period.

    Returns:
        Capacity in individuals.

    Raises:
        ValueError: If the species has no reported total in the period.
    """
    observed = counts_df[
        (counts_df["name"] == species_name) & (counts_df["year"] <= calibration_end)
    ]["total_pop"].dropna()
    if observed.empty:
        raise ValueError(
            f"No reported totals for '{species_name}' up to {calibration_end}."
        )
    return float(observed.max())


def backtest_metrics(comparison: pd.DataFrame, baseline_total: float) -> Dict:
    """Score a projection on level and trend against a naive baseline.

    The naive baseline holds the population constant at ``baseline_total``
    (normally the start-year count). Skill above zero means the model's mean
    absolute error is lower than the baseline's.

    Args:
        comparison: Output of ``compare_to_counts`` restricted to the years
            being scored (exclude the initialization year).
        baseline_total: Constant population used by the naive baseline.

    Returns:
        A dictionary of accuracy metrics.
    """
    observed = comparison["observed_total"].to_numpy(dtype=float)
    projected = comparison["projected_total"].to_numpy(dtype=float)
    errors = projected - observed
    mae = float(np.abs(errors).mean())
    baseline_mae = float(np.abs(baseline_total - observed).mean())
    observed_mean = float(observed.mean())

    step_agreement = np.nan
    if len(observed) > 1:
        step_agreement = float(
            (np.sign(np.diff(observed)) == np.sign(np.diff(projected))).mean()
        )

    return {
        "years_scored": len(observed),
        "normalized_mae_pct": 100.0 * mae / observed_mean,
        "bias_pct": 100.0 * float(errors.mean()) / observed_mean,
        "baseline_normalized_mae_pct": 100.0 * baseline_mae / observed_mean,
        "skill_vs_baseline": (1.0 - mae / baseline_mae if baseline_mae > 0 else np.nan),
        "step_direction_agreement": step_agreement,
        "overall_trend_match": bool(
            np.sign(observed[-1] - observed[0]) == np.sign(projected[-1] - projected[0])
        ),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def run_growth_model(
    species_name: str,
    start_year: int,
    forecast_horizon: int,
    juvenile_fraction: float = 0.5,
    male_ratio: float = 0.5,
    juv_death_rate: Optional[float] = None,
) -> None:
    """
    Run the forecast and validation for one species and print the results.

    Args:
        species_name (str): Species to forecast.
        start_year (int): Year to initialise from.
        forecast_horizon (int): Number of years to project.
        juvenile_fraction: Assumed juvenile share of observed totals.
        male_ratio: Assumed male share of the population.
        juv_death_rate: Annual juvenile mortality override. ``None`` keeps
            the traits table value, which falls back to adult mortality.
    """
    traits_df = load_species_traits("data/raw/species_traits.csv")
    traits_df["male_ratio"] = male_ratio
    if juv_death_rate is not None:
        traits_df["juv_death_rate"] = juv_death_rate
    counts_df = load_species_counts("data/raw/species_counts.csv", juvenile_fraction)

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
@click.option("--species_name", default="Wildebeest", help="Species to forecast.")
@click.option("--start_year", default=1965, type=int, help="Year to initialise from.")
@click.option(
    "--forecast_horizon", default=6, type=int, help="Number of years to project."
)
@click.option(
    "--juvenile_fraction",
    default=0.5,
    type=click.FloatRange(0.0, 1.0),
    show_default=True,
    help="Assumed juvenile share of each observed total.",
)
@click.option(
    "--male_ratio",
    default=0.5,
    type=click.FloatRange(0.0, 1.0),
    show_default=True,
    help="Assumed male share of the population.",
)
@click.option(
    "--juv_death_rate",
    default=None,
    type=click.FloatRange(0.0, 1.0),
    help="Annual juvenile mortality. Defaults to the adult rate.",
)
def growth_model(
    species_name: str,
    start_year: int,
    forecast_horizon: int,
    juvenile_fraction: float,
    male_ratio: float,
    juv_death_rate: Optional[float],
) -> None:
    """Run the population forecast CLI command.

    Args:
        species_name: Species to forecast.
        start_year: Observed year from which to initialize the forecast.
        forecast_horizon: Number of annual steps to project.
        juvenile_fraction: Assumed juvenile share of observed totals.
        male_ratio: Assumed male share of the population.
        juv_death_rate: Annual juvenile mortality override.
    """
    run_growth_model(
        species_name,
        start_year,
        forecast_horizon,
        juvenile_fraction,
        male_ratio,
        juv_death_rate,
    )


if __name__ == "__main__":
    growth_model()  # pylint: disable=no-value-for-parameter
