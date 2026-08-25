#!/usr/bin/env python3
"""Compare juvenile-fraction and sex-ratio assumptions by backtest error."""

from pathlib import Path

import click
import numpy as np
import pandas as pd

from growth_model import (
    ACTIVE_SPECIES,
    compare_to_counts,
    forecast_population,
    load_species_counts,
    load_species_traits,
)


DEFAULT_JUVENILE_FRACTIONS = tuple(np.arange(0.1, 0.51, 0.1))
DEFAULT_MALE_RATIOS = tuple(np.arange(0.3, 0.71, 0.1))


def tune_assumptions(
    traits_path: Path,
    counts_path: Path,
    start_year: int,
    forecast_horizon: int,
    juvenile_fractions: tuple[float, ...] = DEFAULT_JUVENILE_FRACTIONS,
    male_ratios: tuple[float, ...] = DEFAULT_MALE_RATIOS,
) -> pd.DataFrame:
    """Evaluate assumption combinations against reported census totals.

    Args:
        traits_path: Path to the normalized traits CSV.
        counts_path: Path to the normalized counts CSV.
        start_year: First year used to initialize each backtest.
        forecast_horizon: Number of annual steps to project.
        juvenile_fractions: Initial juvenile shares to test.
        male_ratios: Male population shares to test.

    Returns:
        A dataframe with one error summary per assumption combination and
        species.
    """
    traits = load_species_traits(traits_path)
    rows = []

    for juvenile_fraction in juvenile_fractions:
        counts = load_species_counts(counts_path, juvenile_fraction)
        for male_ratio in male_ratios:
            tuned_traits = traits.copy()
            tuned_traits["male_ratio"] = male_ratio
            for species_name in ACTIVE_SPECIES:
                forecast = forecast_population(
                    tuned_traits,
                    counts,
                    species_name,
                    start_year,
                    forecast_horizon,
                )
                comparison = compare_to_counts(forecast, counts, species_name)
                absolute_errors = (
                    comparison["projected_total"] - comparison["observed_total"]
                ).abs()
                rows.append(
                    {
                        "species": species_name,
                        "juvenile_fraction": juvenile_fraction,
                        "male_ratio": male_ratio,
                        "mean_absolute_error": absolute_errors.mean(),
                        "observations_used": len(comparison),
                    }
                )

    return pd.DataFrame(rows)


@click.command()
@click.option(
    "--traits-path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=Path("data/raw/species_traits.csv"),
    show_default=True,
)
@click.option(
    "--counts-path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=Path("data/raw/species_counts.csv"),
    show_default=True,
)
@click.option("--start-year", default=1965, type=int, show_default=True)
@click.option("--forecast-horizon", default=32, type=int, show_default=True)
def main(
    traits_path: Path,
    counts_path: Path,
    start_year: int,
    forecast_horizon: int,
) -> None:
    """Print the best tested assumption combination for each species."""
    results = tune_assumptions(
        traits_path, counts_path, start_year, forecast_horizon
    )
    best = results.loc[
        results.groupby("species")["mean_absolute_error"].idxmin()
    ].sort_values("species")
    click.echo(best.to_string(index=False))


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
