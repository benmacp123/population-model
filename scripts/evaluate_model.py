#!/usr/bin/env python3
"""Evaluate population forecasts against reported census totals."""

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


def evaluate_model(
    traits_path: Path,
    counts_path: Path,
    start_year: int,
    forecast_horizon: int,
    juvenile_fraction: float = 0.5,
    male_ratio: float = 0.5,
) -> pd.DataFrame:
    """Calculate backtest accuracy metrics for each forecastable species.

    Args:
        traits_path: Path to the normalized traits CSV.
        counts_path: Path to the normalized counts CSV.
        start_year: Year from which each forecast is initialized.
        forecast_horizon: Number of annual steps to project.
        juvenile_fraction: Assumed juvenile share of observed totals.
        male_ratio: Assumed male share of the population.

    Returns:
        A dataframe containing one accuracy summary per species.
    """
    traits = load_species_traits(traits_path)
    traits["male_ratio"] = male_ratio
    counts = load_species_counts(counts_path, juvenile_fraction)
    summaries = []

    for species_name in ACTIVE_SPECIES:
        try:
            forecast = forecast_population(
                traits, counts, species_name, start_year, forecast_horizon
            )
        except ValueError:
            continue

        comparison = compare_to_counts(forecast, counts, species_name)
        if comparison.empty:
            continue
        errors = comparison["projected_total"] - comparison["observed_total"]
        observed_mean = comparison["observed_total"].mean()
        summaries.append(
            {
                "species": species_name,
                "observations_used": len(comparison),
                "mae": errors.abs().mean(),
                "rmse": np.sqrt((errors**2).mean()),
                "mape_pct": (
                    100.0 * errors.abs() / comparison["observed_total"]
                ).mean(),
                "normalized_mae_pct": 100.0 * errors.abs().mean() / observed_mean,
                "mean_bias": errors.mean(),
            }
        )

    return pd.DataFrame(summaries).sort_values("species").reset_index(drop=True)


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
@click.option(
    "--juvenile-fraction",
    default=0.5,
    type=click.FloatRange(0.0, 1.0),
    show_default=True,
)
@click.option(
    "--male-ratio",
    default=0.5,
    type=click.FloatRange(0.0, 1.0),
    show_default=True,
)
def main(
    traits_path: Path,
    counts_path: Path,
    start_year: int,
    forecast_horizon: int,
    juvenile_fraction: float,
    male_ratio: float,
) -> None:
    """Print backtest accuracy metrics for the active species."""
    results = evaluate_model(
        traits_path,
        counts_path,
        start_year,
        forecast_horizon,
        juvenile_fraction,
        male_ratio,
    )
    click.echo(results.to_string(index=False))


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
