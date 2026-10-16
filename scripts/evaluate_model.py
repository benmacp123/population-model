#!/usr/bin/env python3
"""Evaluate population forecasts against reported census totals.

Each forecast starts from the observed start-year total. Years up to
``calibration_end`` are only used to set the diagnostic carrying capacity;
accuracy is scored on the later holdout years, so neither the start year nor
the capacity's source years count as independent predictions.
"""

from pathlib import Path
from typing import Sequence

import click
import pandas as pd

from growth_model import (
    ACTIVE_SPECIES,
    backtest_metrics,
    calc_carrying_capacity,
    calibration_peak_capacity,
    compare_to_counts,
    forecast_population,
    load_species_counts,
    load_species_traits,
)

# Annual juvenile mortality scenarios. These span a plausible range for a
# sensitivity test; they are not literature-sourced values.
JUV_DEATH_SCENARIOS = {"low": 0.2, "central": 0.35, "high": 0.5}
CAPACITY_MODES = ("calibration-peak", "community-biomass")


def evaluate_model(
    traits_path: Path,
    counts_path: Path,
    start_year: int,
    calibration_end: int,
    end_year: int,
    capacity_mode: str = "calibration-peak",
    juvenile_fraction: float = 0.5,
    male_ratio: float = 0.5,
    scenarios: Sequence[str] = tuple(JUV_DEATH_SCENARIOS),
) -> pd.DataFrame:
    """Calculate holdout accuracy for each species and mortality scenario.

    Args:
        traits_path: Path to the normalized traits CSV.
        counts_path: Path to the normalized counts CSV.
        start_year: Year from which each forecast is initialized.
        calibration_end: Last year of the calibration period.
        end_year: Last year projected and scored.
        capacity_mode: ``calibration-peak`` caps each species at its highest
            calibration-period total; ``community-biomass`` uses the model's
            start-year biomass proxy.
        juvenile_fraction: Assumed juvenile share of observed totals.
        male_ratio: Assumed male share of the population.
        scenarios: Names from ``JUV_DEATH_SCENARIOS`` to run.

    Returns:
        A dataframe with one row per species and scenario.
    """
    traits = load_species_traits(traits_path)
    traits["male_ratio"] = male_ratio
    counts = load_species_counts(counts_path, juvenile_fraction)
    rows = []

    for scenario in scenarios:
        scenario_traits = traits.copy()
        scenario_traits["juv_death_rate"] = JUV_DEATH_SCENARIOS[scenario]

        for species_name in ACTIVE_SPECIES:
            species_mask = scenario_traits["name"] == species_name
            if capacity_mode == "calibration-peak":
                capacity = calibration_peak_capacity(
                    counts, species_name, calibration_end
                )
            else:
                capacity = calc_carrying_capacity(
                    scenario_traits[species_mask].iloc[0].to_dict(),
                    scenario_traits,
                    counts,
                    start_year,
                )
            scenario_traits.loc[species_mask, "carrying_capacity"] = capacity

            forecast = forecast_population(
                scenario_traits,
                counts,
                species_name,
                start_year,
                end_year - start_year,
            )
            comparison = compare_to_counts(forecast, counts, species_name)
            holdout = comparison[comparison["year"] > calibration_end]
            if holdout.empty:
                continue

            start_total = comparison.loc[
                comparison["year"] == start_year, "observed_total"
            ].iloc[0]
            rows.append(
                {
                    "species": species_name,
                    "juv_mortality": scenario,
                    "capacity": round(capacity),
                    **backtest_metrics(holdout, start_total),
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
@click.option("--calibration-end", default=1979, type=int, show_default=True)
@click.option("--end-year", default=1997, type=int, show_default=True)
@click.option(
    "--capacity-mode",
    type=click.Choice(CAPACITY_MODES),
    default="calibration-peak",
    show_default=True,
)
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
    calibration_end: int,
    end_year: int,
    capacity_mode: str,
    juvenile_fraction: float,
    male_ratio: float,
) -> None:
    """Print holdout accuracy metrics for the active species."""
    if not start_year <= calibration_end < end_year:
        raise click.BadParameter(
            "Require start-year <= calibration-end < end-year.",
            param_hint="--calibration-end",
        )
    results = evaluate_model(
        traits_path,
        counts_path,
        start_year,
        calibration_end,
        end_year,
        capacity_mode,
        juvenile_fraction,
        male_ratio,
    )
    click.echo(
        f"Start {start_year}; calibration to {calibration_end}; holdout "
        f"{calibration_end + 1}-{end_year}; capacity={capacity_mode}; "
        f"juvenile_fraction={juvenile_fraction} (provisional); "
        f"male_ratio={male_ratio} (provisional)"
    )
    click.echo(
        "Juvenile mortality scenarios (sensitivity range, not sourced): "
        + ", ".join(f"{k}={v}" for k, v in JUV_DEATH_SCENARIOS.items())
    )
    click.echo(results.to_string(index=False, float_format=lambda x: f"{x:.2f}"))


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
