#!/usr/bin/env python3
"""Calibrate on an early period and evaluate forecasts on a later holdout.

For each active species:

1. Fit juvenile mortality, a birth multiplier and a capacity multiplier on
   ``start-year`` to ``calibration-end`` with a bounded grid search.
2. Restart the forecast from the last reported total in the calibration
   period and project to ``end-year``.
3. Score the holdout years against a naive persistence baseline that carries
   that same anchor total forward unchanged.

Results are written to ``results/`` as a compact summary table and a
year-by-year projections table used by ``scripts/plot_backtest.py``.
"""

from pathlib import Path

import click
import numpy as np
import pandas as pd

from calibration import PARAM_GRID, calibrate_species, holdout_forecast, near_best
from growth_model import (
    ACTIVE_SPECIES,
    backtest_metrics,
    load_species_counts,
    load_species_traits,
)


def _on_bound(value: float, name: str) -> bool:
    """Whether a fitted value sits on the edge of its search range."""
    return value in (min(PARAM_GRID[name]), max(PARAM_GRID[name]))


def evaluate_model(
    traits_path: Path,
    counts_path: Path,
    start_year: int,
    calibration_end: int,
    end_year: int,
    juvenile_fraction: float = 0.5,
    male_ratio: float = 0.5,
    tolerance: float = 0.1,
):
    """Calibrate, forecast and score each active species.

    Args:
        traits_path: Path to the normalized traits CSV.
        counts_path: Path to the normalized counts CSV.
        start_year: First calibration year (projection start).
        calibration_end: Last calibration year.
        end_year: Last projected year.
        juvenile_fraction: Assumed juvenile share of observed totals.
        male_ratio: Assumed male share of the population.
        tolerance: Relative calibration-error tolerance for the scenario set.

    Returns:
        A ``(summary, projections)`` pair of dataframes.
    """
    traits_df = load_species_traits(traits_path)
    traits_df["male_ratio"] = male_ratio
    counts = load_species_counts(counts_path, juvenile_fraction)
    summaries, projections = [], []

    for species_name in ACTIVE_SPECIES:
        traits = traits_df[traits_df["name"] == species_name].iloc[0].to_dict()
        fits = calibrate_species(
            traits, counts, species_name, start_year, calibration_end
        )
        scenario_fits = near_best(fits, tolerance)
        best = fits.iloc[0]

        observed = (
            counts[counts["name"] == species_name]
            .dropna(subset=["total_pop"])
            .set_index("year")["total_pop"]
        )
        anchor_year = int(observed[observed.index <= calibration_end].index.max())
        anchor_total = float(observed[anchor_year])

        fitted = holdout_forecast(
            traits,
            counts,
            species_name,
            fits.head(1),
            start_year,
            calibration_end,
            calibration_end,
        )
        forecast = holdout_forecast(
            traits,
            counts,
            species_name,
            scenario_fits,
            anchor_year,
            end_year,
            calibration_end,
        )

        holdout_obs = observed[observed.index > anchor_year]
        holdout = pd.DataFrame(
            {
                "observed_total": holdout_obs,
                "projected_total": forecast.loc[holdout_obs.index, "central"],
            }
        )
        within = (holdout_obs >= forecast.loc[holdout_obs.index, "low"]) & (
            holdout_obs <= forecast.loc[holdout_obs.index, "high"]
        )
        fit_obs = observed[
            (observed.index > start_year) & (observed.index <= calibration_end)
        ]
        params = {
            "juv_death_rate": best["juv_death_rate"],
            "birth_multiplier": best["birth_multiplier"],
            "capacity_multiplier": best["capacity_multiplier"],
        }
        summaries.append(
            {
                "species": species_name,
                **params,
                "params_on_bound": ", ".join(
                    name for name, value in params.items() if _on_bound(value, name)
                )
                or "none",
                "scenario_sets": len(scenario_fits),
                "calibration_nmae_pct": 100.0
                * best["calibration_mae"]
                / fit_obs.mean(),
                **backtest_metrics(holdout, anchor_total),
                "band_coverage": float(within.mean()),
            }
        )

        for period, frame in (("calibration", fitted), ("holdout", forecast)):
            projections.append(
                frame.assign(species=species_name, period=period)
                .rename_axis("year")
                .reset_index()
            )
        projections.append(
            pd.DataFrame(
                {
                    "year": observed.index,
                    "species": species_name,
                    "period": "observed",
                    "central": observed.to_numpy(),
                    "low": np.nan,
                    "high": np.nan,
                }
            )
        )

    return pd.DataFrame(summaries), pd.concat(projections, ignore_index=True)


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
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("results"),
    show_default=True,
)
def main(
    traits_path: Path,
    counts_path: Path,
    start_year: int,
    calibration_end: int,
    end_year: int,
    juvenile_fraction: float,
    male_ratio: float,
    output_dir: Path,
) -> None:
    """Print and save holdout accuracy for the active species."""
    if not start_year < calibration_end < end_year:
        raise click.BadParameter(
            "Require start-year < calibration-end < end-year.",
            param_hint="--calibration-end",
        )
    summary, projections = evaluate_model(
        traits_path,
        counts_path,
        start_year,
        calibration_end,
        end_year,
        juvenile_fraction,
        male_ratio,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    summary.round(3).to_csv(output_dir / "backtest_summary.csv", index=False)
    projections.round(1).to_csv(output_dir / "backtest_projections.csv", index=False)

    click.echo(
        f"Calibration {start_year}-{calibration_end}; holdout to {end_year}; "
        f"juvenile_fraction={juvenile_fraction} (provisional); "
        f"male_ratio={male_ratio} (provisional)"
    )
    click.echo(summary.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    click.echo(
        f"\nWrote {output_dir}/backtest_summary.csv and backtest_projections.csv"
    )


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
