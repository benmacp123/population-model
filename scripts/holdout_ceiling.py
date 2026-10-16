#!/usr/bin/env python3
"""Diagnostic: the best holdout fit the model structure can achieve.

Fits the same three parameters directly on the 1980-1997 holdout, which a
real forecast could never do, and reports how far those parameters miss the
1965-1979 calibration period. A good holdout fit with a poor calibration fit
means the structure is not the only problem: the dynamics changed between
the two periods. None of these numbers are forecast scores.
"""

import itertools
from pathlib import Path

import click
import pandas as pd

from calibration import PARAM_GRID
from growth_model import (
    ACTIVE_SPECIES,
    calibration_peak_capacity,
    load_species_counts,
    load_species_traits,
    project_population,
)


def _totals(traits, counts, species_name, params, start_year, end_year, capacity):
    """Projected totals for one parameter set, indexed by year."""
    row = counts[(counts["name"] == species_name) & (counts["year"] == start_year)]
    record = dict(
        traits,
        juv_death_rate=params["juv_death_rate"],
        birth_multiplier=params["birth_multiplier"],
    )
    projection = project_population(
        record,
        float(row.iloc[0]["adult_pop"]),
        float(row.iloc[0]["juvenile_pop"]),
        capacity * params["capacity_multiplier"],
        start_year,
        end_year - start_year,
    )
    return projection.set_index("year")["total_pop"]


def holdout_ceiling(
    traits_path: Path,
    counts_path: Path,
    start_year: int = 1965,
    calibration_end: int = 1979,
    end_year: int = 1997,
) -> pd.DataFrame:
    """Best holdout fit per species and its error on the calibration period."""
    traits_df = load_species_traits(traits_path)
    counts = load_species_counts(counts_path)
    rows = []
    for species_name in ACTIVE_SPECIES:
        traits = traits_df[traits_df["name"] == species_name].iloc[0].to_dict()
        observed = (
            counts[counts["name"] == species_name]
            .dropna(subset=["total_pop"])
            .set_index("year")["total_pop"]
        )
        holdout = observed[observed.index > calibration_end]
        fit = observed[
            (observed.index > start_year) & (observed.index <= calibration_end)
        ]
        capacity = calibration_peak_capacity(counts, species_name, calibration_end)
        baseline_mae = (observed[calibration_end] - holdout).abs().mean()

        best_mae, best_params = None, None
        for values in itertools.product(*PARAM_GRID.values()):
            params = dict(zip(PARAM_GRID, values))
            projected = _totals(
                traits,
                counts,
                species_name,
                params,
                calibration_end,
                end_year,
                capacity,
            )
            mae = (projected[holdout.index] - holdout).abs().mean()
            if best_mae is None or mae < best_mae:
                best_mae, best_params = mae, params

        calibration = _totals(
            traits,
            counts,
            species_name,
            best_params,
            start_year,
            calibration_end,
            capacity,
        )
        rows.append(
            {
                "species": species_name,
                **best_params,
                "holdout_skill_ceiling": 1.0 - best_mae / baseline_mae,
                "calibration_nmae_pct_with_these_params": 100.0
                * (calibration[fit.index] - fit).abs().mean()
                / fit.mean(),
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
def main(traits_path: Path, counts_path: Path) -> None:
    """Print the holdout skill ceiling for the active species."""
    table = holdout_ceiling(traits_path, counts_path)
    click.echo(
        "Diagnostic only: parameters fitted on the 1980-1997 holdout, "
        "then scored on 1965-1979."
    )
    click.echo(table.to_string(index=False, float_format=lambda x: f"{x:.2f}"))


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
