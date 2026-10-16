"""Bounded grid calibration and holdout forecasting for the growth model.

Calibration fits three parameters per species on an early period only:

- ``juv_death_rate``: annual juvenile mortality
- ``birth_multiplier``: share of the one-young-per-interbirth-interval
  maximum actually realised
- ``capacity_multiplier``: carrying capacity as a multiple of the highest
  reported total in the calibration period

The bounds are biologically defensible ranges for a sensitivity search, not
sourced values. Parameter sets whose calibration error is within a tolerance
of the best fit form the low/central/high scenario envelope.
"""

import itertools
from typing import Dict, Sequence

import numpy as np
import pandas as pd

from growth_model import calibration_peak_capacity, project_population

PARAM_GRID: Dict[str, Sequence[float]] = {
    "juv_death_rate": tuple(np.round(np.arange(0.10, 0.601, 0.05), 2)),
    "birth_multiplier": tuple(np.round(np.arange(0.50, 1.001, 0.10), 2)),
    "capacity_multiplier": tuple(np.round(np.arange(0.50, 1.501, 0.10), 2)),
}


def _observed_totals(counts_df: pd.DataFrame, species_name: str) -> pd.Series:
    """Reported totals for one species, indexed by year."""
    species = counts_df[counts_df["name"] == species_name].dropna(subset=["total_pop"])
    return species.set_index("year")["total_pop"].sort_index()


def _start_stages(counts_df: pd.DataFrame, species_name: str, year: int):
    """Seeded adult and juvenile counts for a reported year."""
    row = counts_df[(counts_df["name"] == species_name) & (counts_df["year"] == year)]
    if row.empty or pd.isna(row.iloc[0]["total_pop"]):
        raise ValueError(f"No reported total for '{species_name}' in {year}.")
    return float(row.iloc[0]["adult_pop"]), float(row.iloc[0]["juvenile_pop"])


def _project_totals(
    traits: Dict,
    params: Dict[str, float],
    stages,
    capacity_base: float,
    start_year: int,
    end_year: int,
) -> pd.Series:
    """Projected totals from ``start_year`` to ``end_year`` for one parameter set."""
    record = dict(
        traits,
        juv_death_rate=params["juv_death_rate"],
        birth_multiplier=params["birth_multiplier"],
    )
    projection = project_population(
        record,
        stages[0],
        stages[1],
        capacity_base * params["capacity_multiplier"],
        start_year,
        end_year - start_year,
    )
    return projection.set_index("year")["total_pop"]


def calibrate_species(
    traits: Dict,
    counts_df: pd.DataFrame,
    species_name: str,
    start_year: int,
    calibration_end: int,
    grid: Dict[str, Sequence[float]] = None,
) -> pd.DataFrame:
    """Score every grid parameter set on the calibration period.

    The projection starts from the observed ``start_year`` total and is
    scored on reported years after it, up to ``calibration_end``.

    Args:
        traits: Trait record for the species.
        counts_df: Observed counts with seeded stages.
        species_name: Species to calibrate.
        start_year: Year the calibration projection starts from.
        calibration_end: Last year used for fitting.
        grid: Parameter grid; defaults to ``PARAM_GRID``.

    Returns:
        One row per parameter set with its calibration mean absolute error,
        sorted best first.
    """
    grid = PARAM_GRID if grid is None else grid
    observed = _observed_totals(counts_df, species_name)
    fit_years = observed[
        (observed.index > start_year) & (observed.index <= calibration_end)
    ]
    if fit_years.empty:
        raise ValueError(f"No calibration years for '{species_name}'.")

    stages = _start_stages(counts_df, species_name, start_year)
    capacity_base = calibration_peak_capacity(counts_df, species_name, calibration_end)
    names = list(grid)
    rows = []
    for values in itertools.product(*(grid[name] for name in names)):
        params = dict(zip(names, values))
        projected = _project_totals(
            traits, params, stages, capacity_base, start_year, calibration_end
        )
        mae = float((projected[fit_years.index] - fit_years).abs().mean())
        rows.append({**params, "calibration_mae": mae})

    return pd.DataFrame(rows).sort_values("calibration_mae").reset_index(drop=True)


def near_best(fits: pd.DataFrame, tolerance: float = 0.1) -> pd.DataFrame:
    """Parameter sets within ``tolerance`` of the best calibration error."""
    best = fits["calibration_mae"].iloc[0]
    return fits[fits["calibration_mae"] <= best * (1.0 + tolerance)]


def holdout_forecast(
    traits: Dict,
    counts_df: pd.DataFrame,
    species_name: str,
    fits: pd.DataFrame,
    anchor_year: int,
    end_year: int,
    calibration_end: int,
) -> pd.DataFrame:
    """Forecast the holdout period from the last calibration observation.

    The central forecast uses the best-fitting parameters. Low and high are
    the year-by-year minimum and maximum across the near-best set.

    Args:
        traits: Trait record for the species.
        counts_df: Observed counts with seeded stages.
        species_name: Species to forecast.
        fits: Output of ``near_best`` (best fit first).
        anchor_year: Reported year the forecast restarts from.
        end_year: Last projected year.
        calibration_end: Last calibration year (sets the capacity base).

    Returns:
        Year-indexed ``central``, ``low`` and ``high`` projected totals.
    """
    stages = _start_stages(counts_df, species_name, anchor_year)
    capacity_base = calibration_peak_capacity(counts_df, species_name, calibration_end)
    param_names = [name for name in fits.columns if name != "calibration_mae"]
    runs = pd.DataFrame(
        {
            i: _project_totals(
                traits,
                row[param_names].to_dict(),
                stages,
                capacity_base,
                anchor_year,
                end_year,
            )
            for i, row in fits.iterrows()
        }
    )
    return pd.DataFrame(
        {
            "central": runs.iloc[:, 0],
            "low": runs.min(axis=1),
            "high": runs.max(axis=1),
        }
    )
