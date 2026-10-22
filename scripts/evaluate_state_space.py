#!/usr/bin/env python3
"""Rolling-origin one-step evaluation of the version 2 state-space model.

At each origin year the models are refitted on regional counts up to and
including that year, then forecast the following year for every region
observed in it. Four forecasters are compared on the log scale:

- ``rain``: pooled Gompertz state-space model with lagged dry-season rain
- ``no_rain``: the same model with the rainfall coefficient fixed at zero
- ``random_walk``: last observed log count, sd from past one-year changes
- ``calibration_mean``: each region's mean log count up to the origin

The design was fixed before scoring (see the README). Development species are used to
build the pipeline; the confirmation species are run once, with
``--confirmation``, after the pipeline is frozen.
"""

from pathlib import Path

import click
import numpy as np
import pandas as pd
from scipy.stats import norm

from state_space import fit, forecast_next, prepare_series

DEVELOPMENT_SPECIES = ("Wildebeest", "Waterbuck", "Kudu")
CONFIRMATION_SPECIES = (
    "Eland",
    "Giraffe",
    "Impala",
    "Roan",
    "Sable",
    "Tsessebe",
    "Zebra",
    "Warthog",
)
FIRST_ORIGIN, LAST_ORIGIN = 1979, 1994
RAIN_STANDARDISE_YEARS = (1962, 1979)


def lagged_rain(rainfall: pd.DataFrame, years: np.ndarray, column: str, lag: int):
    """Standardised rainfall acting on the transition into each year.

    Args:
        rainfall: Table from ``data/raw/kruger_rainfall_cru.csv``.
        years: Model years.
        column: Rainfall column to use.
        lag: Seasonal-year lag; 1 means the season ending at the previous
            census, which is known at the forecast origin.

    Returns:
        Array aligned with ``years``.

    Raises:
        ValueError: If a required seasonal year is missing.
    """
    series = rainfall.set_index("seasonal_year")[column]
    start, end = RAIN_STANDARDISE_YEARS
    # Standardising on early years only keeps later rainfall out of the fit.
    baseline = series.loc[start:end]
    needed = years - lag
    missing = set(needed) - set(series.index)
    if missing:
        raise ValueError(f"No rainfall for seasonal years {sorted(missing)}.")
    return ((series.loc[needed] - baseline.mean()) / baseline.std()).to_numpy()


def crps_normal(mean: np.ndarray, sd: np.ndarray, obs: np.ndarray) -> np.ndarray:
    """Continuous ranked probability score of a normal forecast.

    A zero sd is a point forecast, whose CRPS is the absolute error. This
    arises for the calibration-mean baseline in regions with only zero
    counts so far.

    Args:
        mean: Forecast means.
        sd: Forecast sds.
        obs: Observed values.

    Returns:
        CRPS for each forecast (lower is better).
    """
    shape = np.broadcast(mean, sd, obs).shape
    mean, sd, obs = np.broadcast_arrays(*map(np.atleast_1d, (mean, sd, obs)))
    crps = np.abs(obs - mean).astype(float)
    spread = sd > 0
    z = (obs[spread] - mean[spread]) / sd[spread]
    crps[spread] = sd[spread] * (
        z * (2.0 * norm.cdf(z) - 1.0) + 2.0 * norm.pdf(z) - 1.0 / np.sqrt(np.pi)
    )
    return crps.reshape(shape)


def _naive_forecasts(y: np.ndarray, years: np.ndarray, origin_index: int):
    """Random-walk and calibration-mean forecasts for the next year.

    Args:
        y: Log counts, shape ``(T, R)``.
        years: Model years.
        origin_index: Index of the origin year.

    Returns:
        Dict of model name to ``(mean, sd)`` arrays over regions.
    """
    past = y[: origin_index + 1]
    t_next = origin_index + 1
    steps = []
    for region in range(y.shape[1]):
        seen = np.flatnonzero(~np.isnan(past[:, region]))
        gaps = np.diff(years[seen])
        steps.extend(np.diff(past[seen, region]) / np.sqrt(gaps))
    # Drift-free random walk, so the sd is the root mean square step.
    step_sd = float(np.sqrt(np.mean(np.square(steps))))

    rw_mean, rw_sd = np.empty(y.shape[1]), np.empty(y.shape[1])
    for region in range(y.shape[1]):
        seen = np.flatnonzero(~np.isnan(past[:, region]))
        last = seen[-1] if len(seen) else 0
        rw_mean[region] = past[last, region] if len(seen) else np.nan
        rw_sd[region] = step_sd * np.sqrt(t_next - last)

    level_sd = np.nanstd(past, axis=0, ddof=1)
    level_sd = np.where(np.isfinite(level_sd), level_sd, np.nanstd(past, ddof=1))
    return {
        "random_walk": (rw_mean, rw_sd),
        "calibration_mean": (np.nanmean(past, axis=0), level_sd),
    }


def evaluate_species(
    regional: pd.DataFrame,
    rainfall: pd.DataFrame,
    species: str,
    rain_column: str = "dry_season_mm",
    rain_lag: int = 1,
    fixed_o=None,
) -> pd.DataFrame:
    """Score all forecasters for one species across the rolling origins.

    Args:
        regional: Regional counts long table.
        rainfall: Seasonal rainfall table.
        species: Species to evaluate.
        rain_column: Rainfall column used as the covariate.
        rain_lag: Seasonal-year lag of the covariate.
        fixed_o: Optional fixed observation sd (sensitivity runs).

    Returns:
        One row per origin, region and model with the forecast and scores.
    """
    series = prepare_series(regional, species)
    rain = lagged_rain(rainfall, series.years, rain_column, rain_lag)
    rows = []
    for origin in range(FIRST_ORIGIN, LAST_ORIGIN + 1):
        origin_index = int(np.flatnonzero(series.years == origin)[0])
        target = series.y[origin_index + 1]
        if np.isnan(target).all():
            continue
        window = series.y.copy()
        window[origin_index + 1 :] = np.nan
        masked = type(series)(series.years, series.regions, window, series.era)

        forecasts = _naive_forecasts(series.y, series.years, origin_index)
        for name, use_rain in (("rain", True), ("no_rain", False)):
            params = fit(masked, rain, use_rain=use_rain, fixed_o=fixed_o)
            forecasts[name] = forecast_next(params, series, rain, origin_index)

        for name, (mean, sd) in forecasts.items():
            for region_index, region in enumerate(series.regions):
                obs = target[region_index]
                if np.isnan(obs) or np.isnan(mean[region_index]):
                    continue
                m, s = mean[region_index], sd[region_index]
                rows.append(
                    {
                        "species": species,
                        "origin": origin,
                        "region": region,
                        "model": name,
                        "observed": obs,
                        "mean": m,
                        "sd": s,
                        "crps": float(crps_normal(m, s, obs)),
                        "abs_error": abs(m - obs),
                        "in_80": abs(obs - m) <= norm.ppf(0.9) * s,
                        "in_95": abs(obs - m) <= norm.ppf(0.975) * s,
                    }
                )
    return pd.DataFrame(rows)


def summarise(scores: pd.DataFrame) -> pd.DataFrame:
    """Mean scores per species and model, with skill against the random walk.

    Args:
        scores: Output of ``evaluate_species`` for one or more species.

    Returns:
        One row per species and model.
    """
    summary = (
        scores.groupby(["species", "model"])
        .agg(
            n=("crps", "size"),
            mean_crps=("crps", "mean"),
            mae=("abs_error", "mean"),
            coverage_80=("in_80", "mean"),
            coverage_95=("in_95", "mean"),
        )
        .reset_index()
    )
    random_walk = summary[summary["model"] == "random_walk"].set_index("species")
    summary["crps_skill_vs_random_walk"] = 1.0 - summary["mean_crps"] / summary[
        "species"
    ].map(random_walk["mean_crps"])
    return summary


@click.command()
@click.option(
    "--regional-path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=Path("data/raw/regional_counts.csv"),
    show_default=True,
)
@click.option(
    "--rainfall-path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=Path("data/raw/kruger_rainfall_cru.csv"),
    show_default=True,
)
@click.option(
    "--confirmation",
    is_flag=True,
    help="Score the pre-registered confirmation species (run once, frozen).",
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("results"),
    show_default=True,
)
def main(
    regional_path: Path, rainfall_path: Path, confirmation: bool, output_dir: Path
) -> None:
    """Run the primary rolling-origin evaluation and write the results.

    Args:
        regional_path: Regional counts CSV.
        rainfall_path: Seasonal rainfall CSV.
        confirmation: Whether to score the confirmation species.
        output_dir: Directory for the results CSVs.
    """
    regional = pd.read_csv(regional_path)
    rainfall = pd.read_csv(rainfall_path)
    species_set = CONFIRMATION_SPECIES if confirmation else DEVELOPMENT_SPECIES
    label = "confirmation" if confirmation else "development"
    scores = pd.concat(
        [evaluate_species(regional, rainfall, species) for species in species_set],
        ignore_index=True,
    )
    summary = summarise(scores)
    output_dir.mkdir(parents=True, exist_ok=True)
    scores.round(4).to_csv(output_dir / f"v2_{label}_scores.csv", index=False)
    summary.round(3).to_csv(output_dir / f"v2_{label}_summary.csv", index=False)
    click.echo(
        f"Version 2 {label} species, one-step rolling origin {FIRST_ORIGIN}-{LAST_ORIGIN}"
    )
    click.echo(summary.to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
