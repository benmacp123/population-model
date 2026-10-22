#!/usr/bin/env python3
"""Version 3 design sweep: alert reliability by survey interval and precision.

For each rare-antelope species, fit the no-rain version 2 model to the
baseline window, then simulate surveys under every combination of
pre-registered design and scenario, filter them and score a threshold
alert on the core-region index. The design was fixed before the run (see
the README); nothing here scores against real counts.
"""

from itertools import product
from pathlib import Path

import click
import numpy as np
import pandas as pd

from monitoring import BURN_IN_YEARS, alert_metrics, core_regions
from monitoring import crossing_probability, filter_paths, reliability, simulate
from monitoring import survey_years, with_noise_split
from state_space import fit, prepare_series

SPECIES = ("Roan", "Sable", "Tsessebe")
# The first window is primary; 1979 is the sensitivity run.
BASELINE_ENDS = (1985, 1979)
OPERATING_MODELS = ("fitted_b", "near_random_walk")
NEAR_RANDOM_WALK_B = 0.99
ESTIMATORS = ("true", "as_fitted_split")
SHARES = (0.25, 0.5, 0.75)
INTERVALS = (1, 2, 3)
OBS_SDS = (0.1, 0.2, 0.4)
SCENARIOS = {
    "none": (0.0, 0.0),
    "decline_3": (0.03, 0.0),
    "decline_5": (0.05, 0.0),
    "decline_10": (0.10, 0.0),
    "drift_3": (0.0, 0.03),
}
THRESHOLDS = (0.7, 0.5)
ALPHAS = (0.1, 0.25, 0.5)
N_SIMS = 1000
BASE_SEED = 20261022


def fit_operating_model(regional: pd.DataFrame, species: str, baseline_end: int):
    """Fit the no-rain model to the baseline window.

    Args:
        regional: Regional counts long table.
        species: Species name.
        baseline_end: Last year used in the fit.

    Returns:
        ``(params, core, regions)``: fitted parameters, core-region indices
        and the region names.
    """
    series = prepare_series(regional, species)
    masked = series.y.copy()
    masked[series.years > baseline_end] = np.nan
    window = type(series)(series.years, series.regions, masked, series.era)
    params = fit(window, np.zeros(len(series.years)), use_rain=False)
    return params, core_regions(series, baseline_end), list(series.regions)


def _estimator_params(truth, fitted, estimator: str):
    """Parameters the simulated manager filters with.

    Args:
        truth: Operating-model parameters.
        fitted: Parameters as fitted to the baseline window.
        estimator: ``true`` to use the truth, or ``as_fitted_split`` to use
            the fitted process sd while the truth uses the swept split.

    Returns:
        Parameters for ``filter_paths``.
    """
    if estimator == "true":
        return truth
    return type(truth)(m=truth.m, b=truth.b, c=0.0, q=fitted.q, o=truth.o)


def run_sweep(regional: pd.DataFrame, n_sims: int):
    """Run every pre-registered cell.

    Args:
        regional: Regional counts long table.
        n_sims: Simulations per cell.

    Returns:
        ``(metrics, calibration, operating)`` tables.
    """
    metric_rows, operating_rows = [], []
    calibration_store = {}
    cell = 0
    for baseline_end, species in product(BASELINE_ENDS, SPECIES):
        fitted, core, regions = fit_operating_model(regional, species, baseline_end)
        operating_rows.append(
            {
                "baseline_end": baseline_end,
                "species": species,
                "core_regions": " ".join(regions[i] for i in core),
                "b": fitted.b,
                "q": fitted.q,
                "o_era1": fitted.o[0],
                "o_era2": fitted.o[1],
                "total_variance": fitted.q**2 + fitted.o[1] ** 2,
            }
        )
        for opmodel, share, interval, obs_sd, scenario in product(
            OPERATING_MODELS, SHARES, INTERVALS, OBS_SDS, SCENARIOS
        ):
            cell += 1
            truth = with_noise_split(fitted, share, obs_sd)
            if opmodel == "near_random_walk":
                truth.b = NEAR_RANDOM_WALK_B
            decline, drift = SCENARIOS[scenario]
            rng = np.random.default_rng(BASE_SEED + cell)
            true_centred, observed = simulate(
                truth, n_sims, rng, decline=decline, drift=drift
            )
            surveyed = survey_years(observed.shape[1], interval)
            observed[:, ~surveyed] = np.nan
            true_index = true_centred[..., core].mean(axis=-1)

            for estimator in ESTIMATORS:
                d, p = filter_paths(
                    _estimator_params(truth, fitted, estimator), observed
                )
                for threshold in THRESHOLDS:
                    log_threshold = float(np.log(threshold))
                    prob = crossing_probability(d, p, core, log_threshold)
                    primary = (
                        baseline_end == BASELINE_ENDS[0]
                        and opmodel == "fitted_b"
                        and estimator == "true"
                    )
                    if primary:
                        group = "model_true" if scenario == "none" else "misspecified"
                        keep = surveyed & (np.arange(len(surveyed)) >= BURN_IN_YEARS)
                        store = calibration_store.setdefault(
                            (species, threshold, group), ([], [])
                        )
                        store[0].append(prob[:, keep])
                        store[1].append(true_index[:, keep] < log_threshold)
                    for alpha in ALPHAS:
                        metric_rows.append(
                            {
                                "baseline_end": baseline_end,
                                "species": species,
                                "operating_model": opmodel,
                                "estimator": estimator,
                                "process_share": share,
                                "interval": interval,
                                "obs_sd": obs_sd,
                                "scenario": scenario,
                                "threshold": threshold,
                                "alpha": alpha,
                                **alert_metrics(
                                    prob,
                                    true_index,
                                    surveyed,
                                    alpha,
                                    log_threshold,
                                    declining=decline > 0,
                                ),
                            }
                        )

    calibration = []
    for (species, threshold, group), (probs, events) in calibration_store.items():
        table = reliability(
            np.concatenate([x.ravel() for x in probs]),
            np.concatenate([x.ravel() for x in events]),
        )
        table.insert(0, "group", group)
        table.insert(0, "threshold", threshold)
        table.insert(0, "species", species)
        calibration.append(table)
    return (
        pd.DataFrame(metric_rows),
        pd.concat(calibration, ignore_index=True),
        pd.DataFrame(operating_rows),
    )


def calibration_verdict(calibration: pd.DataFrame, min_n: int = 200):
    """Apply the pre-registered calibration rule to each group.

    Args:
        calibration: Output of ``run_sweep``.
        min_n: Smallest bin count that is judged.

    Returns:
        One row per species, threshold and group with the largest gap
        between observed frequency and bin mid-point, and the verdict.
    """
    judged = calibration[calibration["n"] >= min_n].copy()
    judged["gap"] = (judged["observed_freq"] - judged["bin_mid"]).abs()
    verdict = (
        judged.groupby(["species", "threshold", "group"])["gap"].max().reset_index()
    )
    verdict["calibrated"] = verdict["gap"] <= 0.05
    return verdict


@click.command()
@click.option(
    "--regional-path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=Path("data/raw/regional_counts.csv"),
    show_default=True,
)
@click.option("--n-sims", type=int, default=N_SIMS, show_default=True)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("results/v3"),
    show_default=True,
)
def main(regional_path: Path, n_sims: int, output_dir: Path) -> None:
    """Run the version 3 design sweep and write the results.

    Args:
        regional_path: Regional counts CSV.
        n_sims: Simulations per cell.
        output_dir: Directory for the results CSVs.
    """
    regional = pd.read_csv(regional_path)
    metrics, calibration, operating = run_sweep(regional, n_sims)
    verdict = calibration_verdict(calibration)
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics.round(4).to_csv(output_dir / "survey_design.csv", index=False)
    calibration.round(4).to_csv(output_dir / "calibration.csv", index=False)
    verdict.round(4).to_csv(output_dir / "calibration_verdict.csv", index=False)
    operating.round(4).to_csv(output_dir / "operating_models.csv", index=False)
    click.echo(operating.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    click.echo(verdict.to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
