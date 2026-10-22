#!/usr/bin/env python3
"""Version 3 phase 5: replay a fixed alert rule over 1986-1997.

For each species, fit the no-rain version 2 model to the 1969-1985
baseline, choose alpha by the rule fixed before the run (the smallest
grid value with at most 0.05 false alarms per survey-year when the fitted
model is simulated with annual surveys), then run the filter forward
through the real counts and record when the alert would first have fired.
A raw-count rule (the observed core-region index below the threshold) is
reported alongside. The outcome is already known, so results are
descriptive, not a test.
"""

from pathlib import Path

import click
import numpy as np
import pandas as pd

from monitoring import alert_metrics, core_regions, crossing_probability
from monitoring import filter_paths, simulate, survey_years
from state_space import SpeciesSeries, fit, kalman_filter, prepare_series

SPECIES = ("Roan", "Sable", "Tsessebe")
CONTROLS = ("Giraffe", "Impala", "Zebra")
BASELINE = (1969, 1985)
REPLAY = (1986, 1997)
THRESHOLD = 0.7
ALPHAS = (0.1, 0.25, 0.5)
MAX_FALSE_ALARM = 0.05
N_SIMS = 1000
# Offset from the sweep's seeds so no simulated path is shared with it.
SEED = 20261022 + 100_000


def fit_baseline(regional: pd.DataFrame, species: str):
    """Fit the no-rain model and pick core regions from the baseline only.

    Args:
        regional: Regional counts long table.
        species: Species name.

    Returns:
        ``(series, params, core)``: the full series, the fitted parameters
        and the core-region indices.
    """
    series = prepare_series(regional, species)
    outside = (series.years < BASELINE[0]) | (series.years > BASELINE[1])
    masked = series.y.copy()
    masked[outside] = np.nan
    window = SpeciesSeries(series.years, series.regions, masked, series.era)
    params = fit(window, np.zeros(len(series.years)), use_rain=False)
    return series, params, core_regions(window, BASELINE[1])


def false_alarm_rates(params, core, rng: np.random.Generator) -> pd.DataFrame:
    """False alarms per survey-year for each alpha under the fitted model.

    The fitted parameters are both the truth and the estimator, surveys
    are annual and observations use the fitted era-2 sd.

    Args:
        params: Fitted no-rain parameters.
        core: Core-region indices.
        rng: Random generator.

    Returns:
        One row per alpha with its false-alarm rate per survey-year.
    """
    true_centred, observed = simulate(params, N_SIMS, rng)
    surveyed = survey_years(observed.shape[1], 1)
    d, p = filter_paths(params, observed)
    log_threshold = float(np.log(THRESHOLD))
    prob = crossing_probability(d, p, core, log_threshold)
    true_index = true_centred[..., core].mean(axis=-1)
    rows = []
    for alpha in ALPHAS:
        metrics = alert_metrics(
            prob, true_index, surveyed, alpha, log_threshold, declining=False
        )
        rows.append(
            {
                "alpha": alpha,
                "false_alarm_per_survey_year": metrics["false_alarm_per_survey_year"],
            }
        )
    return pd.DataFrame(rows)


def choose_alpha(rates: pd.DataFrame):
    """Apply the pre-registered alpha rule.

    Args:
        rates: Output of ``false_alarm_rates``.

    Returns:
        ``(alpha, rule_met)``. If no grid value meets the bound, the
        largest alpha is returned with ``rule_met`` false; this fallback
        was fixed before the replay was run.
    """
    passing = rates[rates["false_alarm_per_survey_year"] <= MAX_FALSE_ALARM]
    if passing.empty:
        return float(max(ALPHAS)), False
    return float(passing["alpha"].min()), True


def replay_years(series: SpeciesSeries, params, core, alpha: float) -> pd.DataFrame:
    """Crossing probability and both alerts for every replay year.

    Each year's filtered state uses only counts up to that year, so the
    replay is what a manager would have seen at the time. A year counts as
    surveyed for the filter if any core region was counted; the raw index
    needs every core region, so it is missing otherwise.

    Args:
        series: Full observations.
        params: Baseline-fitted parameters.
        core: Core-region indices.
        alpha: Alert level.

    Returns:
        One row per replay year.
    """
    log_threshold = float(np.log(THRESHOLD))
    rain = np.zeros(len(series.years))
    rows = []
    for t in np.flatnonzero((series.years >= REPLAY[0]) & (series.years <= REPLAY[1])):
        window = slice(0, t + 1)
        _, d, p = kalman_filter(
            params, series.y[window], rain[window], series.era[window]
        )
        prob = float(crossing_probability(d, p, core, log_threshold))
        counts = series.y[t, core]
        surveyed = bool(np.any(~np.isnan(counts)))
        complete = bool(np.all(~np.isnan(counts)))
        raw_index = (
            float(np.mean(counts - params.m[core])) if complete else float("nan")
        )
        rows.append(
            {
                "year": int(series.years[t]),
                "crossing_probability": prob,
                "filter_alert": surveyed and prob > alpha,
                "raw_index_ratio": float(np.exp(raw_index)),
                "raw_alert": complete and raw_index < log_threshold,
            }
        )
    return pd.DataFrame(rows)


def _first_year(years: pd.DataFrame, column: str):
    """First year in which an alert column is true.

    Args:
        years: Output of ``replay_years``.
        column: Alert column name.

    Returns:
        The year, or ``pd.NA`` if the alert never fires.
    """
    fired = years.loc[years[column], "year"]
    return int(fired.iloc[0]) if len(fired) else pd.NA


def run_replay(regional: pd.DataFrame):
    """Replay every species and negative control.

    Args:
        regional: Regional counts long table.

    Returns:
        ``(summary, alphas, years)`` tables.
    """
    summary, alpha_tables, year_tables = [], [], []
    for index, species in enumerate(SPECIES + CONTROLS):
        series, params, core = fit_baseline(regional, species)
        rates = false_alarm_rates(params, core, np.random.default_rng(SEED + index))
        alpha, rule_met = choose_alpha(rates)
        years = replay_years(series, params, core, alpha)
        rates.insert(0, "species", species)
        years.insert(0, "species", species)
        alpha_tables.append(rates)
        year_tables.append(years)
        summary.append(
            {
                "species": species,
                "role": "target" if species in SPECIES else "control",
                "core_regions": " ".join(series.regions[i] for i in core),
                "total_variance": params.q**2 + params.o[1] ** 2,
                "alpha": alpha,
                "alpha_rule_met": rule_met,
                "first_filter_alert": _first_year(years, "filter_alert"),
                "first_raw_alert": _first_year(years, "raw_alert"),
                "filter_alert_years": int(years["filter_alert"].sum()),
                "raw_alert_years": int(years["raw_alert"].sum()),
                "raw_index_years": int(years["raw_index_ratio"].notna().sum()),
            }
        )
    return (
        pd.DataFrame(summary),
        pd.concat(alpha_tables, ignore_index=True),
        pd.concat(year_tables, ignore_index=True),
    )


@click.command()
@click.option(
    "--regional-path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=Path("data/raw/regional_counts.csv"),
    show_default=True,
)
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=Path("results/v3"),
    show_default=True,
)
def main(regional_path: Path, output_dir: Path) -> None:
    """Run the replay and write the results.

    Args:
        regional_path: Regional counts CSV.
        output_dir: Directory for the results CSVs.
    """
    summary, alphas, years = run_replay(pd.read_csv(regional_path))
    output_dir.mkdir(parents=True, exist_ok=True)
    summary.round(4).to_csv(output_dir / "replay_summary.csv", index=False)
    alphas.round(4).to_csv(output_dir / "replay_alpha.csv", index=False)
    years.round(4).to_csv(output_dir / "replay_years.csv", index=False)
    click.echo(summary.to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
