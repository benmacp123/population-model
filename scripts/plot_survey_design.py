#!/usr/bin/env python3
"""Plot decline detection against false alarms for the version 3 sweep.

Reads ``results/v3/survey_design.csv`` written by
``scripts/survey_design.py`` and writes one small-multiple figure: species
in rows, survey interval in columns, one line per survey precision whose
three points are the alert levels alpha.
"""

from pathlib import Path

import click
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # pylint: disable=wrong-import-position
import pandas as pd  # pylint: disable=wrong-import-position

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e4e3df"
# Ordinal blue ramp, darkest for the most precise survey; validated with
# the dataviz skill's palette checker. Marker shape repeats the encoding.
PRECISION_STYLE = {
    0.1: ("#104281", "o"),
    0.2: ("#2a78d6", "s"),
    0.4: ("#86b6ef", "^"),
}
DECLINE_SCENARIO = "decline_5"
# The primary cells reported in the README's design-sweep results.
PRIMARY = {
    "baseline_end": 1985,
    "operating_model": "fitted_b",
    "estimator": "true",
    "process_share": 0.5,
    "threshold": 0.7,
}


def primary_cells(design: pd.DataFrame) -> pd.DataFrame:
    """Pair false-alarm and detection metrics for the primary cells.

    Args:
        design: Full sweep table.

    Returns:
        One row per species, interval, precision and alpha with the chance
        of any false alarm in 10 stable years, the share of declines
        detected within the miss window and the median delay.
    """
    rows = design
    for column, value in PRIMARY.items():
        rows = rows[rows[column] == value]
    keys = ["species", "interval", "obs_sd", "alpha"]
    stable = rows[rows["scenario"] == "none"][keys + ["any_false_alarm_10y"]]
    decline = rows[rows["scenario"] == DECLINE_SCENARIO]
    decline = decline[keys + ["miss_rate", "median_delay_years"]]
    paired = stable.merge(decline, on=keys)
    paired["detected"] = 1.0 - paired["miss_rate"]
    return paired.sort_values(keys).reset_index(drop=True)


def _style_axes(ax) -> None:
    """Apply the recessive grid and axis styling shared by every panel.

    Args:
        ax: Matplotlib axes.
    """
    ax.set_facecolor(SURFACE)
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xticks([0, 0.5, 1])
    ax.set_yticks([0, 0.5, 1])
    ax.grid(color=GRID, lw=0.8)
    ax.tick_params(colors=INK_SECONDARY, labelsize=8, length=0)
    percent = matplotlib.ticker.PercentFormatter(xmax=1, decimals=0)
    ax.xaxis.set_major_formatter(percent)
    ax.yaxis.set_major_formatter(percent)
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_survey_design(paired: pd.DataFrame, path: Path) -> None:
    """Draw the species-by-interval grid and save it.

    Args:
        paired: Output of ``primary_cells``.
        path: Output image path.
    """
    species_list = list(dict.fromkeys(paired["species"]))
    intervals = sorted(paired["interval"].unique())
    fig, axes = plt.subplots(
        len(species_list),
        len(intervals),
        figsize=(3.0 * len(intervals) + 0.8, 2.7 * len(species_list) + 1.0),
        dpi=150,
        sharex=True,
        sharey=True,
    )
    fig.patch.set_facecolor(SURFACE)

    for row, species in enumerate(species_list):
        for col, interval in enumerate(intervals):
            ax = axes[row, col]
            _style_axes(ax)
            cell = paired[
                (paired["species"] == species) & (paired["interval"] == interval)
            ]
            for obs_sd, (color, marker) in PRECISION_STYLE.items():
                line = cell[cell["obs_sd"] == obs_sd].sort_values("alpha")
                ax.plot(
                    line["any_false_alarm_10y"],
                    line["detected"],
                    color=color,
                    lw=2,
                    marker=marker,
                    ms=6.5,
                    mec=SURFACE,
                    mew=1.2,
                    label=f"{obs_sd:.1f}",
                    zorder=3,
                )
            if row == 0:
                years = "year" if interval == 1 else "years"
                ax.set_title(
                    f"Survey every {interval} {years}",
                    loc="left",
                    color=INK,
                    fontsize=10,
                    fontweight="bold",
                )
            # Species names sit in the empty lower-right corner of each row's
            # first panel, keeping the left margin for the shared axis label.
            if col == 0:
                ax.text(
                    0.97,
                    0.05,
                    species,
                    transform=ax.transAxes,
                    ha="right",
                    color=INK,
                    fontsize=11,
                    fontweight="bold",
                )

    # One direct note on how to read a line, on the panel where it is clearest.
    example = paired[
        (paired["species"] == species_list[-1])
        & (paired["interval"] == intervals[0])
        & (paired["obs_sd"] == max(PRECISION_STYLE))
    ].sort_values("alpha")
    for _, point in example.iloc[[0, -1]].iterrows():
        axes[-1, 0].annotate(
            f"α {point['alpha']:g}",
            (point["any_false_alarm_10y"], point["detected"]),
            xytext=(6, -4),
            textcoords="offset points",
            color=INK_SECONDARY,
            fontsize=7.5,
        )

    fig.supxlabel(
        "Chance of at least one false alarm in 10 stable years",
        color=INK_SECONDARY,
        fontsize=9,
        y=0.065,
    )
    fig.supylabel(
        "Share of 5%-a-year declines alerted within 5 years of crossing",
        color=INK_SECONDARY,
        fontsize=9,
    )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        title="Survey observation sd (log scale)",
        title_fontsize=8,
        loc="lower center",
        ncol=len(labels),
        frameon=False,
        fontsize=8,
        labelcolor=INK_SECONDARY,
        bbox_to_anchor=(0.5, 0.0),
    )
    fig.tight_layout(rect=(0.01, 0.07, 1, 1))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)


@click.command()
@click.option(
    "--design-path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=Path("results/v3/survey_design.csv"),
    show_default=True,
)
@click.option(
    "--output-path",
    type=click.Path(path_type=Path, dir_okay=False),
    default=Path("results/v3/survey_design.png"),
    show_default=True,
)
def main(design_path: Path, output_path: Path) -> None:
    """Write the version 3 design figure.

    Args:
        design_path: Sweep results CSV.
        output_path: Output PNG path.
    """
    plot_survey_design(primary_cells(pd.read_csv(design_path)), output_path)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
