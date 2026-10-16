#!/usr/bin/env python3
"""Plot observed totals against calibrated and holdout projections.

Reads ``results/backtest_projections.csv`` written by
``scripts/evaluate_model.py`` and writes one small-multiple figure.
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
MODEL = "#2a78d6"


def plot_backtest(projections: pd.DataFrame, summary: pd.DataFrame, path: Path):
    """Draw one panel per species and save the figure.

    Args:
        projections: Year-by-year table from ``evaluate_model``.
        summary: Per-species summary from ``evaluate_model``.
        path: Output image path.
    """
    species_list = list(summary["species"])
    fig, axes = plt.subplots(
        1, len(species_list), figsize=(4.2 * len(species_list), 3.6), dpi=150
    )
    fig.patch.set_facecolor(SURFACE)

    for ax, species in zip(axes, species_list):
        rows = projections[projections["species"] == species]
        observed = rows[rows["period"] == "observed"]
        fitted = rows[rows["period"] == "calibration"]
        holdout = rows[rows["period"] == "holdout"]
        stats = summary[summary["species"] == species].iloc[0]
        anchor_year = holdout["year"].min()
        anchor_total = holdout["central"].iloc[0]

        ax.set_facecolor(SURFACE)
        ax.axvspan(
            fitted["year"].min(), anchor_year, color=GRID, alpha=0.5, lw=0, zorder=0
        )
        ax.fill_between(
            holdout["year"],
            holdout["low"],
            holdout["high"],
            color=MODEL,
            alpha=0.18,
            lw=0,
            label="Scenario range",
        )
        ax.plot(
            fitted["year"],
            fitted["central"],
            color=MODEL,
            lw=2,
            ls=(0, (1, 1.5)),
            label="Calibration fit",
        )
        ax.plot(holdout["year"], holdout["central"], color=MODEL, lw=2, label="Model")
        ax.hlines(
            anchor_total,
            anchor_year,
            holdout["year"].max(),
            color=INK_SECONDARY,
            lw=1.5,
            ls=(0, (4, 3)),
            label="Persistence baseline",
        )
        ax.plot(
            observed["year"],
            observed["central"],
            "o",
            ms=4.5,
            color=INK,
            mec=SURFACE,
            mew=1,
            label="Observed",
            zorder=3,
        )

        ax.set_title(species, loc="left", color=INK, fontsize=11, fontweight="bold")
        ax.set_title(
            f"holdout skill {stats['skill_vs_baseline']:+.2f}",
            loc="right",
            color=INK_SECONDARY,
            fontsize=8,
        )
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", color=GRID, lw=0.8)
        ax.tick_params(colors=INK_SECONDARY, labelsize=8, length=0)
        ax.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(
                lambda v, _: f"{v / 1000:.0f}k" if v else "0"
            )
        )
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color(INK_SECONDARY)

    axes[0].set_ylabel("Census total", color=INK_SECONDARY, fontsize=9)
    axes[0].text(
        0.02,
        0.04,
        "calibration",
        transform=axes[0].transAxes,
        color=INK_SECONDARY,
        fontsize=8,
    )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=5,
        frameon=False,
        fontsize=8,
        labelcolor=INK_SECONDARY,
    )
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, facecolor=SURFACE)
    plt.close(fig)


@click.command()
@click.option(
    "--results-dir",
    type=click.Path(path_type=Path, file_okay=False, exists=True),
    default=Path("results"),
    show_default=True,
)
def main(results_dir: Path) -> None:
    """Write results/backtest.png from the saved evaluation tables."""
    projections = pd.read_csv(results_dir / "backtest_projections.csv")
    summary = pd.read_csv(results_dir / "backtest_summary.csv")
    output = results_dir / "backtest.png"
    plot_backtest(projections, summary, output)
    click.echo(f"Wrote {output}")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
