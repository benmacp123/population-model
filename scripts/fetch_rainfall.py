#!/usr/bin/env python3
"""Fetch CRU TS precipitation and write a Kruger seasonal rainfall table.

This is an optional fetch step, kept separate from the analysis: the output
``data/raw/kruger_rainfall_cru.csv`` is committed, so the model reproduces
without network access. It downloads the decadal CRU TS 4.10 ASCII grids
for 1961-2000, averages the 0.5-degree cells covering Kruger and sums the
months into seasonal years.

A seasonal year is labelled by the calendar year it ends in: seasonal year
``y`` runs from October ``y - 1`` to September ``y``. Its dry season is
April to September of ``y``, which ends around that year's census.
"""

import gzip
import tempfile
import urllib.request
from pathlib import Path

import click
import numpy as np
import pandas as pd

CRU_VERSION = "4.10"
CRU_BASE_URL = (
    "https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.10/cruts.2604091129.v4.10/pre"
)
DECADES = ("1961.1970", "1971.1980", "1981.1990", "1991.2000")
FIRST_MONTH_YEAR = 1961
GRID_ROWS, GRID_COLS = 360, 720
# Cell centres whose 0.5-degree squares lie mostly inside the park, which
# spans roughly 22.3-25.5 S and 31.0-32.0 E. The match is approximate at
# the park's northern tip and western boundary.
KRUGER_LATS = (-22.25, -22.75, -23.25, -23.75, -24.25, -24.75, -25.25)
KRUGER_LONS = (31.25, 31.75)


def cell_index(lat: float, lon: float):
    """Row and column of a CRU ASCII grid cell centre.

    Row 0 is the southernmost row (89.75 S), confirmed by Kruger receiving
    about 520 mm a year where the mirrored row (Aswan) receives about 6 mm.

    Args:
        lat: Cell-centre latitude in degrees (south negative).
        lon: Cell-centre longitude in degrees (west negative).

    Returns:
        A ``(row, column)`` pair of integer indices.
    """
    return int(round((lat + 89.75) / 0.5)), int(round((lon + 179.75) / 0.5))


def read_grid(path: Path, scale: float) -> np.ndarray:
    """Read a gzipped CRU ASCII file into a months x rows x columns array.

    Args:
        path: Path to a ``.dat.gz`` or ``.stn.gz`` CRU ASCII file.
        scale: Divisor converting stored integers to actual values.

    Returns:
        Float array with missing values (``-999``) as NaN.
    """
    with gzip.open(path, "rt") as handle:
        values = np.loadtxt(handle, dtype=np.int32)
    grid = values.reshape(-1, GRID_ROWS, GRID_COLS).astype(float)
    grid[grid == -999] = np.nan
    return grid / scale


def kruger_cells(grid: np.ndarray) -> np.ndarray:
    """Select the Kruger cells from a gridded array.

    Args:
        grid: Months x rows x columns array.

    Returns:
        Months x cells array for the cells in ``KRUGER_LATS`` x ``KRUGER_LONS``.
    """
    cells = [cell_index(lat, lon) for lat in KRUGER_LATS for lon in KRUGER_LONS]
    return np.stack([grid[:, row, col] for row, col in cells], axis=1)


def seasonal_rainfall(
    monthly_mm: np.ndarray, stations: np.ndarray, first_year: int
) -> pd.DataFrame:
    """Sum park-average monthly rainfall into seasonal years.

    Args:
        monthly_mm: Months x cells precipitation in millimetres, starting in
            January of ``first_year``.
        stations: Months x cells contributing station counts, same shape.
        first_year: Calendar year of the first month.

    Returns:
        One row per complete seasonal year with wet-season, dry-season and
        annual totals in millimetres and the fewest contributing stations
        in any cell-month of that year.

    Raises:
        ValueError: If the inputs differ in shape or cover partial years.
    """
    if monthly_mm.shape != stations.shape:
        raise ValueError("Rainfall and station arrays must have the same shape.")
    if monthly_mm.shape[0] % 12:
        raise ValueError("Monthly input must cover whole calendar years.")

    park = monthly_mm.mean(axis=1)
    rows = []
    for start in range(9, len(park) - 11, 12):
        # Index ``start`` is October of the year before the seasonal year.
        window = slice(start, start + 12)
        year = first_year + (start + 3) // 12
        months = park[window]
        rows.append(
            {
                "seasonal_year": year,
                "wet_season_mm": float(months[:6].sum()),
                "dry_season_mm": float(months[6:].sum()),
                "annual_mm": float(months.sum()),
                "min_station_count": int(np.nanmin(stations[window])),
            }
        )
    return pd.DataFrame(rows)


def _download(name: str, cache_dir: Path) -> Path:
    """Download one CRU file into the cache unless it is already there.

    Args:
        name: File name under ``CRU_BASE_URL``.
        cache_dir: Directory holding downloaded files.

    Returns:
        Path to the local copy.
    """
    path = cache_dir / name
    if not path.exists():
        urllib.request.urlretrieve(f"{CRU_BASE_URL}/{name}", path)
    return path


def fetch_rainfall(cache_dir: Path, output_path: Path) -> pd.DataFrame:
    """Download the CRU decades, extract Kruger and write the seasonal table.

    Args:
        cache_dir: Directory for the downloaded CRU files.
        output_path: Destination CSV path.

    Returns:
        The seasonal rainfall table that was written.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    rain, stations = [], []
    for decade in DECADES:
        prefix = f"cru_ts{CRU_VERSION}.{decade}.pre"
        rain.append(
            kruger_cells(read_grid(_download(f"{prefix}.dat.gz", cache_dir), 10.0))
        )
        stations.append(
            kruger_cells(read_grid(_download(f"{prefix}.stn.gz", cache_dir), 1.0))
        )
    table = seasonal_rainfall(
        np.concatenate(rain), np.concatenate(stations), FIRST_MONTH_YEAR
    )
    table["source"] = (
        f"CRU TS {CRU_VERSION} pre, {len(KRUGER_LATS) * len(KRUGER_LONS)} cells"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.round(1).to_csv(output_path, index=False)
    return table


@click.command()
@click.option(
    "--cache-dir",
    type=click.Path(path_type=Path, file_okay=False),
    default=None,
    help="Keep downloads here; defaults to a temporary directory.",
)
@click.option(
    "--output-path",
    type=click.Path(path_type=Path, dir_okay=False),
    default=Path("data/raw/kruger_rainfall_cru.csv"),
    show_default=True,
)
def main(cache_dir: Path, output_path: Path) -> None:
    """Fetch CRU TS rainfall and write the Kruger seasonal table.

    Args:
        cache_dir: Optional directory for the downloaded CRU files.
        output_path: Destination CSV path.
    """
    if cache_dir is None:
        with tempfile.TemporaryDirectory() as tmp:
            table = fetch_rainfall(Path(tmp), output_path)
    else:
        table = fetch_rainfall(cache_dir, output_path)
    click.echo(f"Wrote {len(table)} seasonal years to {output_path}.")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
