#!/usr/bin/env python3
"""Normalize the SANParks Kruger census summary for the model.

Writes two extracts from the same raw file: park totals per species and
year (used by version 1) and a long table of regional counts (used by the
version 2 state-space model).
"""

from pathlib import Path

import click
import pandas as pd

SOURCE_RECORD = "judithk.814.4"
SOURCE_OBJECT = "judithk.815.1"
SOURCE_URL = "https://cn.dataone.org/cn/v2/resolve/judithk.815.1"
SURVEY_METHOD = "fixed-wing aerial total count"
SPECIES_NAMES = {
    "Blue Wildebeest": "Wildebeest",
    "Buffalo": "Buffalo",
    "Eland": "Eland",
    "Giraffe": "Giraffe",
    "Impala": "Impala",
    "Kudu": "Kudu",
    "Roan": "Roan",
    "Sable": "Sable",
    "Tsessebe": "Tsessebe",
    "Warthog": "Warthog",
    "Waterbuck": "Waterbuck",
    "Zebra": "Zebra",
}
REGIONS = ("South", "Central", "North", "FarNorth")
# The DataONE abstract describes 1965-1976 as aerial totals corrected by
# ground counts and 1977-1997 as Ecological Aerial Survey totals.
LAST_CORRECTED_AERIAL_YEAR = 1976


def parse_counts(input_path: Path, output_path: Path) -> pd.DataFrame:
    """Read the raw SANParks summary and write normalized candidate counts.

    Args:
        input_path: Path to the downloaded SANParks tab-separated summary.
        output_path: Destination path for the normalized CSV extract.

    Returns:
        The normalized candidate-species counts, including missing totals.

    Raises:
        ValueError: If required columns are absent or a year is outside the
            comparable 1965-1997 total-count period.
    """
    source = pd.read_csv(input_path, sep="\t", na_values=[""], keep_default_na=True)
    required = {"year", "Species", "TOTAL"}
    missing = required.difference(source.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    selected = source[source["Species"].isin(SPECIES_NAMES)].copy()
    selected["name"] = selected["Species"].map(SPECIES_NAMES)
    selected["year"] = pd.to_numeric(selected["year"], errors="raise").astype(int)
    if not selected["year"].between(1965, 1997).all():
        raise ValueError("Census counts must be limited to 1965-1997.")

    counts = selected[["name", "year", "TOTAL"]].rename(columns={"TOTAL": "total_pop"})
    counts["total_pop"] = pd.to_numeric(counts["total_pop"], errors="coerce")
    counts["source_record"] = SOURCE_RECORD
    counts["source_object"] = SOURCE_OBJECT
    counts["source"] = SOURCE_URL
    counts["survey_method"] = SURVEY_METHOD
    counts = counts.sort_values(["name", "year"]).reset_index(drop=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_path, index=False, na_rep="")
    return counts


def parse_regional_counts(input_path: Path, output_path: Path) -> pd.DataFrame:
    """Write one row per species, year and region with a reported count.

    Each row also records the survey era and how the published park total
    relates to the regional counts for that species-year, so that years
    whose total includes estimates for unsurveyed regions stay visible.

    Args:
        input_path: Path to the downloaded SANParks tab-separated summary.
        output_path: Destination path for the regional long-table CSV.

    Returns:
        The regional counts, one row per non-missing regional cell.

    Raises:
        ValueError: If required columns are absent or a year is outside the
            comparable 1965-1997 total-count period.
    """
    source = pd.read_csv(input_path, sep="\t", na_values=[""], keep_default_na=True)
    required = {"year", "Species", "TOTAL", *REGIONS}
    missing = required.difference(source.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    selected = source[source["Species"].isin(SPECIES_NAMES)].copy()
    selected["name"] = selected["Species"].map(SPECIES_NAMES)
    selected["year"] = pd.to_numeric(selected["year"], errors="raise").astype(int)
    if not selected["year"].between(1965, 1997).all():
        raise ValueError("Census counts must be limited to 1965-1997.")

    selected["regions_reported"] = selected[list(REGIONS)].notna().sum(axis=1)
    selected["regional_sum"] = selected[list(REGIONS)].sum(axis=1, min_count=1)
    # A published total above the regional sum means unsurveyed regions were
    # estimated; a mismatch with all four present is a source inconsistency.
    selected["total_status"] = "no total"
    has_total = selected["TOTAL"].notna()
    matches = has_total & (selected["TOTAL"] == selected["regional_sum"])
    selected.loc[matches, "total_status"] = "equals regional sum"
    selected.loc[
        has_total & ~matches & (selected["regions_reported"] < 4), "total_status"
    ] = "includes unsurveyed-region estimate"
    selected.loc[
        has_total & ~matches & (selected["regions_reported"] == 4), "total_status"
    ] = "differs from regional sum"

    regional = selected.melt(
        id_vars=["name", "year", "total_status"],
        value_vars=list(REGIONS),
        var_name="region",
        value_name="count",
    ).dropna(subset=["count"])
    regional["count"] = regional["count"].astype(int)
    regional["era"] = regional["year"].map(
        lambda year: (
            "corrected aerial total"
            if year <= LAST_CORRECTED_AERIAL_YEAR
            else "ecological aerial survey"
        )
    )
    regional["source_record"] = SOURCE_RECORD
    regional["source_object"] = SOURCE_OBJECT
    regional["source"] = SOURCE_URL
    regional["region"] = pd.Categorical(regional["region"], categories=REGIONS)
    regional = regional.sort_values(["name", "region", "year"]).reset_index(drop=True)
    regional = regional[
        [
            "name",
            "year",
            "region",
            "count",
            "era",
            "total_status",
            "source_record",
            "source_object",
            "source",
        ]
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    regional.to_csv(output_path, index=False)
    return regional


@click.command()
@click.option(
    "--input-path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=Path("data/raw/KNPCountsSummary.txt"),
    show_default=True,
)
@click.option(
    "--output-path",
    type=click.Path(path_type=Path, dir_okay=False),
    default=Path("data/raw/species_counts.csv"),
    show_default=True,
)
@click.option(
    "--regional-output-path",
    type=click.Path(path_type=Path, dir_okay=False),
    default=Path("data/raw/regional_counts.csv"),
    show_default=True,
)
def main(input_path: Path, output_path: Path, regional_output_path: Path) -> None:
    """Build the normalized park-total and regional count CSVs.

    Args:
        input_path: Path to the downloaded SANParks summary.
        output_path: Destination path for the normalized counts CSV.
        regional_output_path: Destination path for the regional counts CSV.
    """
    counts = parse_counts(input_path, output_path)
    click.echo(
        f"Wrote {len(counts)} rows to {output_path}; "
        f"{counts['total_pop'].notna().sum()} have reported totals."
    )
    regional = parse_regional_counts(input_path, regional_output_path)
    click.echo(f"Wrote {len(regional)} regional rows to {regional_output_path}.")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
