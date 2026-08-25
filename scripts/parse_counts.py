#!/usr/bin/env python3
"""Normalize the SANParks Kruger census summary for the model."""

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

    counts = selected[["name", "year", "TOTAL"]].rename(
        columns={"TOTAL": "total_pop"}
    )
    counts["total_pop"] = pd.to_numeric(counts["total_pop"], errors="coerce")
    counts["source_record"] = SOURCE_RECORD
    counts["source_object"] = SOURCE_OBJECT
    counts["source"] = SOURCE_URL
    counts["survey_method"] = SURVEY_METHOD
    counts = counts.sort_values(["name", "year"]).reset_index(drop=True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_path, index=False, na_rep="")
    return counts


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
def main(input_path: Path, output_path: Path) -> None:
    """Build the normalized census-count CSV from the command line.

    Args:
        input_path: Path to the downloaded SANParks summary.
        output_path: Destination path for the normalized counts CSV.
    """
    counts = parse_counts(input_path, output_path)
    click.echo(
        f"Wrote {len(counts)} rows to {output_path}; "
        f"{counts['total_pop'].notna().sum()} have reported totals."
    )


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
