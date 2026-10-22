"""Tests for normalising the SANParks census summary."""

import pandas as pd
import pytest

from parse_counts import parse_regional_counts

HEADER = "year\tSpecies\tSouth\tCentral\tNorth\tFarNorth\tTOTAL\n"


@pytest.fixture
def raw_path(tmp_path):
    """Tiny census summary covering each total-status case."""
    path = tmp_path / "summary.txt"
    path.write_text(
        HEADER
        + "1975\tKudu\t10\t20\t30\t40\t100\n"
        + "1979\tKudu\t\t20\t30\t40\t120\n"
        + "1980\tKudu\t10\t20\t30\t40\t101\n"
        + "1981\tKudu\t10\t\t\t\t\n",
        encoding="utf-8",
    )
    return path


def test_regional_rows_only_for_reported_cells(raw_path, tmp_path):
    """Missing regional cells are dropped, never written as zero."""
    regional = parse_regional_counts(raw_path, tmp_path / "out.csv")
    assert len(regional) == 4 + 3 + 4 + 1
    assert (regional["count"] > 0).all()


def test_total_status_flags_infill_and_mismatch(raw_path, tmp_path):
    """Totals above a partial regional sum are flagged as estimates."""
    regional = parse_regional_counts(raw_path, tmp_path / "out.csv")
    status = regional.groupby("year")["total_status"].first()
    assert status[1975] == "equals regional sum"
    assert status[1979] == "includes unsurveyed-region estimate"
    assert status[1980] == "differs from regional sum"
    assert status[1981] == "no total"


def test_era_splits_at_1976(raw_path, tmp_path):
    """Counts up to 1976 and from 1977 carry different survey eras."""
    regional = parse_regional_counts(raw_path, tmp_path / "out.csv")
    eras = regional.groupby("year")["era"].first()
    assert eras[1975] == "corrected aerial total"
    assert eras[1979] == "ecological aerial survey"


def test_rejects_years_outside_total_count_period(tmp_path):
    """Distance-sampling years after 1997 must not enter the series."""
    path = tmp_path / "summary.txt"
    path.write_text(HEADER + "1998\tKudu\t1\t1\t1\t1\t4\n", encoding="utf-8")
    with pytest.raises(ValueError, match="1965-1997"):
        parse_regional_counts(path, tmp_path / "out.csv")


def test_committed_regional_extract_matches_parser(tmp_path):
    """The committed regional CSV is reproducible from the raw file."""
    rebuilt = parse_regional_counts(
        "data/raw/KNPCountsSummary.txt", tmp_path / "out.csv"
    )
    committed = pd.read_csv("data/raw/regional_counts.csv")
    pd.testing.assert_frame_equal(
        rebuilt.astype({"region": str}).reset_index(drop=True), committed
    )
