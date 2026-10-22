"""Tests for the CRU TS rainfall extraction helpers."""

import numpy as np
import pandas as pd
import pytest

from fetch_rainfall import cell_index, seasonal_rainfall


def test_cell_index_counts_rows_from_the_south():
    """Row 0 is 89.75 S and column 0 is 179.75 W."""
    assert cell_index(-89.75, -179.75) == (0, 0)
    assert cell_index(-24.75, 31.25) == (130, 422)


def test_seasonal_year_runs_october_to_september():
    """Seasonal year y sums October y-1 to September y, split wet and dry."""
    months = np.arange(36, dtype=float)
    rain = np.column_stack([months, months])
    table = seasonal_rainfall(rain, np.full_like(rain, 8.0), first_year=1961)
    # Two complete seasonal years fit in three calendar years.
    assert table["seasonal_year"].tolist() == [1962, 1963]
    first = table.iloc[0]
    assert first["wet_season_mm"] == pytest.approx(sum(range(9, 15)))
    assert first["dry_season_mm"] == pytest.approx(sum(range(15, 21)))
    assert first["annual_mm"] == pytest.approx(sum(range(9, 21)))


def test_seasonal_rainfall_reports_fewest_stations():
    """The station column exposes any cell-month that fell back to climatology."""
    rain = np.ones((24, 2))
    stations = np.full((24, 2), 8.0)
    stations[12, 1] = 0.0
    table = seasonal_rainfall(rain, stations, first_year=1961)
    assert table.iloc[0]["min_station_count"] == 0


def test_seasonal_rainfall_rejects_partial_years():
    """Input must start in January and cover whole years."""
    with pytest.raises(ValueError, match="whole calendar years"):
        seasonal_rainfall(np.ones((18, 2)), np.ones((18, 2)), 1961)


def test_committed_rainfall_marks_known_droughts():
    """The committed table shows the 1982/83 and 1991/92 droughts as driest."""
    table = pd.read_csv("data/raw/kruger_rainfall_cru.csv")
    driest = set(table.nsmallest(2, "annual_mm")["seasonal_year"])
    assert driest == {1983, 1992}
    assert table["seasonal_year"].min() <= 1964
