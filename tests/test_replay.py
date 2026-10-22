"""Tests for the version 3 retrospective replay."""

import numpy as np
import pandas as pd
import pytest

from replay import ALPHAS, REPLAY, choose_alpha, fit_baseline, replay_years
from state_space import Params, SpeciesSeries


def _rates(values):
    """False-alarm table with one rate per grid alpha.

    Args:
        values: False-alarm rates, one per entry of ``ALPHAS``.

    Returns:
        A table shaped like ``false_alarm_rates`` output.
    """
    return pd.DataFrame({"alpha": ALPHAS, "false_alarm_per_survey_year": values})


def test_choose_alpha_takes_smallest_qualifying_value():
    """The most sensitive alpha within the false-alarm bound is chosen."""
    assert choose_alpha(_rates([0.08, 0.04, 0.01])) == (0.25, True)


def test_choose_alpha_falls_back_to_largest_when_none_qualify():
    """With no qualifying alpha, the fixed fallback is flagged."""
    assert choose_alpha(_rates([0.3, 0.2, 0.1])) == (max(ALPHAS), False)


def _series():
    """Two-region series spanning the replay with a gap in one region.

    Returns:
        A series of 101 animals per region that halves from 1990, with
        region 1 uncounted in 1992.
    """
    years = np.arange(1980, REPLAY[1] + 1)
    y = np.full((len(years), 2), np.log(101.0))
    # A 50% drop from 1990 on, and region 1 uncounted in 1992.
    y[years >= 1990] = np.log(51.0)
    y[years == 1992, 1] = np.nan
    return SpeciesSeries(years, ["A", "B"], y, np.ones(len(years), dtype=int))


def test_replay_raw_rule_needs_every_core_region():
    """A year missing a core region has no raw index or raw alert."""
    params = Params(
        m=np.full(2, np.log(101.0)), b=0.7, c=0.0, q=0.1, o=np.array([0.1, 0.1])
    )
    years = replay_years(_series(), params, [0, 1], alpha=0.5)
    gap = years[years["year"] == 1992].iloc[0]
    assert np.isnan(gap["raw_index_ratio"])
    assert not gap["raw_alert"]
    assert gap["filter_alert"]


def test_replay_alerts_fire_only_after_the_drop():
    """Both rules stay quiet before a halving and fire once it arrives."""
    params = Params(
        m=np.full(2, np.log(101.0)), b=0.7, c=0.0, q=0.1, o=np.array([0.1, 0.1])
    )
    years = replay_years(_series(), params, [0, 1], alpha=0.5)
    before = years[years["year"] < 1990]
    assert not before["filter_alert"].any()
    assert not before["raw_alert"].any()
    assert years.loc[years["year"] == 1990, "raw_alert"].iloc[0]
    assert years["year"].tolist() == list(range(REPLAY[0], REPLAY[1] + 1))


@pytest.mark.parametrize("species", ["Roan", "Sable", "Tsessebe"])
def test_baseline_fit_matches_design_sweep(species):
    """Masking pre-1969 years leaves the sweep's operating models unchanged."""
    regional = pd.read_csv("data/raw/regional_counts.csv")
    operating = pd.read_csv("results/v3/operating_models.csv")
    row = operating[
        (operating["species"] == species) & (operating["baseline_end"] == 1985)
    ].iloc[0]
    series, params, core = fit_baseline(regional, species)
    assert params.b == pytest.approx(row["b"], abs=1e-4)
    assert params.q == pytest.approx(row["q"], abs=1e-4)
    assert " ".join(series.regions[i] for i in core) == row["core_regions"]
