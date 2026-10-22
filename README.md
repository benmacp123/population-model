# Population Model
[![Python Applications test with GitHub Actions](https://github.com/benmacp123/population-model/actions/workflows/main.yml/badge.svg)](https://github.com/benmacp123/population-model/actions/workflows/main.yml)

A simple two-stage (juvenile/adult) population model for Kruger National Park herbivores, calibrated on 1965-1979 census totals and tested against the 1980-1995 counts it never saw (1996-1997 have no reported totals).

## Approach

The model is fitted on an early period of observed counts, restarted from the last calibration observation, projected forward year-by-year, and scored against the later observed counts. The gap between projected and observed is the point - a forecast on its own can't be checked, a backtest can. Every score is compared with two naive baselines, one carrying the last known count forward and one using the mean of the calibration-period counts, so "better than doing nothing" is measured rather than assumed.

## Current work

The model runs end to end from committed data. It:

- keeps species traits and year-by-year counts as separate datasets
- derives `data/raw/species_traits.csv` from the committed PanTHERIA WR05 file
- normalizes the SANParks 1965-1997 fixed-wing total-count series into
  `data/raw/species_counts.csv`, with a source citation on every row and
  missing totals left missing rather than zero
- applies birth, mortality and maturation dynamics with density dependence
- calibrates three bounded parameters per species on 1965-1979 and reports
  low/central/high scenarios from the near-best parameter sets
- scores the holdout on level and trend against persistence and
  calibration-mean baselines, and writes a results table and figure to `results/`
- exposes a single-species projection through a `click` CLI

The active species are wildebeest, waterbuck and kudu. Zebra is excluded
until its missing reproductive traits are sourced. The census series stops at
1997 because Kruger moved from total counts to distance sampling in 1998.

**Status: complete, result negative.** The model reproduces the
calibration period to within 11-13% but has no forecasting skill on the
holdout: it never meaningfully beats carrying the 1979 count forward, is
well behind simply using the 1965-1979 mean, and misses the 1980s rise and
1990s collapse. That finding holds across the stage split, sex ratio, grid
resolution, search bounds and restart rule; the exact per-species skill
values do not. See
[Results](#results). This is the model's honest finding, not a bug; the most
likely missing ingredient is rainfall.

## Results

`make results` reproduces the table and figure below
(`results/backtest_summary.csv`, `results/backtest.png`).

![Observed census totals against calibrated and holdout projections](results/backtest.png)

| Species | Calibration error | Holdout error | Persistence baseline | Skill | Step direction agreement | Observed inside scenario range |
|---|---|---|---|---|---|---|
| Wildebeest | 11.0% | 42.2% | 39.5% | -0.07 | 53% | 0% |
| Waterbuck | 10.6% | 33.9% | 33.5% | -0.01 | 67% | 6% |
| Kudu | 12.8% | 45.9% | 36.4% | -0.26 | 40% | 0% |

A second naive baseline, the mean of the reported 1965-1979 totals, is
harder to beat than persistence and uses no holdout data:

| Species | Calibration-mean baseline error | Model skill against it |
|---|---|---|
| Wildebeest | 15.4% | -1.73 |
| Waterbuck | 30.0% | -0.13 |
| Kudu | 34.0% | -0.35 |

Errors are mean absolute error as a percentage of the observed mean.
Skill is `1 - model error / baseline error`; above zero beats the baseline.
Scores cover the 16 reported holdout years, 1980-1995. The 1979 anchor
shared by the forecast and the persistence baseline is partly an estimate
(see [Data caveats](#data-caveats)).
Step direction agreement is the share of consecutive reported years where
the projection moves the same way as the census. The persistence baseline
never moves, so the model trivially beats it on direction; level skill is
the meaningful test.

Fitted parameters (central scenario). These are grid- and stage-split-specific
fitted values, not vital-rate estimates: halving the grid step moves them,
and changing the provisional juvenile fraction moves fitted juvenile
mortality between 0.25 and 0.60.

| Species | Juvenile mortality | Birth multiplier | Capacity (x calibration peak) | On search bound |
|---|---|---|---|---|
| Wildebeest | 0.40 | 0.9 | 0.5 | capacity (grid edge only) |
| Waterbuck | 0.40 | 0.6 | 0.9 | none |
| Kudu | 0.55 | 0.5 | 1.5 | birth, capacity |

What the results show:

1. **The model fits the past but cannot forecast it.** Fitted on 1965-1979,
   it settles towards an equilibrium below its capacity parameter. For
   wildebeest that equilibrium sits close to the 1979 anchor, so the
   projection is nearly flat and scores close to persistence by
   construction. The census instead rose sharply in the early
   1980s and then fell 38-81% from peak by 1995.
2. **The dynamics changed between periods.** As a diagnostic, fitting the
   same three parameters directly on 1980-1997 (`make holdout-ceiling`)
   beats persistence comfortably (skill 0.86 wildebeest, 0.31 waterbuck,
   0.42 kudu): the stage structure can overshoot capacity and fall back,
   following the 1980s rise and a gentle decline. But those parameters miss
   1965-1979 by 33-73%, and none reproduce the steep 1991-95 crash. No
   single parameter set describes both periods. A changing external driver
   is one explanation; three compensating parameters fitted to 8 noisy points
   would also transfer poorly, so this diagnostic does not separate the two.
   The diagnostic restarts from a 50/50 split in 1979 and has not been
   rerun with the carried stage share used by the forecast. Published Kruger
   work links these declines to rainfall and the 1991-92 drought (for
   example Ogutu & Owen-Smith 2003, *Ecology Letters* 6: 412-419), and the
   1970s were unusually wet. These diagnostic numbers are not forecast
   scores.
3. **The scenario range is too narrow.** Near-best parameter sets agree with
   each other, so the band almost never contains the observed counts. The
   uncertainty that matters here is structural, not parametric.
4. **Kudu fits on a search bound.** It needs the lowest birth rate and
   highest capacity allowed. Widening the capacity range to 0.3-3.0 improves
   its calibration error (12.8% to 12.6%) but worsens holdout skill (-0.26 to
   -0.38), a small illustration of overfitting. Wildebeest's lowest-capacity
   fit is a grid artefact: over the same wider range its optimum stays at an
   interior 0.5.
5. **The restart carries the model's own stage split.** The 1979 restart
   uses the observed total but splits it using each parameter set's own
   1979 juvenile share (0.56 wildebeest, 0.39 waterbuck, 0.24 kudu for the
   central fit), not the provisional 50/50 seed. An earlier version re-seeded
   50/50, which caused an artificial first-year kudu drop and moved skill by
   up to 0.15 (published then as -0.04, -0.16 and -0.18). The switch was made
   after those holdout results had been seen; it is a consistency fix, not a
   tuning choice, and it does not change the finding.

## Data caveats

- **Some totals include regional infill.** In 1970, 1978, 1979, 1994 and
  1995 the source `TOTAL` exceeds the sum of the regional counts by a round
  number, so unsurveyed regions were estimated. The 1979 anchor includes
  1200 of 7412 wildebeest (16%), 2500 of 8974 kudu (28%) and 450 of 2605
  waterbuck (17%) for the uncounted South region. 1965 has no regional
  breakdown for kudu or waterbuck.
- **The series spans more than one method.** The DataONE record describes
  1965-1976 as aerial totals corrected by ground counts and 1977-1997 as
  Ecological Aerial Survey totals; the KNB copy states park coverage was
  complete only in 1980-1993. The 1998 switch to distance sampling is a
  further break, which is why the series stops at 1997.
- **The calibration window was managed.** Wildebeest and zebra were culled
  in the late 1960s and early 1970s, and lions in 1974-75, so 1965-1979 is
  not unmanaged dynamics. These dates come from secondary sources.

## Data structure

Two inputs are held separately, matching how the data actually gets sourced:

| Dataset | Grain | Source in practice |
|---|---|---|
| Species traits | one row per species | PanTHERIA and other published literature |
| Species counts | one row per species per year | Kruger aerial census / survey data |

Traits are stable constants; counts vary year to year. Keeping them apart means each carries its own provenance, and it avoids duplicating the biology across every count row.

Three trait fields are optional overrides, left empty by default:

- `juv_death_rate` - falls back to the adult mortality rate
- `carrying_capacity` - falls back to the biomass proxy below
- `birth_multiplier` - falls back to 1 (one young per interbirth interval)

Where a real value is available from literature, set it in the traits table; no code change needed.

The traits table can be regenerated with `make build-traits`. Its provenance
and missing-value conventions are recorded in `data/README.md`.

## Model formulation

- **Carrying capacity** - explicit value where supplied, otherwise a biomass proxy from the start year:
  ```text
  K = (Σ_i N_i * w_i) / w_target
  ```
  where `N_i` is total population of species `i`, `w_i` its average weight, and `w_target` the target species' average weight.
  The calibrated evaluation instead uses `K = capacity_multiplier * peak 1965-1979 total`.

- **Births per year** (females only):
  ```text
  B = N_adult * (1 -  male_ratio) * r_birth * birth_multiplier
  r_birth = 365.25 / interbirth_interval_days   (fallback: 12 / gestation_period_months)
  ```

- **Mortality**:
  ```text
  d_adult = 1 / life_span_years
  d_juvenile = juv_death_rate  (falls back to d_adult)
  D_adult = N_adult * d_adult
  D_juvenile = N_juvenile * d_juvenile
  ```

- **Maturation from juveniles to adults**, sex-weighted, density-limited, and capped at surviving juveniles:
  ```text
  M = min(N_juvenile * r_mat * f(K, N), N_juvenile * (1 - d_juvenile))
  r_mat = (1 / m_male) * male_ratio + (1 / m_female) * (1 -  male_ratio)
  f(K, N) = max(0, 1 -  (N_adult + N_juvenile) / K)
  ```

- **Population update**, floored at zero:
  ```text
  N_adult'    = max(0, N_adult + M -  D_adult)
  N_juvenile' = max(0, N_juvenile + B -  D_juvenile -  M)
  ```

- **Calibration**: grid search over juvenile mortality 0.10-0.60, birth
  multiplier 0.5-1.0 and capacity multiplier 0.5-1.5, minimising mean
  absolute error on the 8 reported 1966-1979 totals (1969-1973, 1975, 1978
  and 1979). Parameter sets within 10% of
  the best error form the low/high scenario envelope. These bounds are a
  defensible search range, not sourced values.

## Known simplifications

Recorded openly rather than buried - these are the model's honest limits.

1. **Juvenile mortality is fitted, not sourced.** The trait-table fallback is the adult rate; the evaluation fits it within 0.10-0.60.
2. **Carrying capacity is circular.** Both the biomass proxy and the calibration-peak capacity are inferred from counted animals rather than habitat or forage.
3. **Birth rate has no seasonality or age structure.** One young per interbirth interval, scaled by a fitted multiplier. Species without an interval (tsessebe, warthog) fall back to `12 / gestation_months`; warthog also lacks a maturity age, so it cannot be run at all.
4. **Mortality is a flat annual hazard.** `1 / life_span_years` is not an age-structured survival curve, and PanTHERIA's lifespan is the maximum recorded longevity, so the 4-5% a year it gives probably understates field adult mortality. Newborns also face no mortality in the step they are born, and maturation is a constant yearly fraction rather than a delay.
5. **Density dependence applies only to maturation**, not to births or mortality.
6. **The adult/juvenile split is assumed, not observed.** Census totals have no stage breakdown, so every total is split 50/50 as a provisional placeholder until a literature-sourced stable stage proportion is selected. Override with `--juvenile_fraction`.
7. **Sex-specific maturity is copied from one value.** PanTHERIA gives a single sexual maturity age, which is used for both males and females.
8. **The male ratio is a 50/50 placeholder**, not a sourced value. Override with `--male_ratio`.
9. **No environmental drivers.** Rainfall, drought, predation, disease and management removals are absent, so the model cannot produce the observed boom-and-bust.

## Known issues

1. **No holdout skill** for any species (see [Results](#results)).
2. **The single-species CLI still uses the community-biomass proxy.** It treats
   all counted herbivore biomass as available to one species, giving
   capacities 8-30 times the observed peak, and silently drops species with no
   start-year total (buffalo in 1965). Use `make evaluate-model` for the
   calibrated model.
3. **Kudu fits on search bounds**, a sign of structural misfit (wildebeest's apparent bound is a grid artefact).

## What I would do next

The project is closed, but a review of the method identified what a second
iteration should change:

1. **Model observation error explicitly.** Single-year swings of +82% (kudu)
   and +55% (waterbuck) exceed what these species can biologically achieve,
   so much of the year-to-year variation is survey noise. A state-space
   model (density-dependent growth with separate process and observation
   error) separates the two - the same problem as state estimation with a
   Kalman filter in power-system monitoring.
2. **Use fewer, identifiable parameters.** PanTHERIA traits are
   species-level life history, not population vital rates, and three
   parameters fitted to 8 calibration points are poorly identified. A
   Gompertz or Ricker model with traits used only as a sanity check on the
   fitted growth rate is more defensible.
3. **Test forecasts more fairly.** One-step-ahead and rolling-origin
   forecasts give far more evaluation points than a single 15-year
   projection, plus a check that uncertainty bands cover the observations at
   their stated rate.
4. **Add rainfall as a covariate.** SANParks publishes Kruger rainfall
   (1911-2003) on the same DataONE portal as the census; published Kruger
   studies link these declines to rainfall.

## Getting started

1. Create and activate a Python 3.12+ virtual environment.
2. Install requirements:
   ```bash
   make install
   ```
3. Reproduce the calibration, holdout scores and figure:
   ```bash
   make results
   ```
4. Run a single-species projection:
   ```bash
   ./src/growth_model.py --species_name Wildebeest --start_year 1965 --forecast_horizon 6 --juv_death_rate 0.35
   ```
   `--start_year` must exist in the counts dataset for the chosen species, since the projection is initialised from observed values.
5. Run the tests, lint and formatting check:
   ```bash
   make test lint format-check
   ```

## Notes

- This repository is a research prototype, not a management tool.
- Data provenance for every committed file is recorded in [`data/README.md`](data/README.md).
- Dependencies are pinned in `requirements.txt` so CI and local runs use the same tools.

## Licence

Code, documentation and results are released under the MIT Licence; see
[`LICENSE`](LICENSE). Third-party data in `data/raw/` stays under its
source's terms, recorded in [`data/README.md`](data/README.md).
