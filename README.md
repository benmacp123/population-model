# Population Model
[![Python Applications test with GitHub Actions](https://github.com/benmacp123/population-model/actions/workflows/main.yml/badge.svg)](https://github.com/benmacp123/population-model/actions/workflows/main.yml)

A simple two-stage (juvenile/adult) population model for Kruger National Park herbivores, calibrated on 1965-1979 census totals and tested against the 1980-1997 counts it never saw.

## Approach

The model is fitted on an early period of observed counts, restarted from the last calibration observation, projected forward year-by-year, and scored against the later observed counts. The gap between projected and observed is the point - a forecast on its own can't be checked, a backtest can. Every score is compared with a naive baseline that simply carries the last known count forward, so "better than doing nothing" is measured rather than assumed.

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
- scores the 1980-1997 holdout on level and trend against a persistence
  baseline, and writes a results table and figure to `results/`
- exposes a single-species projection through a `click` CLI

The active species are wildebeest, waterbuck and kudu. Zebra is excluded
until its missing reproductive traits are sourced. The census series stops at
1997 because Kruger moved from total counts to distance sampling in 1998.

**Status: complete, result negative.** The model reproduces the
calibration period to within 11-13% but has no forecasting skill on the
holdout: it is slightly worse than carrying the 1979 count forward for all
three species, and misses the 1980s rise and 1990s collapse. See
[Results](#results). This is the model's honest finding, not a bug; the most
likely missing ingredient is rainfall.

## Results

`make results` reproduces the table and figure below
(`results/backtest_summary.csv`, `results/backtest.png`).

![Observed census totals against calibrated and holdout projections](results/backtest.png)

| Species | Calibration error | Holdout error | Persistence baseline | Skill | Step direction agreement | Observed inside scenario range |
|---|---|---|---|---|---|---|
| Wildebeest | 11.0% | 41.1% | 39.5% | -0.04 | 27% | 0% |
| Waterbuck | 10.6% | 38.8% | 33.5% | -0.16 | 27% | 6% |
| Kudu | 12.8% | 43.0% | 36.4% | -0.18 | 33% | 0% |

Errors are mean absolute error as a percentage of the observed mean.
Skill is `1 - model error / baseline error`; above zero beats persistence.
Step direction agreement is the share of consecutive reported years where
the projection moves the same way as the census. The persistence baseline
never moves, so the model trivially beats it on direction; level skill is
the meaningful test.

Fitted parameters (central scenario):

| Species | Juvenile mortality | Birth multiplier | Capacity (x calibration peak) | On search bound |
|---|---|---|---|---|
| Wildebeest | 0.40 | 0.9 | 0.5 | capacity |
| Waterbuck | 0.40 | 0.6 | 0.9 | none |
| Kudu | 0.55 | 0.5 | 1.5 | birth, capacity |

What the results show:

1. **The model fits the past but cannot forecast it.** A single-equilibrium
   model settles towards capacity. The census instead rose sharply in the
   early 1980s and then fell 38-81% from peak by 1995. Published Kruger work
   links these declines to rainfall and the 1991-92 drought (for example
   Ogutu & Owen-Smith 2003, *Ecology Letters* 6: 412-419).
2. **The scenario range is too narrow.** Near-best parameter sets agree with
   each other, so the band almost never contains the observed counts. The
   uncertainty that matters here is structural, not parametric.
3. **Two species fit on a search bound.** Kudu needs the lowest birth rate
   and highest capacity allowed, and wildebeest the lowest capacity. Widening
   the bounds would fit calibration better without fixing the holdout.
4. **The restart split shows.** Restarting in 1979 re-seeds a 50/50
   adult/juvenile split, which produces the sharp first-year kudu drop.

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
  absolute error on reported 1966-1979 totals. Parameter sets within 10% of
  the best error form the low/high scenario envelope. These bounds are a
  defensible search range, not sourced values.

## Known simplifications

Recorded openly rather than buried - these are the model's honest limits.

1. **Juvenile mortality is fitted, not sourced.** The trait-table fallback is the adult rate; the evaluation fits it within 0.10-0.60.
2. **Carrying capacity is circular.** Both the biomass proxy and the calibration-peak capacity are inferred from counted animals rather than habitat or forage.
3. **Birth rate has no seasonality or age structure.** One young per interbirth interval, scaled by a fitted multiplier. Species without an interval (tsessebe, warthog) fall back to `12 / gestation_months`.
4. **Mortality is a flat annual hazard.** `1 / life_span_years` is not an age-structured survival curve.
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
3. **Kudu and wildebeest fit on search bounds**, a sign of structural misfit.

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
   parameters fitted to ~9 calibration points are poorly identified. A
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
