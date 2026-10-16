# Population Model
[![Python Applications test with GitHub Actions](https://github.com/benmacp123/population-model/actions/workflows/main.yml/badge.svg)](https://github.com/benmacp123/population-model/actions/workflows/main.yml)

A simple two-stage (juvenile/adult) population model with a backtest step: project a species forward from a historical starting year, then compare the projection against the observed counts that followed.

## Approach

The model is initialised from one year of observed counts, projected forward year-by-year, and then validated against any later observed counts in the dataset. The gap between projected and observed is the point - a forecast on its own can't be checked, a backtest can.

## Current work

The prototype runs end to end from committed data. It:

- keeps species traits and year-by-year counts as separate datasets
- derives `data/raw/species_traits.csv` from the committed PanTHERIA WR05 file
- normalizes the SANParks 1965-1997 fixed-wing total-count series into
  `data/raw/species_counts.csv`, with a source citation on every row and
  missing totals left missing rather than zero
- computes carrying capacity, either from an explicit per-species value or from a total-biomass proxy
- projects adult and juvenile populations forward from a chosen start year
- applies birth, mortality and maturation dynamics with density dependence
- compares the projection against observed counts and reports absolute and percentage error
- reports backtest accuracy (`make evaluate-model`) and a diagnostic
  assumption grid (`make tune-assumptions`)
- exposes all of this through a `click` CLI

The active species are wildebeest, waterbuck and kudu. Zebra is excluded
until its missing reproductive traits are sourced. The census series stops at
1997 because Kruger moved from total counts to distance sampling in 1998.

**Status: exploratory.** With interbirth-interval births, juvenile-mortality
scenarios and a calibration-period capacity, the 1980-1997 holdout error is
about 16-38% of observed totals, roughly the same as assuming the population
never changed. The model does not reproduce the observed 1980s-1990s
declines (see [Known issues](#known-issues)). Forecasts should not be read as
ecological predictions. Progress against [`PLAN.md`](PLAN.md) is recorded
there.

## Data structure

Two inputs are held separately, matching how the data actually gets sourced:

| Dataset | Grain | Source in practice |
|---|---|---|
| Species traits | one row per species | PanTHERIA and other published literature |
| Species counts | one row per species per year | Kruger aerial census / survey data |

Traits are stable constants; counts vary year to year. Keeping them apart means each carries its own provenance, and it avoids duplicating the biology across every count row.

Two trait fields are optional overrides, left as `None` by default:

- `juv_death_rate` - falls back to the adult mortality rate
- `carrying_capacity` - falls back to the biomass proxy below

Where a real value is available from literature, set it in the traits table; no code change needed.

The traits table can be regenerated with `make build-traits`. Its provenance
and missing-value conventions are recorded in `data/README.md`.

## Model formulation

- **Carrying capacity** - explicit value where supplied, otherwise a biomass proxy from the start year:
  ```text
  K = (Σ_i N_i * w_i) / w_target
  ```
  where `N_i` is total population of species `i`, `w_i` its average weight, and `w_target` the target species' average weight.

- **Births per year** (females only):
  ```text
  B = N_adult * (1 -  male_ratio) * r_birth
  r_birth = 365.25 / interbirth_interval_days   (fallback: 12 / gestation_period_months)
  ```

- **Mortality**:
  ```text
  d_adult = 1 / life_span_years
  d_juvenile = juv_death_rate  (falls back to d_adult)
  D_adult = N_adult * d_adult
  D_juvenile = N_juvenile * d_juvenile
  ```

- **Maturation from juveniles to adults**, sex-weighted and density-limited:
  ```text
  M = N_juvenile * r_mat * f(K, N)
  r_mat = (1 / m_male) * male_ratio + (1 / m_female) * (1 -  male_ratio)
  f(K, N) = max(0, 1 -  (N_adult + N_juvenile) / K)
  ```

- **Population update**, floored at zero:
  ```text
  N_adult'    = max(0, N_adult + M -  D_adult)
  N_juvenile' = max(0, N_juvenile + B -  D_juvenile -  M)
  ```

- **Validation**, for years present in both projection and observations:
  ```text
  error     = N_projected -  N_observed
  pct_error = 100 * error / N_observed
  ```

## Known simplifications

Recorded openly rather than buried - these are the model's honest limits and the first targets for improvement.

1. **Juvenile mortality defaults to the adult rate.** In reality juvenile mortality is usually much higher and often drives the dynamics. With the fallback in place the model tends to overshoot.
2. **Carrying capacity is partly circular.** The biomass proxy infers capacity from the population already present rather than from habitat or resource limits, and it shifts depending on which species are included in the dataset.
3. **Birth rate has no seasonality or age structure.** One young per interbirth interval for every adult female. Species without an interval (tsessebe, warthog) fall back to `12 / gestation_months`, i.e. back-to-back breeding.
4. **Mortality is a flat annual hazard.** `1 / life_span_years` is not an age-structured survival curve.
5. **Density dependence applies only to maturation**, not to births or mortality.
6. **The adult/juvenile split is assumed, not observed.** Census totals have no stage breakdown, so every total is split 50/50 as a provisional placeholder until a literature-sourced stable stage proportion is selected. Override with `--juvenile_fraction`.
7. **Sex-specific maturity is copied from one value.** PanTHERIA gives a single sexual maturity age, which is used for both males and females.
8. **The male ratio is a 50/50 placeholder**, not a sourced value. Override with `--male_ratio`.

## Known issues

Holdout results from `make evaluate-model` (start 1965, calibration
1965-1979, holdout 1980-1997, provisional 50/50 stage split and sex ratio,
calibration-peak capacity). Juvenile mortality scenarios are 0.2/0.35/0.5
per year, a sensitivity range rather than sourced values.

| Species | Holdout error, % of observed mean | Naive constant baseline | Skill vs baseline |
|---|---|---|---|
| Wildebeest | 16-24% | 13% | -0.77 to -0.22 |
| Waterbuck | 29-35% | 38% | +0.08 to +0.24 |
| Kudu | 33-38% | 34% | -0.13 to +0.01 |

Skill above zero means the model beats holding the 1965 count constant.
Before this pass, with gestation-based births and community-biomass
capacity, error was 1,800-4,200% of the observed mean.

Open issues:

1. **No useful skill yet.** The model roughly matches the naive baseline: it
   is slightly better for waterbuck and worse for wildebeest.
2. **Observed declines are missed.** All three species fell substantially
   through the 1980s-1990s (waterbuck from ~4,000 to ~960). The model
   settles near capacity and is biased high by 5-24%. Rainfall and drought,
   predation and management are not represented.
3. **The capacity that works is circular.** The calibration-peak capacity is
   each species' own highest 1965-1979 count. It keeps the holdout honest but
   is a diagnostic, not a habitat estimate.
4. **The community-biomass proxy is still inflated.** It remains available
   (`--capacity-mode community-biomass`) and is still used by the
   single-species CLI. It treats all counted herbivore biomass as available
   to one species, giving capacities 8-30 times the observed peak, and it
   silently drops species with no start-year total (buffalo in 1965).
5. **Juvenile mortality and the stage split are unsourced.** The scenario
   range and the 50/50 split are placeholders pending literature values.
6. **Tooling versions are unpinned.** `requirements.txt` does not pin
   versions, so CI may install a different `black` or `pylint` from the one
   used locally, and the formatting check can fail on unchanged code.

## Planned next steps

The detailed, checkpointed plan is in [`PLAN.md`](PLAN.md). A lean pass
has covered the fair baseline (holdout scoring, naive baseline, trend
metrics), interbirth-interval births, juvenile-mortality scenarios and a
diagnostic capacity. Remaining, if this phase is revisited:

1. Observed-versus-projected plots.
2. Literature-sourced stage split and juvenile mortality.
3. A habitat- or forage-based carrying capacity.
4. Bounded calibration with low/central/high scenarios.
5. An environmental driver (for example rainfall) to test whether it
   explains the 1980s-1990s declines.

Later portfolio phases, deliberately out of scope here: richer age/sex
structure, predation, and spatial or habitat-based modelling.

## Getting started

1. Create and activate a Python virtual environment.
2. Install requirements:
   ```bash
   make install
   ```
3. Run the model:
   ```bash
   ./src/growth_model.py --species_name Wildebeest --start_year 1965 --forecast_horizon 6
   ```
4. Evaluate the 1980-1997 holdout across juvenile-mortality scenarios:
   ```bash
   make evaluate-model
   ```
5. Run the tests, lint and formatting check:
   ```bash
   make test lint format-check
   ```

`--start_year` must exist in the counts dataset for the chosen species, since the projection is initialised from observed values.

Example output (single-species CLI with the community-biomass capacity
proxy, which still overshoots; see Known issues 4):

```text
$ ./src/growth_model.py --species_name Wildebeest --start_year 1965 --forecast_horizon 6 --juv_death_rate 0.35

Projection for Wildebeest from 1965:
  1965: adults=6518 juveniles=6518 total=13035
  1966: adults=8747 juveniles=4965 total=13712
  ...
  1971: adults=18454 juveniles=10810 total=29264

Validation against observed counts:
  1965: projected=13035 observed=13035 error=+0 (+0.0%)
  1969: projected=21687 observed=16072 error=+5615 (+34.9%)
  1970: projected=25223 observed=13950 error=+11273 (+80.8%)
  1971: projected=29264 observed=12557 error=+16707 (+133.0%)
```

## Notes

- This repository is a research prototype and should be extended before use in production analysis.
- Data provenance for every committed file is recorded in [`data/README.md`](data/README.md).
