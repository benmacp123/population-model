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

**Status: exploratory.** The current backtest overshoots observed totals by
an order of magnitude for all three species (see [Known issues](#known-issues)).
Forecasts should not be read as ecological predictions yet. The improvement
path is set out in [`PLAN.md`](PLAN.md).

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
  r_birth = 12 / gestation_period_months
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
3. **Birth rate assumes continuous breeding.** `12 / gestation_months` implies back-to-back reproduction with no inter-birth interval and no seasonality.
4. **Mortality is a flat annual hazard.** `1 / life_span_years` is not an age-structured survival curve.
5. **Density dependence applies only to maturation**, not to births or mortality.
6. **The adult/juvenile split is assumed, not observed.** Census totals have no stage breakdown, so every total is split 50/50 as a provisional placeholder until a literature-sourced stable stage proportion is selected. Override with `--juvenile_fraction`.
7. **Sex-specific maturity is copied from one value.** PanTHERIA gives a single sexual maturity age, which is used for both males and females.
8. **The male ratio is a 50/50 placeholder**, not a sourced value. Override with `--male_ratio`.

## Known issues

Observed in the current 1965-1997 backtest (default 50/50 assumptions,
start year 1965). These are the targets for the next round of work in
[`PLAN.md`](PLAN.md).

| Species | Normalized MAE | Capacity proxy (K) | Highest observed total | Implied births per female per year |
|---|---|---|---|---|
| Wildebeest | ~1,800% | ~124,700 | 16,072 | 1.45 |
| Waterbuck | ~4,200% | ~121,200 | 5,042 | 1.35 |
| Kudu | ~3,500% | ~120,200 | 10,760 | 1.65 |

1. **Runaway growth.** All three species grow roughly exponentially and
   overshoot observed totals by an order of magnitude or more. Simplifications
   1, 3 and 5 compound: births are too high, juveniles die too slowly, and
   density dependence never restrains births.
2. **The birth rate is implausibly high.** `12 / gestation_months` gives
   1.3-1.7 young per female per year. These species typically produce about
   one per year, and `interbirth_interval_days` (365 days for all three) is
   already in the traits table but unused.
3. **Carrying capacity barely binds.** The biomass proxy treats the entire
   counted herbivore community's biomass as available to each single species,
   so K is roughly 8-30 times the highest observed total and the density
   factor stays near 1 until the projection is far too large.
4. **The capacity proxy silently drops species with no start-year count.**
   Missing totals are skipped in the biomass sum, so K depends on which
   species happened to be reported in the start year (buffalo, for example,
   has no 1965 total). Changing the start year can redefine K for reasons
   unrelated to the target species.
5. **The start year is scored as a prediction.** The initialization row has
   zero error by construction but is included in the accuracy metrics, which
   flatters them slightly.
6. **Accuracy metrics measure magnitude only.** There is no trend-direction
   score and no comparison against a naive constant-population baseline, so
   the report cannot yet say whether the model captures increase or decline.
7. **Tooling versions are unpinned.** `requirements.txt` does not pin
   versions, so CI may install a different `black` or `pylint` from the one
   used locally, and the formatting check can fail on unchanged code.

## Planned next steps

The detailed, checkpointed plan is in [`PLAN.md`](PLAN.md). In summary:

1. **Establish a fair baseline** - separate initialization from evaluation,
   add trend metrics and a constant-population baseline, plot observed versus
   projected totals.
2. **Fix the largest biological pressure points** - interbirth-interval
   births, explicit juvenile mortality, sensitivity-tested stage split.
3. **Improve carrying capacity** - make its source explicit and stop it
   depending silently on the species set.
4. **Calibrate for trend, not exact counts** - bounded fits with a holdout
   period and low/central/high scenarios.
5. **Make the result presentable** - reproducible evaluation command, results
   table and plots.

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
4. Evaluate the full 1965-1997 backtest:
   ```bash
   make evaluate-model
   ```
5. Run the tests, lint and formatting check:
   ```bash
   make test lint format-check
   ```

`--start_year` must exist in the counts dataset for the chosen species, since the projection is initialised from observed values.

Example output (current model, illustrating the overshoot described above):

```text
Projection for Wildebeest from 1965:
  1965: adults=6518 juveniles=6518 total=13035
  1966: adults=8747 juveniles=8415 total=17162
  ...
  1971: adults=30124 juveniles=36709 total=66834

Validation against observed counts:
  1965: projected=13035 observed=13035 error=+0 (+0.0%)
  1969: projected=39456 observed=16072 error=+23384 (+145.5%)
  1970: projected=51570 observed=13950 error=+37620 (+269.7%)
  1971: projected=66834 observed=12557 error=+54277 (+432.2%)
```

## Notes

- This repository is a research prototype and should be extended before use in production analysis.
- Data provenance for every committed file is recorded in [`data/README.md`](data/README.md).
