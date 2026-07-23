# Population Model
[![Python Applications test with GitHub Actions](https://github.com/benmacp123/population-model/actions/workflows/main.yml/badge.svg)](https://github.com/benmacp123/population-model/actions/workflows/main.yml)

A simple two-stage (juvenile/adult) population model with a backtest step: project a species forward from a historical starting year, then compare the projection against the observed counts that followed.

## Approach

The model is initialised from one year of observed counts, projected forward year-by-year, and then validated against any later observed counts in the dataset. The gap between projected and observed is the point — a forecast on its own can't be checked, a backtest can.

## Current work

`growth_model.py` is a small executable model that:

- loads two separate placeholder datasets — species traits and year-by-year counts
- computes carrying capacity, either from an explicit per-species value or from a total-biomass proxy
- projects adult and juvenile populations forward from a chosen start year
- applies birth, mortality and maturation dynamics with density dependence
- compares the projection against observed counts and reports absolute and percentage error
- exposes all of this through a `click` CLI

This is a prototype framework rather than a finished ecological model. The simplifications below are deliberate and documented.

## Data structure

Two inputs are held separately, matching how the data actually gets sourced:

| Dataset | Grain | Source in practice |
|---|---|---|
| Species traits | one row per species | published literature — gestation, lifespan, maturation age, sex ratio, average weight |
| Species counts | one row per species per year | aerial census / survey / monitoring data |

Traits are stable constants; counts vary year to year. Keeping them apart means each carries its own provenance, and it avoids duplicating the biology across every count row.

Two trait fields are optional overrides, left as `None` by default:

- `juv_death_rate` — falls back to the adult mortality rate
- `carrying_capacity` — falls back to the biomass proxy below

Where a real value is available from literature, set it in the traits table; no code change needed.

## Model formulation

- **Carrying capacity** — explicit value where supplied, otherwise a biomass proxy from the start year:
  ```text
  K = (Σ_i N_i * w_i) / w_target
  ```
  where `N_i` is total population of species `i`, `w_i` its average weight, and `w_target` the target species' average weight.

- **Births per year** (females only):
  ```text
  B = N_adult * (1 - male_ratio) * r_birth
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
  r_mat = (1 / m_male) * male_ratio + (1 / m_female) * (1 - male_ratio)
  f(K, N) = max(0, 1 - (N_adult + N_juvenile) / K)
  ```

- **Population update**, floored at zero:
  ```text
  N_adult'    = max(0, N_adult + M - D_adult)
  N_juvenile' = max(0, N_juvenile + B - D_juvenile - M)
  ```

- **Validation**, for years present in both projection and observations:
  ```text
  error     = N_projected - N_observed
  pct_error = 100 * error / N_observed
  ```

## Known simplifications

Recorded openly rather than buried — these are the model's honest limits and the first targets for improvement.

1. **Juvenile mortality defaults to the adult rate.** In reality juvenile mortality is usually much higher and often drives the dynamics. With the fallback in place the model tends to overshoot.
2. **Carrying capacity is partly circular.** The biomass proxy infers capacity from the population already present rather than from habitat or resource limits, and it shifts depending on which species are included in the dataset.
3. **Birth rate assumes continuous breeding.** `12 / gestation_months` implies back-to-back reproduction with no inter-birth interval and no seasonality.
4. **Mortality is a flat annual hazard.** `1 / life_span_years` is not an age-structured survival curve.
5. **Density dependence applies only to maturation**, not to births or mortality.

## Planned next steps

1. **Connect to real datasets**
   - replace placeholder traits and counts with real survey and literature values, each with a source citation
   - load from versioned CSVs rather than in-code definitions
   - normalise incoming data to the common schema

2. **Calibrate against observed counts**
   - use the validation output to fit juvenile mortality and carrying capacity rather than assuming them
   - report fit quality across species

3. **Add richer demographic structure**
   - age classes beyond adult/juvenile
   - sex-specific populations and fertility
   - explicit stage transitions rather than a single maturation term

4. **Incorporate predation and density dependence**
   - predator-driven mortality or harvest pressure
   - resource limits affecting birth and survival rates
   - Allee effects or density-dependent reproduction

5. **Move toward spatial / GIS-aware modelling**
   - habitat quality and resource availability as spatial layers
   - carrying capacity mapped across cells or landscape units
   - movement, dispersal and local density dynamics

## Getting started

1. Create and activate a Python virtual environment.
2. Install requirements:
   ```bash
   pip install -r requirements.txt
   ```
3. Run the model:
   ```bash
   ./growth_model.py --species_name Zebra --start_year 2018 --forecast_horizon 6
   ```

`--start_year` must exist in the counts dataset for the chosen species, since the projection is initialised from observed values.

Example output:

```text
Projection for Zebra from 2018:
  2018: adults=700 juveniles=300 total=1000
  ...

Validation against observed counts:
  2020: projected=1586 observed=1180 error=+406 (+34.4%)
  2022: projected=2112 observed=1330 error=+782 (+58.8%)
```

## Notes

- This repository is a lightweight prototype and should be extended before use in production analysis.
- Placeholder figures are illustrative only and are not sourced from real survey data.
