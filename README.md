# Population Model
[![Python Applications test with GitHub Actions](https://github.com/benmacp123/simple-population-model/actions/workflows/main.yml/badge.svg)](https://github.com/benmacp123/simple-population-model/actions/workflows/main.yml)

Simple population modelling based on maturation rates of individual species, sex-specific dynamics, and carrying capacity.

## Current work

This project currently includes a small executable model in `growth_model.py` that:

- loads placeholder species data via `pandas`
- computes carrying capacity from total biomass and average weight
- projects adult and juvenile populations year-by-year
- applies birth, death, and maturation dynamics at the population level
- supports a CLI interface using `click`

The current implementation is designed to be a prototype framework rather than a final ecological model.

## What it currently does

The current model uses a simple demographic projection with the following core equations:

- Carrying capacity proxy:
  $$K = \frac{\sum_{i} N_i \cdot w_i}{w_{\text{target}}}$$
  where $N_i$ is population size, $w_i$ is average weight, and $w_{\text{target}}$ is the target species average weight.

- Births per year:
  $$B = N_{\text{adult}} \cdot r_{\text{birth}}$$
  with $r_{\text{birth}} = \frac{12}{\text{gestation\_period\_months}}$.

- Mortality rates:
  $$d = \frac{1}{\text{life\_span\_years}}$$
  $$D_{\text{adult}} = N_{\text{adult}} \cdot d$$
  $$D_{\text{juvenile}} = N_{\text{juvenile}} \cdot d$$

- Maturation from juveniles to adults:
  $$M = N_{\text{juvenile}} \cdot r_{\text{mat}} \cdot f(K, N)$$
  where $r_{\text{mat}}$ is the sex-weighted maturation rate and $f(K, N)$ is a carrying-capacity factor.

- Population updates:
  $$N'_{\text{adult}} = N_{\text{adult}} + M - D_{\text{adult}}$$
  $$N'_{\text{juvenile}} = N_{\text{juvenile}} + B - D_{\text{juvenile}} - M$$

The executable model then projects these values year-by-year for the selected species and forecast horizon.

## Planned next steps

1. **Connect to real datasets**
   - replace placeholder data with real species survey or monitoring data
   - normalize incoming data to a common schema
   - support alternate input structures for age, sex, and population counts

2. **Add richer demographic structure**
   - separate age classes beyond adult/juvenile
   - track sex-specific populations and fertility
   - model stage transitions explicitly rather than a single maturation term

3. **Incorporate predation and density dependence**
   - add predator-driven mortality or harvest pressure
   - allow resource limits to affect birth and survival rates
   - support Allee effects or density-dependent reproduction

4. **Move toward spatial / GIS-aware modelling**
   - integrate spatial data for habitat quality and resource availability
   - map carrying capacity across cells or landscape units
   - support movement, dispersal, and local density dynamics

## Getting started

1. Create and activate a Python virtual environment.
2. Install requirements:
   ```bash
   pip install -r requirements.txt
   ```
3. Run the model:
   ```bash
   ./growth_model.py --species_name Zebra --start_year 2024 --forecast_horizon 10
   ```

## Notes

- This repository is a lightweight prototype and should be extended before use in production analysis.
- The current model is intentionally simple to allow future enhancement toward real datasets, age structure, sex-specific dynamics, predation, and spatial modelling.
