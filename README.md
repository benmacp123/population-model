# simple-population-model

Simple population modelling based on maturation rates of individual species, sex-specific dynamics, and carrying capacity.

[![Python Applications test with GitHub Actions](https://github.com/benmacp123/simple-population-model/actions/workflows/main.yml/badge.svg)](https://github.com/benmacp123/simple-population-model/actions/workflows/main.yml)
## Current work

This project currently includes a small executable model in `growth_model.py` that:

- loads placeholder species data via `pandas`
- computes carrying capacity from total biomass and average weight
- projects adult and juvenile populations year-by-year
- applies birth, death, and maturation dynamics at the population level
- supports a CLI interface using `click`

The current implementation is designed to be a prototype framework rather than a final ecological model.

## What it currently does

- `get_species_df()` returns example species records for testing
- `calc_carrying_capacity()` computes a simple carrying-capacity proxy
- `population_change()` applies scaled birth, death, and maturation changes
- `run_growth_model()` performs the yearly forecast and prints results
- CLI options include `species_name`, `start_year`, and `forecast_horizon`

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
