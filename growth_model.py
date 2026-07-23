# !/usr/bin/env python3

from typing import Any, Dict, Tuple


def adult_pop_change(juv_maturation_adult: float, adult_death_rate: float) -> float:
    """
    Calculate the change in adult population based on the current adult and juvenile populations,
    carrying capacity, and adult death rate.

    Args:
        juv_maturation_adult (float): The rate at which juveniles mature into adults.
        adult_death_rate (float): The rate at which adults die.

    Returns:
        float: The change in adult population.
    """
    return juv_maturation_adult - adult_death_rate

def juv_pop_change(birth_rate: float, juv_death_rate: float, juv_maturation_adult: float) -> float:
    """
    Calculate the change in juvenile population based on the birth rate, juvenile death rate,
    and the rate at which juveniles mature into adults.

    Args:
        birth_rate (float): The rate at which new juveniles are born.
        juv_death_rate (float): The rate at which juveniles die.
        juv_maturation_adult (float): The rate at which juveniles mature into adults.

    Returns:
        float: The change in juvenile population.
    """

    return birth_rate - juv_death_rate - juv_maturation_adult

def calc_carrying_capacity(species_name: str) -> float:
    """
    Calculate the carrying capacity for a given species based on its average weight and the total biomass.

    Args:
        species_name (str): The name of the species for which to calculate carrying capacity.

    Returns:
        float: The carrying capacity for the specified species.
    """
    total_kg = 0
    species_df = get_species_df()  # Assuming this function retrieves a list of species objects
    for _, species in species_df.iterrows():
        species_pop = species.population
        total_kg += species_pop * species.avg_weight_kg  # Assuming each species object has an avg_weight_kg attribute

    species_avg_weight_kg = species_df[species_df['name'] == species_name]['avg_weight_kg'].values[0]
    return total_kg / species_avg_weight_kg

def get_species_df():
    """
    Placeholder function to retrieve a list of species objects.
    Replace this with actual implementation to fetch species data.

    Returns:
        list: A list of species objects.
    """
    # Example implementation, replace with actual data retrieval logic
    return []

def get_species_details(species_name: str, year: int) -> Dict[str, Any]:
    """
    Placeholder function to retrieve details for a specific species.
    Replace this with actual implementation to fetch species details.

    Args:
        species_name (str): The name of the species for which to retrieve details.
        year (int): The year for which to retrieve species details.

    Returns:
        dict: A dictionary containing species details.
    """
    # Example implementation, replace with actual data retrieval logic
    species_details: Dict[str, Any] = {
        "name": species_name,
        "avg_weight_kg": 10.0,  # Example average weight
        "adult_pop": 1000,      # Example population
        "juvenile_pop": 500,     # Example juvenile population
        "total_pop": 1500,        # Example total population
        "gestation_period_months": 6,  # Example gestation period
        "life_span_years": 5,    # Example life span
        "sex_matur_m_years": 2,  # Example male sexual maturity age
        "sex_matur_f_years": 2,  # Example female sexual maturity age
        "sex_ratio": 0.5         # Example sex ratio        
    }
    return species_details

def population_change(
    species_details: Dict[str, Any],
    carrying_capacity: float = None
    ) -> Tuple[float, float]:
    """
    Calculate the change in population in a year for a given species based on its birth rate, death rate, and carrying capacity.

    Args:
        species_details (Dict[str, Any]): A dictionary containing details about the species.
        carrying_capacity (float, optional): The carrying capacity for the species. Defaults to None.

    Returns:
        Tuple[float, float]: A tuple containing the new adult and juvenile populations.
    """
    # Calculate the birth and death rates based on species details
    birth_rate = 12 / species_details["gestation_period_months"]
    death_rate = 1 / species_details["life_span_years"]

    # Calculate the rate the avg rate of sexual maturation for both males and females, adjusted for the sex ratio
    sex_mat_rate_m = 1 / species_details["sex_matur_m_years"]
    sex_mat_rate_f = 1 / species_details["sex_matur_f_years"]
    sex_mat_rate_years = (sex_mat_rate_m + sex_mat_rate_f) * species_details["sex_ratio"]

    # Calculate the rate at which juveniles mature into adults, adjusted for carrying capacity
    juv_maturation_adult = sex_mat_rate_years * (1 - (species_details["adult_pop"]+species_details["juvenile_pop"]) / carrying_capacity)

    # Calculate the changes in adult and juvenile populations
    adult_pop_change = adult_pop_change(juv_maturation_adult=juv_maturation_adult, adult_death_rate=death_rate)
    juv_pop_change = juv_pop_change(birth_rate=birth_rate, juv_death_rate=death_rate, juv_maturation_adult=juv_maturation_adult)
    
    new_juvenile_pop = species_details["juvenile_pop"] + juv_pop_change
    new_adult_pop = species_details["adult_pop"] + adult_pop_change

    return new_adult_pop, new_juvenile_pop

def main():
    species_name = "ExampleSpecies"
    start_year = 2024
    forecast_horizon = 10  # Number of years to forecast

    species_details = get_species_details(species_name, start_year)
    carrying_capacity = calc_carrying_capacity(species_name=species_name)

    species_population_per_year = {}

    for year in range(start_year, start_year + forecast_horizon):
        new_adult_pop, new_juvenile_pop = population_change(species_details, carrying_capacity)
        species_population_per_year[year] = (new_adult_pop, new_juvenile_pop)
        species_details["adult_pop"] = new_adult_pop
        species_details["juvenile_pop"] = new_juvenile_pop
        species_details["total_pop"] = new_adult_pop + new_juvenile_pop

    print(f"Population forecast for {species_name}:")
    for year, (adult_pop, juv_pop) in species_population_per_year.items():
        print(f"  {year}: Adults = {adult_pop}, Juveniles = {juv_pop}")