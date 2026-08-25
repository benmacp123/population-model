#!/usr/bin/env python3
"""Build the normalized candidate-species traits table from PanTHERIA."""

from pathlib import Path

import click
import pandas as pd

SOURCE_URL = "https://esapubs.org/archive/ecol/E090/184/"
TARGETS = {
    "Zebra": "Equus quagga",
    "Wildebeest": "Connochaetes taurinus",
    "Waterbuck": "Kobus ellipsiprymnus",
    "Kudu": "Tragelaphus strepsiceros",
}


def build_traits(input_path: Path, output_path: Path) -> pd.DataFrame:
    """Read PanTHERIA and write the model-facing candidate trait table."""
    source = pd.read_table(input_path).replace(-999, pd.NA)
    selected = source[source["MSW05_Binomial"].isin(TARGETS.values())].copy()
    selected["name"] = selected["MSW05_Binomial"].map(
        {binomial: common for common, binomial in TARGETS.items()}
    )

    traits = pd.DataFrame(
        {
            "name": selected["name"],
            "scientific_name": selected["MSW05_Binomial"],
            "avg_weight_kg": selected["5-1_AdultBodyMass_g"] / 1000,
            "gestation_period_months": selected["9-1_GestationLen_d"] / 30.44,
            "interbirth_interval_days": selected["14-1_InterbirthInterval_d"],
            "life_span_years": selected["17-1_MaxLongevity_m"] / 12,
            # PanTHERIA supplies one maturity age, so both sex fields share it
            # until a sex-specific source is added.
            "sex_matur_m_years": selected["23-1_SexualMaturityAge_d"] / 365.25,
            "sex_matur_f_years": selected["23-1_SexualMaturityAge_d"] / 365.25,
            "male_ratio": pd.NA,
            "juv_death_rate": pd.NA,
            "carrying_capacity": pd.NA,
            "source": SOURCE_URL,
        }
    ).sort_values("name")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    traits.to_csv(output_path, index=False, na_rep="")
    return traits


@click.command()
@click.option(
    "--input-path",
    type=click.Path(path_type=Path, exists=True, dir_okay=False),
    default=Path("data/raw/PanTHERIA_1-0_WR05_Aug2008.txt"),
    show_default=True,
)
@click.option(
    "--output-path",
    type=click.Path(path_type=Path, dir_okay=False),
    default=Path("data/raw/species_traits.csv"),
    show_default=True,
)
def main(input_path: Path, output_path: Path) -> None:
    traits = build_traits(input_path, output_path)
    click.echo(f"Wrote {len(traits)} species to {output_path}")


if __name__ == "__main__":
    main()
