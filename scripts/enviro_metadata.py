import yaml
import pandas as pd
from pathlib import Path

ROOT_DIR = Path(__file__).parents[1]
SRC_DIR = ROOT_DIR / "src"

with open(SRC_DIR.joinpath("config.yaml"), "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

# Target dir
ENVIRO_DIR = ROOT_DIR / "data" / config["geoclim"]["dir"]["main"] / config["geoclim"]["dir"]["environment_data"]

# Get Environmental metadata and save
enviro_meta = pd.read_html("https://madaclim.cirad.fr/environ-data-format/", index_col="Layer")[0]
enviro_meta.to_csv(ENVIRO_DIR / "enviro_metadata.csv")