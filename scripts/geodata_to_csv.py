import csv
import pandas as pd
import yaml
from pathlib import Path

# Get ROOT_DIR (as in definitions.py)
ROOT_DIR = Path(__file__).parents[1]

# Get dir paths + filenames for I/O
with open(ROOT_DIR.joinpath("config.yaml"), "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

# Geospatial data dirs
GEOCLIM_DIR = ROOT_DIR.joinpath("data") / config["geoclim"]["dir"]["main"]
GPS_DATA_DIR = GEOCLIM_DIR / config["geoclim"]["dir"]["gps_data"]

# file names
raw_in = config["geoclim"]["files"]["raw_table"]
csv_out = config["geoclim"]["files"]["gps_in"]

# read excel file
df = pd.read_excel(GPS_DATA_DIR / raw_in)
print(df)
