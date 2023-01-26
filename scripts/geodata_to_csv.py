import pandas as pd
import yaml
from pathlib import Path
import numpy as np

# Get ROOT_DIR (as in definitions.py)
ROOT_DIR = Path(__file__).parents[1]

# Get dir paths + filenames for I/O
with open(ROOT_DIR.joinpath("config.yaml"), "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

# Geospatial data dirs
GEOCLIM_DIR = ROOT_DIR.joinpath("data") / config["geoclim"]["dir"]["main"]
GPS_DATA_DIR = GEOCLIM_DIR / config["geoclim"]["dir"]["gps_data"]

# File names
raw_in = config["geoclim"]["files"]["raw_table"]
gps_all = config["geoclim"]["files"]["gps_all"]
gbs_filtered = config["geoclim"]["files"]["gbs_filtered"]

if __name__ == "__main__":
    # Read from raw table
    df = pd.read_excel(GPS_DATA_DIR / raw_in, decimal=",")

    # Remove blank rows
    df = df.dropna(how="all")
    df["GBS sequence"] = df["GBS sequence"].str.rstrip()    #* \n still present after stripping

    # Remove \n in Botanical series
    df = df.replace(r"\n","", regex=True)

    # Converting to proper NaN and reordering
    df = df.replace({"-": np.nan})
    df = df[["Species", "Species code", "Population code", "GBS sequence", "Botanical series", "Genome size (2C. pg)", "Latitude", "Longitude"]]
    df.head()

    # CSV for all species WITH GPS positions available
    df_gps_all = df
    print(df_gps_all)
    df_gps_all.to_csv(GPS_DATA_DIR / gps_all, index=False)

    # CSV for species with BOTH sequencing + GPS data ONLY
    df_gbs_only = df_gps_all.dropna()
    df_gbs_only.to_csv(GPS_DATA_DIR / gbs_filtered, index=False)  