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
GEOSPATIAL_DIR = GEOCLIM_DIR / config["geoclim"]["dir"]["gps_data"]

# File names
METADATA_RAW = config["geoclim"]["files"]["raw_table"]
GPS_ALL = config["geoclim"]["files"]["gps_all"]
GPS_GBS_ONLY = config["geoclim"]["files"]["gbs_filtered"]
MADACLIM_CURRENT_TIF = config["geoclim"]["files"]["madaclim_current"]

if __name__ == "__main__":
    # Read from raw table
    df = pd.read_excel(GEOSPATIAL_DIR / METADATA_RAW, decimal=",")

    # Remove blank rows
    df = df.dropna(how="all")
    df["GBS sequence"] = df["GBS sequence"].str.rstrip()    #* \n still present after stripping

    # Remove \n in Botanical series
    df = df.replace(r"\n","", regex=True)

    # Converting to proper NaN and reordering
    df = df.replace({"-": np.nan, "NO POSITION": np.nan})
    df = df[["Species", "Species code", "Population code", "GBS sequence", "Botanical series", "Genome size (2C. pg)", "Latitude", "Longitude"]]
    # Remove NO from AND2 species
    df["Longitude"] = df["Longitude"].str.replace("NO ", "")
    df["Longitude"] = df["Longitude"].str.replace(",", ".").astype(float)

    # CSV for all species WITH GPS positions available
    df_gps_all = df
    df_gps_all = df_gps_all.dropna(subset=["Latitude", "Longitude"])
    print(f"Number of samples for all specimens WITH GPS : {len(df_gps_all)}")
    df_gps_all = df_gps_all.reset_index().drop(columns=["index"])
    df_gps_all.to_csv(GEOSPATIAL_DIR / GPS_ALL, index=False)

    # CSV for species with BOTH sequencing + GPS data ONLY
    df_gbs_only = df_gps_all.dropna()
    print(f"Number of samples for specimens with BOTH GPS AND genetic data : {len(df_gbs_only)}")
    df_gbs_only = df_gbs_only.reset_index().drop(columns=["index"])
    df_gbs_only.to_csv(GEOSPATIAL_DIR / GPS_GBS_ONLY, index=False)