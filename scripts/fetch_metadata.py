import pathlib
from pathlib import Path
import json 
import yaml
from typing import List
import argparse
import pandas as pd

# Get main project and package directory
ROOT_DIR = Path(__file__).parents[1]
SRC_DIR = ROOT_DIR / "src"

with open(SRC_DIR.joinpath("config.yaml"), "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

# Setup main data and target directories
DATA_DIR = ROOT_DIR / config["outside_dirs"]["data"]
GEOCLIM_DIR = DATA_DIR / config["geoclim"]["dir"]["main"]

CLIMATE_DIR = GEOCLIM_DIR / config["geoclim"]["dir"]["climate_data"]
ENVIRO_DIR =  GEOCLIM_DIR / config["geoclim"]["dir"]["environment_data"]

# URLs containing metadata
CLIM_META_URL = config["urls"]["clim_methods"]
CLIM_FORMAT_URL = config["urls"]["clim_data_format"]
ENV_META_URL = config["urls"]["env_methods"]
ENV_FORMAT_URL = config["urls"]["env_data_format"]


# Filenames for metadata json
CLIM_METADATA = config["geoclim"]["files"]["clim_metadata"]
CLIM_DATA_FORMAT = config["geoclim"]["files"]["clim_data_format"]
ENV_METADATA = config["geoclim"]["files"]["env_metadata"]
ENV_DATA_FORMAT = config["geoclim"]["files"]["env_data_format"]

class WebPageError(Exception):
    pass

# Read metadata tables from madaclim webpage
def get_data_tables(url: str) -> List[pd.DataFrame]:
    """Retrieves tables from a webpage and returns them as a list of DataFrames.

    Args:
        url (str): URL associated with webpage containing the metadata

    Raises:
        WebPageError: If failing to retrieve html tables 

    Returns:
        List[pd.DataFrame]: A list of DataFrames, where each DataFrame represents a table from the webpage.
    """
    try:
        tables = pd.read_html(url)
    except Exception as e:
        raise WebPageError(f"Failed to read HTML tables from {url} \nError output: {e}")
    
    print(f"Retrieved {len(tables)} table(s) from {url}")
    
    return tables

def save_to_json(tables: List[pd.DataFrame], outdir: pathlib.PosixPath, filename: str):
    """Saves a list of DataFrames to a single JSON file.

    Args:
        tables (List[pd.DataFrame]): A list of DataFrames to save.
        outdir (pathlib.PosixPath): The directory to save the JSON file in.
        filename (str): The name of the JSON file.

    """
    # Initialize dict
    data = {}

    # Save each table in dict as json
    for i, table in enumerate(tables):
        data[f"table_{i}}"] = table.to_json()
    
    # Dump to single output file
    with open(outdir / filename, "w") as jsonfile:
        json.dump(obj=data, fp=jsonfile, indent=4)

    print(f"Saved {filename} to {outdir}")

def build_arg_parser():
    """Builds and returns an argument parser for fetching metadata from the Madaclim database.

    This function creates an argument parser with a `--metadata` argument that allows the user to
    specify which metadata to fetch from the Madaclim database (https://madaclim.cirad.fr/). The
    `--metadata` argument accepts one of three values: "clim" to fetch climate data, "env" to fetch
    environmental data, or "all" to fetch both. If the user does not specify a value for the
    `--metadata` argument, the default value is "all".

    Returns:
        argparse.Namespace: An object containing the values for the defined arguments.
    """
    parser = argparse.ArgumentParser(
        description=
        """
        Fetch the metadata for both the climate and environmental data from the
        Madaclim database (https://madaclim.cirad.fr/)
        """
    )

    parser.add_argument(
        "--metadata",
        choices=["clim", "env", "all"],
        default="all",
        help="Select which of the metadata you want to fetch"
    )
    
    # Parse each arg and return the args objects
    args = parser.parse_args()
    return args

def main():
    """
    Fetches and saves metadata and data format information from the Madaclim database.

    This function fetches metadata and data format information from the Madaclim database
    (https://madaclim.cirad.fr/) based on the value of the `metadata_choice` variable, which is
    determined by the `--metadata` command line argument. The fetched data is then saved to JSON
    files in the specified directories.

    The `--metadata` argument accepts one of three values: "clim" to fetch climate data, "env" to
    fetch environmental data, or "all" to fetch both. If the user does not specify a value for the
    `--metadata` argument, the default value is "all".
    """
    # Select which tables to fetch
    args = build_arg_parser()
    metadata_choice = args.metadata

    # Fetch current_climate data
    if metadata_choice == "clim":
        metadata_dfs = get_data_tables(CLIM_META_URL)
        dataformat_dfs = get_data_tables(CLIM_FORMAT_URL)
        
        # Associated file and dir names
        outdir = [CLIMATE_DIR]
        outfile_meta = CLIM_METADATA
        outfile_format = CLIM_DATA_FORMAT

    # Fetch environmental data
    elif metadata_choice == "env":
        metadata_dfs = get_data_tables(ENV_META_URL)
        dataformat_dfs = get_data_tables(ENV_FORMAT_URL)
        
        # Associated file and dir names
        outdir = [ENVIRO_DIR]
        outfile_meta = ENV_METADATA
        outfile_format = ENV_DATA_FORMAT

    # Fetch both clim/env data
    else:
        metadata_dfs = {"clim_dfs": get_data_tables(CLIM_META_URL), "env_dfs": get_data_tables(ENV_META_URL)}
        dataformat_dfs = {"clim_dfs": get_data_tables(CLIM_FORMAT_URL), "env_dfs": get_data_tables(ENV_FORMAT_URL)}
        
        # Multiple outdirs for each type of fetched data
        outdir = {"clim_dir": CLIMATE_DIR, "env_dir": ENVIRO_DIR}
        outfile_meta = {"clim_name": CLIM_METADATA, "env_name": ENV_METADATA}
        outfile_format = {"clim_name": CLIM_DATA_FORMAT, "env_name": ENV_DATA_FORMAT}

    
    # Convert retrieved dataframes to json file in respective data dirs according to metadata arg
    if len(outdir) == 1:
        save_to_json(tables=metadata_dfs, outdir=outdir[0], filename=outfile_meta)
        save_to_json(tables=dataformat_dfs, outdir=outdir[0], filename=outfile_format)

    else:
        save_to_json(tables=metadata_dfs["clim_dfs"], outdir=outdir["clim_dir"], filename=outfile_meta["clim_name"])
        save_to_json(tables=metadata_dfs["env_dfs"], outdir=outdir["env_dir"], filename=outfile_meta["env_name"])
        save_to_json(tables=dataformat_dfs["clim_dfs"], outdir=outdir["clim_dir"], filename=outfile_format["clim_name"])
        save_to_json(tables=dataformat_dfs["env_dfs"], outdir=outdir["env_dir"], filename=outfile_format["env_name"])



if __name__ == "__main__":
    main()
