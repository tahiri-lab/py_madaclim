import pathlib
from pathlib import Path
import json 
import yaml
from typing import List

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

# urls containing metadata
url_clim = config["urls"]["clim_methods"]
url_env = config["urls"]["env_methods"]

class WebPageError(Exception):
    pass

# Read metadata tables from madaclim webpage
def get_metadata_tables(url: str) -> List[pd.DataFrame]:
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
    print("Displaying first table retrieved...")
    print(tables[0])
    
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
        data[f"table{i}_{table.columns[0]}"] = table.to_json()
    
    # Dump to single output file
    with open(outdir / filename, "w") as jsonfile:
        json.dump(obj=data, fp=jsonfile, indent=4)


if __name__ == "__main__":
    clim_tables = get_metadata_tables(url_clim)
    save_to_json(clim_tables, CLIMATE_DIR, "test.json")