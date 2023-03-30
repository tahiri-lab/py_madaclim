import yaml
import pandas as pd
from pathlib import Path
from typing import List

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


if __name__ == "__main__":
    get_metadata(url_clim)