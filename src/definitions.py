# Get cross-system compatible paths using pathlib
from pathlib import Path

# Relative main folder path
ROOT_DIR = Path(__file__).parents[1]

# Package src dir
SRC_DIR = ROOT_DIR / "src"

# Import YAML config variables + params
import yaml
with open(SRC_DIR / "config.yaml", "r") as yf:
    config = yaml.safe_load(yf)

# Main data dir path
DATA_DIR = ROOT_DIR / "data"

# Genetic data dirs
GENETIC_DIR = DATA_DIR / config["genetic"]["dir"]["main"]
SNP_DIR = GENETIC_DIR / config["genetic"]["dir"]["gbs_full"]
SNP_TRIM_DIR = GENETIC_DIR / config["genetic"]["dir"]["gbs_trimmed"]

# Geoclimatic data dirs
GEOCLIM_DIR = DATA_DIR / config["geoclim"]["dir"]["main"]
GEOSPATIAL_DIR = GEOCLIM_DIR / config["geoclim"]["dir"]["gps_data"]
CLIMATE_DIR = GEOCLIM_DIR / config["geoclim"]["dir"]["climate_data"]

# Modules path
GENETIC_MOD = ROOT_DIR.joinpath("src", "genetic")
GEOCLIM_MOD = ROOT_DIR.joinpath("src", "climate")

# File names
METADATA_RAW = config["geoclim"]["files"]["raw_table"]
GPS_ALL = config["geoclim"]["files"]["gps_all"]
GPS_GBS_ONLY = config["geoclim"]["files"]["gbs_filtered"]
MADACLIM_CURRENT_TIF = config["geoclim"]["files"]["madaclim_current"]

# Madaclim climate variables
class MadaclimLayers:
    """
    #TODO docstring CLASS FOR MADACLIM DB PARAMS
    """
    def __init__(self):
        pass
    
if __name__ == "__main__":
    print(f"SNP path : {SNP_DIR}")
    print(MADACLIM_CURRENT_TIF)
    