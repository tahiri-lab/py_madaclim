from pathlib import Path

# Get user-specific directory path for the main project/root folder
ROOT_DIR = Path(__file__).parent

# Data-relevant paths
SNP_DIR = ROOT_DIR.joinpath("data", "GBS")
GPS_DIR = ROOT_DIR.joinpath("data", "geospatial")

# Module paths
GENETIC_DIR = ROOT_DIR.joinpath("src", "genetic")
CLIMATE_DIR = ROOT_DIR.joinpath("src", "climate")