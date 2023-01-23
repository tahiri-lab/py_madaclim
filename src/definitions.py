# Get cross-system compatible paths using pathlib
from pathlib import Path

# Main project/root folder
ROOT_DIR = Path(__file__).parents[1]

# Import YAML config variables + params
import yaml
with open(ROOT_DIR / "config.yaml", "r") as yf:
    config = yaml.safe_load(yf)

# Data-relevant paths
SNP_DIR = ROOT_DIR.joinpath("data", "GBS")
GPS_DIR = ROOT_DIR.joinpath("data", "geospatial")
SNP_TRIM_DIR = ROOT_DIR.joinpath("data", config["trimmed_dirs"]["GBS_in"])

# Module paths
GENETIC_DIR = ROOT_DIR.joinpath("src", "genetic")
CLIMATE_DIR = ROOT_DIR.joinpath("src", "climate")
