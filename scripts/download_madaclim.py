import yaml
from pathlib import Path
import subprocess

# Get ROOT_DIR (as in definitions.py)
ROOT_DIR = Path(__file__).parents[1]

# Get dir paths + filenames for I/O
with open(ROOT_DIR.joinpath("config.yaml"), "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

# Madaclim current climate and environmental data
CLIM_DATA_DIR = ROOT_DIR.joinpath("data") / config["geoclim"]["dir"]["main"] / config["geoclim"]["dir"]["climate_data"]
current_madaclim_url = config["URLs"]["madaclim_current"]
enviro_madaclim_url = config["URLs"]["environment"]
current_madaclim_tif = config["geoclim"]["files"]["madaclim_current"]
enviro_madaclim_tif = config["geoclim"]["files"]["madaclim_enviro"]

# wget with subprocess
def runcmd(cmd, verbose = False, *args, **kwargs):

    process = subprocess.Popen(
        cmd,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        text = True,
        shell = True
    )
    std_out, std_err = process.communicate()
    if verbose:
        print(std_out.strip(), std_err)


if __name__ == "__main__":
    if not CLIM_DATA_DIR.exists():
        CLIM_DATA_DIR.mkdir()
        # Current climate db download
    if not (CLIM_DATA_DIR / current_madaclim_tif).is_file():
        runcmd(f"wget -P {CLIM_DATA_DIR} -O {current_madaclim_tif} {current_madaclim_url}", verbose=True)
    else:
        print(f"{current_madaclim_tif} file already exists.")
    
    # Environment db download
    #TODO IMPLEMENT PURE PYTHON DOWNLOAD
    # if not (CLIM_DATA_DIR / enviro_madaclim_tif).is_file():
    #     runcmd(f"wget -P {CLIM_DATA_DIR} -O {enviro_madaclim_tif} {enviro_madaclim_url}", verbose=True)
    # else:
    #     print(f"{enviro_madaclim_tif} file already exists.")
