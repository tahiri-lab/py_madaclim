import yaml
from pathlib import Path
import subprocess

# Get ROOT_DIR (as in definitions.py)
ROOT_DIR = Path(__file__).parents[1]

# Get dir paths + filenames for I/O
with open(ROOT_DIR.joinpath("config.yaml"), "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

CLIM_DATA_DIR = ROOT_DIR.joinpath("data") / config["geoclim"]["dir"]["main"] / config["geoclim"]["dir"]["climate_data"]
current_madaclim_url = config["URLs"]["madaclim_current"]

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
    pass

if __name__ == "__main__":
    if not CLIM_DATA_DIR.exists():
        CLIM_DATA_DIR.mkdir()
        runcmd(f"wget -P {CLIM_DATA_DIR} {current_madaclim_url}", verbose=True)

    else:
        runcmd(f"wget -P {CLIM_DATA_DIR} {current_madaclim_url}", verbose=True)
