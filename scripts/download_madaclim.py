import yaml
from pathlib import Path
import subprocess
import requests
import argparse

ROOT_DIR = Path(__file__).parents[1]
SRC_DIR = ROOT_DIR / "src"
PACKAGE_DIR = SRC_DIR / "coffeaphylogeo"

with open(PACKAGE_DIR / "config.yaml", "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

# Madaclim current climate and environmental data dirs/urls/names
CLIM_DATA_DIR = ROOT_DIR.joinpath("data") / config["geoclim"]["dir"]["main"] / config["geoclim"]["dir"]["climate_data"]
ENVIRO_DATA_DIR = ROOT_DIR / "data" / config["geoclim"]["dir"]["main"] / config["geoclim"]["dir"]["environment_data"]

current_madaclim_url = config["urls"]["madaclim_current"]
enviro_madaclim_url = config["urls"]["environment"]

current_madaclim_tif = config["geoclim"]["files"]["madaclim_current_raster"]
enviro_madaclim_tif = config["geoclim"]["files"]["madaclim_enviro_raster"]

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

def download_file(url, local_filename):
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return local_filename

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Download Madaclim data.')
    parser.add_argument('--wget', action='store_true', help='Use wget to download the file')
    args = parser.parse_args()

    if not CLIM_DATA_DIR.exists():
        CLIM_DATA_DIR.mkdir()
    if not (CLIM_DATA_DIR / current_madaclim_tif).is_file():
        if args.wget:
            runcmd(f"wget -P {CLIM_DATA_DIR} -O {current_madaclim_tif} {current_madaclim_url}", verbose=True)
        else:
            print(f"Downloading {current_madaclim_tif}...")
            download_file(current_madaclim_url, CLIM_DATA_DIR / current_madaclim_tif)
            print("Done")
    else:
        print(f"{current_madaclim_tif} file already exists.")

    if not (ENVIRO_DATA_DIR / enviro_madaclim_tif).is_file():
        if args.wget:
            runcmd(f"wget -P {ENVIRO_DATA_DIR} -O {enviro_madaclim_tif} {enviro_madaclim_url}", verbose=True)
        else:
            print(f"Downloading {enviro_madaclim_tif}...")
            download_file(enviro_madaclim_url, ENVIRO_DATA_DIR / enviro_madaclim_tif)
            print("Done")
    else:
        print(f"{enviro_madaclim_tif} file already exists.")

