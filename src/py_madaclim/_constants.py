from pathlib import Path
import importlib.resources as pkg_resources
import pyproj

class Constants:
    # Define relevant directories
    MAIN_DIR = Path(pkg_resources.files("py_madaclim"))
    _DATA_DIR = Path(pkg_resources.files("py_madaclim"), "_data")

    # check directories
    assert _DATA_DIR.is_dir(), f"Directory not found: {_DATA_DIR}"
    assert MAIN_DIR.is_dir(), f"Directory not found: {MAIN_DIR}"

    # Define dataformat and metadata file paths
    CLIM_DATAFORMAT_FILE = _DATA_DIR / "clim_data_format.json"
    CLIM_METADATA_FILE = _DATA_DIR / "clim_metadata.json"
    ENV_DATAFORMAT_FILE = _DATA_DIR / "env_data_format.json"
    ENV_METADATA_FILE = _DATA_DIR / "env_metadata.json"

    # check files
    assert CLIM_DATAFORMAT_FILE.is_file(), f"File not found: {CLIM_DATAFORMAT_FILE}"
    assert CLIM_METADATA_FILE.is_file(), f"File not found: {CLIM_METADATA_FILE}"
    assert ENV_DATAFORMAT_FILE.is_file(), f"File not found: {ENV_DATAFORMAT_FILE}"
    assert ENV_METADATA_FILE.is_file(), f"File not found: {ENV_METADATA_FILE}"

    # Default raster filenames
    DEFAULT_CLIM_RASTER_FILENAME = "madaclim_current.tif"
    DEFAULT_ENV_RASTER_FILENAME = "madaclim_enviro.tif"
    # Raster metadata to avoid instantiation of MadaclimRasters
    MADACLIM_CRS = pyproj.CRS.from_epsg(32738)
    MADACLIM_RASTERS_BOUNDS = (298000.0, 7155000.0, 1101000.0, 8683000.0)

    # URLs
    MADACLIM_URLS = {
        "clim_raster": "http://madaclim.cirad.fr/climate/current.tif",
        "env_raster": "http://madaclim.cirad.fr/environ/environ.tif",
        "clim_methods": "https://madaclim.cirad.fr/methods-climate/",
        "env_methods": "https://madaclim.cirad.fr/methods-environ-data/",
        "clim_data_format": "https://madaclim.cirad.fr/climate-data-format/",
        "env_data_format": "https://madaclim.cirad.fr/environ-data-format/"
    }

    # .env required keys for GBIF API auth for download request
    GBIF_DOTENV_KEYS = ["GBIF_USERNAME", "GBIF_PASSWORD"]

    # GBIF API base urls
    GBIF_BASEURLS = {
        "occurrence": "https://api.gbif.org/v1/occurrence",
        "species": "https://api.gbif.org/v1/species"
    }

    # Download format in occurence/download request endpoint
    GBIF_DOWNLOAD_FORMAT = ["SIMPLE_CSV", "DWCA", "SPECIES_LIST"]

    GBIF_MATCH_PARAMS = [
        "rank", 
        "strict", 
        "verbose", 
        "kingdom", 
        "phylum", 
        "class", 
        "order", 
        "family", 
        "genus"
    ]