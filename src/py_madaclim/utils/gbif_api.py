import pathlib
from pathlib import Path
import os
from zipfile import ZipFile
from typing import Optional, Union
import re
import requests
import requests.exceptions
import time

from dotenv import load_dotenv, dotenv_values
from dwca.read import DwCAReader
import pandas as pd

from py_madaclim._constants import Constants

def get_taxon_key_by_species_match(
    name: Optional[str] = None, 
    return_full_on_match=False, 
    **kwargs: Union[str, bool]
) -> Union[dict, int]:
    """
    Fetches a taxon key by matching species name and other parameters.

    Args:
        name (str, optional): The scientific name of the species to match. If not provided, at least one other match parameter must be provided.
        return_full_on_match (bool, optional): If set to True, returns the full match data. If False, returns only the taxon key.
        **kwargs (Union[str, bool]): Additional match parameters. These could include 'rank', 'strict', 'verbose', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus'.

    Returns:
        dict or int: If return_full_on_match is True, returns a dictionary with full match data. Otherwise, returns the taxon key (integer).

    Raises:
        TypeError: If 'name' is provided and is not a str, or if 'strict' is provided and is not a bool, or if any other kwargs are provided and are not str.
        ValueError: If an invalid argument is provided, or if 'name' is not provided and no other match parameters are provided.
        Exception: If an error occurs while making the HTTP request.
    """
    if name is not None and not isinstance(name, str):
        raise TypeError(f"'name' must be a str if specified, but received {type(name).__name__}")

    params = {k: v for k, v in kwargs.items()}
    possible_params = Constants.GBIF_MATCH_PARAMS

    for param, value in params.items():
        if param not in possible_params:
            raise ValueError(f"Argument '{param}' is not valid. It must be one of {possible_params}")
        if param == "strict" and not isinstance(value, bool):
            raise TypeError(f"The argument 'strict' must be of type bool, but received {type(value).__name__}")
        if param != "strict" and not isinstance(value, str):
            raise TypeError(f"Keyword arguments values must be of type str, but received {type(value).__name__} for {param}")

    if name is not None:
        params["name"] = name

    if len(params) == 0:
        raise ValueError(f"Provide at least one match param argument from {possible_params} if 'name' is not provided.")

    url = Constants.GBIF_BASEURLS["species"] + "/match"
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
    except requests.RequestException as e:
        raise Exception(f"An error occurred while making the request: {e}")

    try:
        data = response.json()
    except:
        raise ValueError("Error: Response is not a valid JSON.")

    if data["matchType"] == "NONE":
        print("No matches found.")
    else:    
        print(
            f"{data['matchType']} match type found with {data['confidence']}% confidence!\n"
            f"canonical name of match: {data['canonicalName']}\nGBIF_taxon_key: {data['usageKey']}"
        )
    if return_full_on_match:
        return data 
    elif data["matchType"] == "NONE":
        return None
    else:
        return data["usageKey"]

def search_occ_mdg_valid_coordinates(taxon_key: int, year_range: Optional[tuple[int]]=None) -> list:
    """
    Searches for occurrences with valid coordinates for a given taxon key using the GBIF API '/occurrence/search' endpoint. 

    Uses a set of predetermined search params for the Madagascar GADM geographic identifier and geographic coordinates requirements.    
    Provides limited customization for taxon key and year range.
        
    Args:
        taxon_key (int): The GBIF taxon key to search for.
        year_range (tuple[int], optional): A tuple specifying the range of years to search for occurrences in.

    Returns:
        list: A list of occurrences with valid coordinates for the given taxon key.

    Raises:
        TypeError: If 'taxon_key' is not an integer, or if 'year_range' is provided and is not a tuple of integers.
        ValueError: If 'year_range' is provided and does not contain exactly 2 elements, or if the first element of 'year_range' is larger than the second.
        Exception: If an error occurs while making the HTTP request.
    """
    if not isinstance(taxon_key, int):
        raise TypeError(f"'taxon_key' must be an int, but received {type(taxon_key).__name__}")

    if year_range is not None:
        if not isinstance(year_range, tuple):
            raise TypeError("'year_range' must be a tuple")
        if len(year_range) != 2:
            raise ValueError("'year_range' must contain exactly 2 elements.")
        if not all([isinstance(year, int) for year in year_range]):
            raise TypeError("all elements of year range must be integers")
        if year_range[0] > year_range[1]:
            raise ValueError("First year element in range must be smaller than the second.")
    
    all_results = []    # Occurences accross selected year range container

    url = Constants.GBIF_BASEURLS['occurrence'] + "/search"
    params = {
        "gadmGid": "MDG",
        "has_coordinate": "true",
        "has_geospatial_issue": "false",
        "taxon_key": taxon_key,
        "year": f"{year_range[0]},{year_range[1]}" if year_range is not None else "",
    }

    if year_range is None:    # Get all possible years of occurences when not specified
        params.pop("year")
    
    # Check for server and proper year range params
    try:
        response = requests.get(url=url, params=params)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while making the request: {e}")
        return
    try:
        data = response.json()
    except ValueError:
        print("Error: Response is not a valid JSON.")
        return
    
    # Warnings/info with valid response
    if data["count"] == 0:
        print("Could not retrieve any records in the specified year range")
        return
    else:
        if year_range:
            print(f"Fetching all {data['count']} occurences in year range {year_range[0]}-{year_range[1]}...")
        else:
            print(f"Fetching all {data['count']} occurences accross all possible years...")
        params["limit"] = 300    # Maximum ratelimit per request
    time.sleep(1)

    # Fetch all GBIF occurences of Madagascar coffea species within specified year range
    while True:
        response = requests.get(url=url, params=params)
        response.raise_for_status()
        data = response.json()

        min_limit = data["offset"]
        max_limit = data["offset"] + data["limit"] if (data["offset"] + data["limit"]) < data["count"] else data["count"]
        print(f"Extracting occurences {min_limit} to {max_limit}...")
        all_results.extend(data["results"])
        
        if data["endOfRecords"]:    
            break  
        
        params["offset"] = data["offset"] + data["limit"]  # set the offset to the next page
        time.sleep(1)    # Safety even with no specified ratelimit/time

    print(f"Total records retrieved: {len(all_results)}")
    
    return all_results

def search_occ_by_gbif_id(gbif_id: Union[str, int]) -> Union[dict, None]:
    """
    Gets details for a single occurrence using its gbifID (key of single record) from the GBIF API '/occurrence/search' endpoint.

    Args:
        gbif_id (Union[str, int]): The unique identifier for an occurrence record in GBIF.

    Returns:
        data (Union[dict, None]): If the gbif_id is valid, returns a dictionary containing the details of the record.
            Otherwise, it returns None.
    """

    # gbif_id conversion/typechecking
    try:
        gbif_id = str(gbif_id)
    except TypeError:
        print(f"Could not convert {gbif_id} to str.")
        return

    # Occurence/key request endpoint
    url = Constants.GBIF_BASEURLS["occurrence"] + f"/{gbif_id}"
    
    try:
        response = requests.get(url=url)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"An error occurred while making the request: {e}")
        return
    
    try:
        data = response.json()
    except ValueError:
        print("Error: Response is not a valid JSON.")
        data = None
    
    return data 

def request_occ_download_mdg_valid_coordinates(
    taxon_key: int,
    email: str,
    dotenv_filepath: Optional[Union[pathlib.Path, str]]=None,
    reponse_format: str="DWCA",
    year_range: Optional[tuple[int]]=None,
) -> Union[requests.models.Response, str]:
    """
    Makes a POST request to GBIF API to download occurrence data for a specific taxon key.

    Uses a set of predetermined payload information such as Madagascar GADM geographic identifier and geographic coordinates requirements.    
    Allows limited customization for taxon key and year range for the payload.

    Args:
        taxon_key (int): The taxon key to request occurrence data for.
        email (str): The email address to send the download link to.
        dotenv_filepath (str or pathlib.Path, optional): The path to the .env file containing GBIF credentials. Defaults to current directory.
        reponse_format (str, optional): The format of the download file. Choices are "DWCA", "SIMPLE_CSV", or "SPECIES_LIST". Defaults to "DWCA".
        year_range (tuple, optional): A tuple of two integers specifying the range of years to request data for. Defaults to None.

    Returns:
        requests.models.Response or str: The response object or the error message if an error occurred.

    Raises:
        TypeError: If 'taxon_key' is not an integer, or if 'dotenv_filepath' is not a valid type for pathlib.Path.
        ValueError: If 'email' is not a valid email address
        ValueError: If required keys are missing in '.env' file
        ValueError: If 'response_format' is not one of the allowed choices.
        FileNotFoundError: If '.env' file is not found.
        Exception: If an error occurred while making the request.
    """
    if not isinstance(taxon_key, int):
        raise TypeError(f"taxon_key must be an int, but received {type(taxon_key).__name__}")
    
    # email format validation
    email_pattern = r"\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,7}\b"
    if not re.match(email_pattern, email):
        raise ValueError(f"Invalid email: {email}")

    # .env checker and gbif credentials to environment vars
    if dotenv_filepath is None:
        dotenv_filepath = Path.cwd() / ".env"
    try:
        dotenv_filepath = Path(dotenv_filepath)
    except TypeError:
        raise TypeError("'dotenv_filepath' must be a valid type for the pathlib.Path constructor (str or system-specific pathlib.Path/PurePath)")
    
    if dotenv_filepath.is_file():
        print(".env file found.")
    else:
        raise FileNotFoundError(f"Could not find .env file: {dotenv_filepath}")
    
    dotenv_keys = Constants.GBIF_DOTENV_KEYS
    config = dotenv_values(dotenv_filepath)

    if not all([key in config.keys() for key in dotenv_keys]):
        raise ValueError(f"Missing credentials keys in '.env' file. The file must contain both {dotenv_keys} keys")
    
    for key, val in config.items():
        if val is None or val == "":
            raise ValueError(f"{key} cannot be null or empty value.")
    
    load_dotenv(dotenv_filepath)
    gbif_username = os.getenv(dotenv_keys[0])
    gbif_password = os.getenv(dotenv_keys[1])

    # Predicates checker
    if not all([isinstance(param, str) for param in [reponse_format, email]]):
        raise TypeError("'response_format' and 'email' must be of type str.")
    
    possible_formats = Constants.GBIF_DOWNLOAD_FORMAT
    if reponse_format not in possible_formats:
        raise ValueError(f"'reponse_format' must be one of {possible_formats}")

    if year_range is not None:
        if not isinstance(year_range, tuple):
            raise TypeError("'year_range' must be a tuple")
        if len(year_range) != 2:
            raise ValueError("'year_range' must contain exactly 2 elements.")
        if not all([isinstance(year, int) for year in year_range]):
            raise TypeError("all elements of year range must be integers")
        if year_range[0] > year_range[1]:
            raise ValueError("First year element in range must be smaller than the second.")
    
    # Query parameters for payload specific to mascarocoffea with valid coordinates
    url = Constants.GBIF_BASEURLS["occurence"] + "/download/request"
    predicates = [
        {"type": "equals", "key": "GADM_GID", "value": "MDG"},
        {"type": "equals", "key": "HAS_COORDINATE", "value": "true"},
        {"type": "equals", "key": "HAS_GEOSPATIAL_ISSUE", "value": "false"},
        {"type": "equals", "key": "TAXON_KEY", "value": taxon_key},
    ]

    if year_range is not None:    # Defaults to all occurrences accross all years
        predicates.append({"type": "equals", "key": "YEAR", "value": f"{year_range[0]},{year_range[1]}"})

    payload = {
        "creator": gbif_username,
        "notificationAddresses": [email],
        "predicate": {
            "type": "and",
            "predicates": predicates
        },
        "format": reponse_format
    }

    try:
        response = requests.post(url=url, json=payload, auth=(gbif_username, gbif_password))
        return response.text
    except requests.HTTPError:
        print(f"Request failed with status code {response.status_code}")
        return response
    
def download_extract_read_occ(
        download_id: str, 
        target_dir: Optional[Union[str, pathlib.Path]]=None
    ) -> Union[None, pd.DataFrame]:
    """
    Downloads, extracts and read the content of an occurrence data file from the GBIF API

    Args:
        download_id (str): The ID of the download to fetch
        target_dir (Optional[Union[str, pathlib.Path]], optional): The directory path to save the download and extracted files to. 
            Defaults to the current working directory

    Raises:
        TypeError: If the 'download_id' is not a string.
        TypeError: If the 'target_dir' is not a string or a pathlib.Path object.
        NotADirectoryError: If the provided argument for 'target_dir' is not a valid directory.
        Exception: If an error occurred while making the request.
        ValueError: If the JSON response could not be parsed.
        ValueError: If the 'format' object from the response is not valid.
        FileNotFoundError: If no '.csv' files are found in files from a non-DWCA download format.

    Returns:
        Union[None, pd.DataFrame]: A pandas dataframe containing all the occurrences data with columns according to the requested format.
            If data could not be properly extracted from the download, it returns None.
    """
    
    if not isinstance(download_id, str):
        raise TypeError(f"'download_id' must be a str, but received {type(download_id).__name__}.")

    # target_dir default and validation
    if target_dir is None:
        target_dir = Path.cwd()
    try:
        target_dir = Path(target_dir)
    except TypeError:
        raise TypeError("'target_dir' must be a valid type for the pathlib.Path constructor (str or system-specific pathlib.Path/PurePath)")
    
    if not target_dir.is_dir():
        raise NotADirectoryError(f"{target_dir} is not a valid directory.")
    
    url = Constants.GBIF_BASEURLS["occurrence"] + f"/download/{download_id}"

    try:
        response = requests.get(url=url)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        raise Exception(f"An error occurred while making the request: {e}")
    try:
        data = response.json()
    except ValueError:
        raise ValueError("Error: Response is not a valid JSON.")
    print(f"Response OK from endpoint {url.rsplit('/', 1)[0]} with {download_id=}.")

    # Get download properties for status bar and I/O name
    zipfile_name = Path(data["downloadLink"].split("/")[-1])
    total_size = data["size"]  # total_size is in bytes

    # Make a new request to get the actual file data
    response = requests.get(data["downloadLink"], stream=True)

    with open(target_dir / zipfile_name, 'wb') as f:
        chunk_size = round(1024 * 1024 / 10)    # approx 0.1MB chunks
        start_time = time.time()
        for n, chunk in enumerate(response.iter_content(chunk_size=chunk_size)):
            percent = (n * chunk_size / total_size) * 100
            now_time = time.time()
            current_speed = (n * chunk_size) / (now_time - start_time) / (1024 * 1024)  # current_speed is in MB/s
            print(
                f"Progress for {zipfile_name} : {percent:.1f}% completed of {total_size / (1024 * 1024):.2f} MB downloaded [ current speed of {current_speed:.2f} MB/s ]",
                end="\r"
            )
            f.write(chunk)
        end_time = time.time()
        print(
            f"Progress for {zipfile_name} : 100.0% completed of {total_size / (1024 * 1024):.2f} MB downloaded [ average speed of {total_size / (1024 * 1024) / (end_time - start_time):.2f} MB/s ]",
            end="\r"
        )
        print()
    
    # Unzipping file to target location
    with ZipFile(target_dir / zipfile_name, mode="r") as zObj:
        zipped_files = zObj.infolist()
        download_container = target_dir / f"download_{zipfile_name.stem}"
        Path.mkdir(download_container, exist_ok=True)    
        print(f"Extracting all {len(zipped_files)} files to target location: {download_container}...")
        zObj.extractall(download_container)

    # Read to pandas df according to download format
    download_format = data["request"]["format"]
    possible_formats = Constants.GBIF_DOWNLOAD_FORMAT
    if download_format not in possible_formats:
        raise ValueError(f"download_format from parsed HTTP response is not one of {possible_formats}")

    if download_format == "DWCA":
        with DwCAReader(target_dir / zipfile_name) as dwca:
            df = dwca.pd_read("occurrence.txt", parse_dates=True)
            print(f"Read and saved core data into pandas df: {dwca.descriptor.core.file_location}") # 'occurrence.txt'

    elif download_format == "SIMPLE_CSV" or download_format == "SPECIES_LIST":
        csv_files = []
        for file in download_container.iterdir():
            if file.suffix == ".csv":
                csv_files.append(file)
        if len(csv_files) == 0:
            raise FileNotFoundError(f"No '.csv' files in the extracted {zipfile_name}")
        elif len(csv_files) > 1:
            print(f"BEWARE! More than 1 '.csv' file in the extracted {zipfile_name}. Could not read into pandas df")
            df = None 
        else:
            print(f"Read and saved {csv_files[0].name} into pandas df")
            df = pd.read_csv(download_container / csv_files[0], delimiter="\t")
    return df