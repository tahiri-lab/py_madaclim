import pathlib
from pathlib import Path
import yaml

class Definitions:
    """
    A class to manage configurations for a project.

    This class reads configuration data from a YAML file and provides
    access to the data through its attributes and methods.

    Attributes:
        root_dir: A pathlib.Path object representing the root directory of the project.
        package_dir: A pathlib.Path object representing the package directory.
        config_filename: The name of the configuration file.
        config: A dictionary containing the configuration data.
        data_dir: The path to the data directory.
        scripts_dir: The path to the scripts directory.
        geoclim_dirs: A dictionary containing paths to directories related to geoclim data.
        geoclim_files: A dictionary containing paths to files related to geoclim data.
        genetic_dirs: A dictionary containing paths to directories related to genetic data.
        genetic_files: A dictionary containing paths to files related to genetic data.
        urls: A dictionary containing URLs for Madaclim database.
        params: A dictionary containing tunable parameters.

    >>> from coffeaphylogeo.definitions import Definitions
    
    # Fetch data from the config.yaml file
    >>> defs = Definitions()
    >>> defs.geoclim_dirs
    {'main': 'geoclim', 'gps_data': 'geospatial', 'climate_data': 'climate', 'environment_data': 'enviro'}
    
    # Change an attribute on the fly
    >>> defs.params
    {'bp_to_keep': 1000}
    >>> defs.params["bp_to_keep"] = 500
    >>> defs.params
    {'bp_to_keep': 500}

    """
    def __init__(self, config_filename: str="config.yaml"):
        """
        Initializes a Definitions object.

        Args:
            config_filename: The name of the configuration file. Defaults to "config.yaml".
        """
        
        # Main project and package dirs
        self.root_dir = Path(__file__).parents[2]
        self.package_dir = Path(__file__).parent
        
        # Get configuration yaml file data
        self.config_filename = config_filename
        self.config = self.read_config(self.package_dir, self.config_filename)

        # Data and scripts dir
        self.data_dir = self.select_nested_dicts(config=self.config, most_upstream_key="outside_dirs", single_key_value="data")
        self.scripts_dir = self.select_nested_dicts(config=self.config, most_upstream_key="outside_dirs", single_key_value="scripts")

        # Geoclim-related dirs and files
        self.geoclim_dirs = self.select_nested_dicts(config=self.config, most_upstream_key="geoclim", single_key_value="dir")
        self.geoclim_files = self.select_nested_dicts(config=self.config, most_upstream_key="geoclim", single_key_value="files")

        # Genetic-related dirs and files
        self.genetic_dirs = self.select_nested_dicts(config=self.config, most_upstream_key="genetic", single_key_value="dir")
        self.genetic_files = self.select_nested_dicts(config=self.config, most_upstream_key="genetic", single_key_value="files")

        # URLs for Madaclim db
        self.urls = self.select_nested_dicts(config=self.config, most_upstream_key="urls")

        # Tunable parameters
        self.params = self.select_nested_dicts(config=self.config, most_upstream_key="params")

    # Read the YAML configuration file
    def read_config(self, package_dir: pathlib.Path, config_filename: str) -> dict:
        """
        Reads a YAML configuration file.

        Args:
            package_dir: The directory containing the configuration file.
            config_filename: The name of the configuration file.

        Returns:
            A dictionary containing the configuration data.
        """
        with open(package_dir / config_filename, "r") as yamlfile:
            config = yaml.safe_load(yamlfile)
        return config
    
    def select_nested_dicts(self, config: dict, most_upstream_key: str, single_key_value: str=None) -> dict:
        """
        Selects nested dictionaries from a configuration dictionary.

        Args:
            config: A dictionary containing the configuration data.
            most_upstream_key: The key of the top-level dictionary to select.
            single_key_value: The key of a single entry to select from the nested dictionary. If None, the entire nested dictionary is returned.

        Returns:
            A dictionary containing the selected data.
        """
        # Nested dicts under the most upstream key
        nested_dicts = {k: v for k, v in config.items() if k == most_upstream_key}
        
        # Get all dictionary or a single entry from nested dict
        if single_key_value is None:
            return nested_dicts[most_upstream_key]
        else:
                return nested_dicts[most_upstream_key][single_key_value]
        
    def get_geoclim_path(self, subdirectory: str) -> pathlib.Path:
        """
        Get the path to the desired subdirectory within the geoclim datasets

        Args:
            subdirectory (str): The name of the subdirectory to retrieve.

        Returns:
            pathlib.Path: Path to the specified geoclim subdirectory.
            
        Raises:
            ValueError: If the specified subdirectory is not a valid geoclim subdirectory.
            FileNotFoundError: If the generated path does not exists.
        """
        # Validate directory walk for geoclim-related data
        possible_subdirectories = [key for key in self.geoclim_dirs.keys()]
        if subdirectory not in possible_subdirectories:
            raise ValueError(f"Datatype directory must be one of {possible_subdirectories}")
        
        # Retrieve main geoclim
        if subdirectory == "main":
            path = self.root_dir / self.data_dir / self.geoclim_dirs[subdirectory]
            if path.exists():
                return path
            else:
                FileNotFoundError(f"{path} does not exists.")
        
        # Retrieve geoclim subdirs
        else:
            path = self.root_dir / self.data_dir / self.geoclim_dirs["main"] / self.geoclim_dirs[subdirectory]
            if path.exists():
                return path
            else:
                FileNotFoundError(f"{path} does not exists.")
    
    def get_genetic_path(self, subdirectory: str, concat: bool=False) -> pathlib.Path:
        """
        Get the path to the desired subdirectory within the genetic datasets
        
        Args:
            subdirectory (str): The name of the subdirectory to retrieve.
            concat (bool): If True, retrieve the concatenated files subdirectory. Defaults to False.
        
        Returns:
            pathlib.Path: The path to the specified genetic subdirectory.
        
        Raises:
            ValueError: If the specified subdirectory is not a valid genetic subdirectory.
            FileNotFoundError: If the generated path does not exist.
        """
        
        # Validate directory walk for genetic-related data
        possible_subdirectories = [key for key in self.genetic_dirs.keys() if not key.endswith("_concat")]
        if subdirectory not in possible_subdirectories:
            raise ValueError(f"Datatype directory must be one of {possible_subdirectories}")
        
        # Retrieve concatened files subdirectory (lower depth)
        if concat:
            path = self.root_dir / self.data_dir / self.genetic_dirs["main"] / self.genetic_dirs[subdirectory] / self.genetic_dirs[subdirectory + "_concat"]
            
            if path.exists():
                return path
            else:
                raise FileNotFoundError(f"{path} does not exists.")
        
        # Retrieve main genetic
        if subdirectory == "main":
            path = self.root_dir / self.data_dir / self.genetic_dirs[subdirectory]
            if path.exists():
                return path
            else:
                raise FileNotFoundError(f"{path} does not exists.")
        
        # Retrieve other genetic subdirs
        else:
            path = self.root_dir / self.data_dir / self.genetic_dirs["main"] / self.genetic_dirs[subdirectory]
            if path.exists():
                return path
            else:
                raise FileNotFoundError(f"{path} does not exists.")


    @property
    def config_filename(self):
        """The name of the configuration file."""
        return self._config_filename
    
    @config_filename.setter
    def config_filename(self, value):
        """
        Sets the name of the configuration file and updates the config attribute.

        Args:
            value: The new name of the configuration file.

        Raises:
            TypeError: If value is not a string.
        """
        if not isinstance(value, str):
            raise TypeError(f"Configuration filename must be a string (config_filename={value} is of type({type(value)})")
        self._config_filename = value

        # Update the config attribute
        self.config = self.read_config(self.package_dir, self._config_filename)
