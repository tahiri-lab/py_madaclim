import json
import re
import pathlib
from pathlib import Path
from typing import List

import pandas as pd

from coffeaphylogeo.definitions import Definitions

defs = Definitions()

class MadaclimLayers:
    """
    A class that represents all the names of the climate and environmental variable layers that can be found from the rasters of the Madaclim database.

    This class provides methods to generate the layers name, select subsets of layers and manipulate their formatting.
    """
    # Assign format and metadata dirs
    
    
    def __init__(self):
        """
        Initializes a new instance of the MadaclimLayers class.

        This constructor #TODO FINISH DOCSTRING
        """
        self.climate_dir = defs.get_geoclim_path("climate_data")    # Path dir for clim-related data
        self.enviro_dir = defs.get_geoclim_path("environment_data")    # Path dir for env-related data
        
        self.clim_data_file = defs.geoclim_files["clim_data_format"]
        self.clim_meta_file = defs.geoclim_files["clim_metadata"]
        self.env_data_file = defs.geoclim_files["env_data_format"]
        self.env_meta_file = defs.geoclim_files["env_metadata"]
        
        # self.all_layers = self._get_madaclim()
        # self.geology_raster_vars = {}

    def _get_madaclim(self, climate_dir, enviro_dir, clim_data_file, clim_meta_file, env_data_file, env_meta_file) -> pd.DataFrame :
        
        # Open data and metadata tables
        with open(climate_dir / clim_data_file, "r") as f:
            clim_format = json.load(f)
        with open(climate_dir /  clim_meta_file,"r") as f:
            clim_meta = json.load(f)
        with open(enviro_dir / env_data_file, "r") as f:
            env_format = json.load(f)
        with open(enviro_dir / env_meta_file, "r") as f:
            env_meta = json.load(f)

        # Extract climate data and format it using the metadata
        df_clim = pd.read_json(clim_format["table_0"])
        df_clim["data_type"] = "clim"    # Tag for latter merge

        
        


        
    
    # def __str__(self):
    #     all_layers = ""
    #     for layer_num, layer_name in self.madaclim_layers.items():
    #         all_layers+= f"{layer_num}: {layer_name}\n"
    #     return all_layers

    @property
    def climate_dir(self):
        """Path: The directory path for climate data."""
        return self._climate_dir
    
    @climate_dir.setter
    def climate_dir(self, value):
        """Sets the directory path for climate data.

        Args:
            value (Path): The directory path for climate data.

        Raises:
            ValueError: If the provided value is not a valid directory path.
        """
        if not value.is_dir():
            raise ValueError(f"{value} is not a valid directory path.")
        
        self._climate_dir = value
        
    @property
    def enviro_dir(self):
        """Path: The directory path for environmental data."""
        return self._enviro_dir
    
    @enviro_dir.setter
    def enviro_dir(self, value):
        """Sets the directory path for environmental data.

        Args:
            value (Path): The directory path for environmental data.

        Raises:
            ValueError: If the provided value is not a valid directory path.
        """
        if not value.is_dir():
            raise ValueError(f"{value} is not a valid directory path.")
        
        self._enviro_dir = value

    @property
    def clim_data_file(self):
        """str: The file name of the climate data file."""
        return self._clim_data_file
    
    @clim_data_file.setter
    def clim_data_file(self, value):
        """Sets the file name of the climate data file.

        Args:
            value (str): The file name of the climate data file.

        Raises:
            TypeError: If the provided value is not a string.
            ValueError: If the provided file does not exist in the climate data directory.
        """
        # Validate type
        if not isinstance(value, str):
            raise TypeError("clim_data_file attribute nust be a string.")
        
        # Validate path and file
        valid_file = self.climate_dir / value
        if not valid_file.exists():
            raise ValueError(f"{valid_file} does not exists")
        
        else:
            self._clim_data_file = value

    @property
    def clim_meta_file(self):
        """str: The file name of the climate metadata file."""
        return self._clim_meta_file
    
    @clim_meta_file.setter
    def clim_meta_file(self, value):
        """Sets the file name of the climate metadata file.

        Args:
            value (str): The file name of the climate metadata file.

        Raises:
            TypeError: If the provided value is not a string.
            ValueError: If the provided file does not exist in the climate data directory.
        """
        # Validate type
        if not isinstance(value, str):
            raise TypeError("clim_meta_file attribute nust be a string.")
        
        # Validate path and file
        valid_file = self.climate_dir / value
        if not valid_file.exists():
            raise ValueError(f"{valid_file} does not exists")
        
        else:
            self._clim_meta_file = value
    
    @property
    def env_data_file(self):
        """str: The file name of the environmental data file."""
        return self._env_data_file
    
    @env_data_file.setter
    def env_data_file(self, value):
        """Sets the file name of the environmental data file.

        Args:
            value (str): The file name of the environmental data file.

        Raises:
            TypeError: If the provided value is not a string.
            ValueError: If the provided file does not exist in the environmental data directory.
        """
        # Validate type
        if not isinstance(value, str):
            raise TypeError("env_data_file attribute nust be a string.")
        
        # Validate path and file
        valid_file = self.enviro_dir / value
        if not valid_file.exists():
            raise ValueError(f"{valid_file} does not exists")
        
        else:
            self._env_data_file = value

    @property
    def env_meta_file(self):
        """str: The file name of the environmental metadata file."""
        return self._env_meta_file
    
    @env_meta_file.setter
    def env_meta_file(self, value):
        """Sets the file name of the environmental metadata file.

        Args:
            value (str): The file name of the environmental metadata file.

        Raises:
            TypeError: If the provided value is not a string.
            ValueError: If the provided file does not exist in the environmental data directory.
        """
        # Validate type
        if not isinstance(value, str):
            raise TypeError("env_meta_file attribute nust be a string.")
        
        # Validate path and file
        valid_file = self.enviro_dir / value
        if not valid_file.exists():
            raise ValueError(f"{valid_file} does not exists")
        
        else:
            self._env_meta_file = value
        

if __name__ == "__main__":
    
    test = MadaclimLayers()

    print(test)
    