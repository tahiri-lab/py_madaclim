import json
import re
import pathlib
import requests
import time
import inspect
from calendar import month_name
from pathlib import Path
from typing import List, Union, Optional, Tuple, Dict
import importlib.resources as pkg_resources

import pandas as pd
import geopandas as gpd
import numpy as np
import rasterio
import pyproj
from shapely.geometry import Point

from coffeaphylogeo._constants import Constants

class MadaclimLayers:
    """A class that represents all of the information and data from the climate and environmental variable layers that can be found from the rasters of the Madaclim database.

    The main metadata retrieval tool for the Madaclim database. Access all layers information with the `all_layers` attribute.
    Also provides methods to filter, generate unique labels from all_layers and also access the crs and band number from the climate and environmental rasters.
    
    Attributes:
        clim_raster (pathlib.Path): The path to the Madaclim climate raster GeoTiff file. Defaults to None if not specified.
        env_raster (pathlib.Path): The path to the Madaclim environmental raster GeoTif file. Defaults to None if not specified.
        all_layers (pd.DataFrame): A DataFrame containing a complete and formatted version of all Madaclim layers.
    
    """
    def __init__(self, clim_raster: Optional[pathlib.Path]=None, env_raster: Optional[pathlib.Path]=None):
        """Initializes a new instance of the MadaclimLayers class.

        This constructor sets the directory paths and file names for climate and environment data,
        and generates a DataFrame containing all Madaclim layers.
        
        Examples:
            >>> from coffeaphylogeo.madaclim_info import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> # Access the all layers df
            >>> all_layers_df = madaclim_info.all_layers
            >>> all_layers_df.info()
            <class 'pandas.core.frame.DataFrame'>
            Int64Index: 79 entries, 0 to 8
            Data columns (total 5 columns):
            #   Column             Non-Null Count  Dtype 
            ---  ------             --------------  ----- 
            0   layer_number       79 non-null     int64 
            1   geoclim_feature    79 non-null     object
            2   geoclim_type       79 non-null     object
            3   layer_name         79 non-null     object
            4   layer_description  71 non-null     object
            dtypes: int64(1), object(4)

            >>> # 'clim_raster' and 'env_raster' attributes are empty by default
            >>> print(madaclim_info.clim_raster)
            None
            >>> print(madaclim_info.env_raster)
            None
            
            >>> # Certain attributes and methods need a valid 'clim_raster' or 'env_raster' attribute a-priori
            >>> madaclim_info.clim_crs
            Traceback (most recent call last):
            ...
                raise ValueError(f"Undefined attribute: '{raster_attr_name}'. You need to assign a valid pathlib.Path to the related raster attribute first.")
            ValueError: Undefined attribute: 'clim_raster'. You need to assign a valid pathlib.Path to the related raster attribute first.

            >>> # You can download the rasters using the 'download_data' method
            >>> madaclim_info.download_data()    # Defaults to current working dir otherwise specify save_dir pathlib.Path
            
            >>> madaclim_info.clim_raster = Path("madaclim_current.tif")
            >>> madaclim_info.env_raster = Path("madaclim_enviro.tif")
            >>> madaclim_info.clim_raster
            PosixPath('madaclim_current.tif')
            >>> madaclim_info.clim_crs
            <Derived Projected CRS: EPSG:32738>
            Name: WGS 84 / UTM zone 38S
            Axis Info [cartesian]:
            - E[east]: Easting (metre)
            - N[north]: Northing (metre)
            Area of Use:
            - name: Between 42°E and 48°E, southern hemisphere between 80°S and equator, onshore and offshore. Madagascar.
            - bounds: (42.0, -80.0, 48.0, 0.0)
            Coordinate Operation:
            - name: UTM zone 38S
            - method: Transverse Mercator
            Datum: World Geodetic System 1984 ensemble
            - Ellipsoid: WGS 84
            - Prime Meridian: Greenwich
        """
        self._clim_dataformat = self._load_dataformat(Constants.CLIM_DATAFORMAT_FILE)
        self._clim_metadata = self._load_metadata(Constants.CLIM_METADATA_FILE)

        self._env_dataformat = self._load_dataformat(Constants.ENV_DATAFORMAT_FILE)
        self._env_metadata = self._load_metadata(Constants.ENV_METADATA_FILE)
        
        self.all_layers = self._get_madaclim_layers()

        self.clim_raster = clim_raster
        self.env_raster = env_raster

    

    @property
    def clim_raster(self) -> pathlib.Path:
        """pathlib.Path: Get or set the path to the climate raster file.

        This property allows you to get the current path to the climate raster file, or set a new
        path. If setting a new path, the value must be a pathlib.Path object or a str. If the value 
        is a str, it will be converted to a pathlib.Path object. The path must exist, otherwise a 
        FileNotFoundError will be raised.

        Returns:
            The current path to the climate raster file.

        Raises:
            TypeError: If the new path is not a pathlib.Path object or str.
            ValueError: If the new path cannot be converted to a pathlib.Path object.
            FileNotFoundError: If the new path does not exist.
        """
        return self._clim_raster
    
    @clim_raster.setter
    def clim_raster(self, value: Optional[pathlib.Path]) -> None:
        if value is None:
            self._clim_raster = value
        else:
            # Validate type
            if not isinstance(value, (pathlib.Path, str)):
                raise TypeError("clim_raster must be a pathlib.Path object or str.")
            
            # Validate path and file
            try:
                value = Path(value)
            except:
                raise ValueError(f"Could not create a pathlib.Path object from {value}")
            
            if not value.exists():
                raise FileNotFoundError(f"{value} does not exists.")
            
            self._clim_raster = value

    @property
    
    def env_raster(self) -> pathlib.Path:
        """pathlib.Path: Get or set the path to the environment raster file.

        This property allows you to get the current path to the environment raster file, or set a new
        path. If setting a new path, the value must be a pathlib.Path object or a str. If the value 
        is a str, it will be converted to a pathlib.Path object. The path must exist, otherwise a 
        FileNotFoundError will be raised.

        Returns:
            The current path to the environment raster file.

        Raises:
            TypeError: If the new path is not a pathlib.Path object or str.
            ValueError: If the new path cannot be converted to a pathlib.Path object.
            FileNotFoundError: If the new path does not exist.
        """
        return self._env_raster
    
    @env_raster.setter
    def env_raster(self, value: Optional[pathlib.Path]) -> None:
        if value is None:
            self._env_raster = value
        else:       
            # Validate type
            if not isinstance(value, (pathlib.Path, str)):
                raise TypeError("env_raster must be a pathlib.Path object or str.")
            
            # Validate path and file
            try:
                value = Path(value)
            except:
                raise ValueError(f"Could not create a pathlib.Path object from {value}")
            
            if not value.exists():
                raise FileNotFoundError(f"{value} does not exists.")
            
            self._env_raster = value

    @property
    def clim_crs(self) -> pyproj.crs.CRS:
        """
        Retrieves the Coordinate Reference System (CRS) from the Madaclim climate raster.

        This property first validates the `clim_raster` attribute, ensuring its integrity and existence. 
        It then opens the raster file and retrieves the CRS in EPSG format. The EPSG code is used to 
        create and return a pyproj CRS object.

        Returns:
            pyproj.crs.CRS: The CRS object derived from the EPSG code of the climate raster.

        Examples:
            >>> # You need to have a valid 'clim_raster' attribute before
            >>> madaclim_info = MadaclimLayers()
            >>> madaclim_info.clim_crs
            Traceback (most recent call last):
            ...
                raise ValueError(f"Undefined attribute: '{raster_attr_name}'. You need to assign a valid pathlib.Path to the related raster attribute first.")
            ValueError: Undefined attribute: 'clim_raster'. You need to assign a valid pathlib.Path to the related raster attribute first.

            >>> madaclim_info.clim_raster = Path("madaclim_current.tif")
            >>> madaclim_info.clim_crs
            <Derived Projected CRS: EPSG:32738>
            Name: WGS 84 / UTM zone 38S
            Axis Info [cartesian]:
            - E[east]: Easting (metre)
            - N[north]: Northing (metre)
            Area of Use:
            - name: Between 42°E and 48°E, southern hemisphere between 80°S and equator, onshore and offshore. Madagascar.
            - bounds: (42.0, -80.0, 48.0, 0.0)
            Coordinate Operation:
            - name: UTM zone 38S
            - method: Transverse Mercator
            Datum: World Geodetic System 1984 ensemble
            - Ellipsoid: WGS 84
            - Prime Meridian: Greenwich
        """
        
        # Validate raster attr value, IO path, integrity
        self._validate_raster("clim_raster")
        
        # Get epsg from clim_raster
        with rasterio.open(self.clim_raster) as clim_raster:
            clim_epsg = clim_raster.crs.to_epsg()  # Get the EPSG code of the CRS
            clim_crs = pyproj.CRS.from_epsg(clim_epsg)  # Create a pyproj CRS object
        return clim_crs

    @property
    def env_crs(self):
        """
        Retrieves the Coordinate Reference System (CRS) from the Madaclim environmental raster.

        This property first validates the `env_raster` attribute, ensuring its integrity and existence. 
        It then opens the raster file and retrieves the CRS in EPSG format. The EPSG code is used to 
        create and return a pyproj CRS object.

        Returns:
            pyproj.crs.CRS: The CRS object derived from the EPSG code of the environmental raster.
        
        Examples:
            >>> # You need to have a valid 'env_raster' attribute
            >>> madaclim_info = MadaclimLayers()
            >>> madaclim_info.env_crs
            Traceback (most recent call last):
            ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ...                
            ValueError: Undefined attribute: 'env_raster'. You need to assign a valid pathlib.Path to the related raster attribute first.
            
            >>> madaclim_info.env_raster = Path("madaclim_enviro.tif")
            >>> madaclim_info.env_crs
            <Derived Projected CRS: EPSG:32738>
            Name: WGS 84 / UTM zone 38S
            Axis Info [cartesian]:
            - E[east]: Easting (metre)
            - N[north]: Northing (metre)
            Area of Use:
            - name: Between 42°E and 48°E, southern hemisphere between 80°S and equator, onshore and offshore. Madagascar.
            - bounds: (42.0, -80.0, 48.0, 0.0)
            Coordinate Operation:
            - name: UTM zone 38S
            - method: Transverse Mercator
            Datum: World Geodetic System 1984 ensemble
            - Ellipsoid: WGS 84
            - Prime Meridian: Greenwich
        """
        
        # Validate raster attr value, IO path, integrity
        self._validate_raster("env_raster")
        
        # Get epsg from env_raster
        with rasterio.open(self.env_raster) as env_raster:
            env_epsg = env_raster.crs.to_epsg()  # Get the EPSG code of the CRS
            env_crs = pyproj.CRS.from_epsg(env_epsg)  # Create a pyproj CRS object
        return env_crs
    
    def __str__(self) -> str:
        """Prints the `MadaclimLayers` instance's attributes.

        Returns:
            str: All the object's attributes as attr_name for keys and attr_value for values.
        """
        info = (
            f"all_layers = \n{self.all_layers.head().to_string(index=False)}\n...\n"
            f"{len(self.all_layers)} rows x {len(self.all_layers.columns)} columns\n\n"  
            f"clim_raster = {self.clim_raster}\n"
            f"clim_crs = {self.clim_crs if self.clim_raster else None}\n"
            f"env_raster = {self.env_raster}\n"
            f"env_crs = {self.env_crs if self.env_raster else None}\n"   
        )
        info = "MadaclimLayers\n" + info
        return info


    def select_geoclim_type_layers(self, geoclim_type: str) -> pd.DataFrame:
        """Method that selects the desired geoclimatic type layers as a dataframe.

        Args:
            geoclim_type (str): The desired geoclimatic layers type to extract.

        Returns:
            pd.DataFrame: A slice of the `all_layers` dataframe containing the desired geoclimatic type layers.

        Raises:
            TypeError: If geoclim_type is not a string.
            ValueError: If geoclim_type does not corresponds to a valid geoclim type.
        
        Examples:
            >>> from coffeaphylogeo.madaclim_layers import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> clim_df = madaclim_info.select_geoclim_type_layers(geoclim_type="clim")
            >>> clim_df.head()
            layer_number                        geoclim_feature geoclim_type layer_name                                 layer_description
            0             1  Monthly minimum temperature (°C x 10)         clim      tmin1   Monthly minimum temperature (°C x 10) - January
            1             2  Monthly minimum temperature (°C x 10)         clim      tmin2  Monthly minimum temperature (°C x 10) - February
            2             3  Monthly minimum temperature (°C x 10)         clim      tmin3     Monthly minimum temperature (°C x 10) - March
            3             4  Monthly minimum temperature (°C x 10)         clim      tmin4     Monthly minimum temperature (°C x 10) - April
            4             5  Monthly minimum temperature (°C x 10)         clim      tmin5       Monthly minimum temperature (°C x 10) - May
        
        """
        all_layers_df = self.all_layers.copy()
        # Validate geoclim_type
        if not isinstance(geoclim_type, str):
            raise TypeError("geoclim_type must be a string.")
        
        possible_geoclim_types = all_layers_df["geoclim_type"].unique() 
        if geoclim_type not in possible_geoclim_types:
            raise ValueError(f"geoclim_type must be one of {possible_geoclim_types}")
        
        select_df = all_layers_df[all_layers_df["geoclim_type"] == geoclim_type]

        return select_df
    
    def get_layers_labels(self, layers_subset: Optional[Union[str, List[int]]]=None, as_descriptive_labels: bool=False) -> list:
        """
        Retrieves unique layer labels based on the provided subset of layers.

        This method fetches the unique labels from the `all_layers` dataframe, 
        given a subset of layers (specified as either layer numbers, a geoclim_type, or a single layer number).
        The layer labels can be returned in a descriptive format if `as_descriptive_labels` is set to True.

        Args:
            layers_subset (Optional[Union[str, List[int]]], optional): A list of layer numbers or a geoclim_type string 
                to subset the labels from, or a single layer number as a string or int. Defaults to None, which will 
                select all layers (no subset).
            as_descriptive_labels (bool, optional): If True, returns the descriptive layer labels. Otherwise, returns 
                the "layer_<num>" format. Defaults to False.

        Raises:
            TypeError: If elements of `layers_subset` cannot be converted to int.
            ValueError: If `layers_subset` is a string not in `possible_geoclim_types`, cannot be converted to int, 
                or if the 'layer_number' and 'layer_name' columns in the `all_layers` dataframe have non-unique entries.

        Returns:
            list: A list of unique layer labels. These labels are either in the "layer_<num>" format or the descriptive format, 
            based on `as_descriptive_labels`.
        
        Example:
            >>> # Get labels for all layers
            >>> from coffeaphylogeo.madaclim_layers import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> all_layers = madaclim_info.get_layers_labels()
            
            >>> # Specify a geoclim subset
            >>> env_layers = madaclim_info.get_layers_labels(layers_subset="env")
            >>> env_layers
            ['layer_71', 'layer_72', 'layer_73', 'layer_74', 'layer_75', 'layer_76', 'layer_77', 'layer_78', 'layer_79']
            
            >>> # Extract more information (can also be the input for 'fetch_specific_layers' method)
            >>> informative_labels = madaclim_info.get_layers_labels(as_descriptive_labels=True)
            >>> informative_labels[:5]
            ['clim_1_tmin1 (Monthly minimum temperature (°C x 10) - January)', 'clim_2_tmin2 (Monthly minimum temperature (°C x 10) - February)', 'clim_3_tmin3 (Monthly minimum temperature (°C x 10) - March)', 'clim_4_tmin4 (Monthly minimum temperature (°C x 10) - April)', 'clim_5_tmin5 (Monthly minimum temperature (°C x 10) - May)']

            >>> # Specify a single layer or a subset of layers
            >>> madaclim_info.get_layers_labels(37, as_descriptive_labels=True)
            ['clim_37_bio1 (Annual mean temperature)']
            >>> madaclim_info.get_layers_labels([68, 75], as_descriptive_labels=True)
            ['clim_68_pet (Annual potential evapotranspiration from the Thornthwaite equation (mm))', 'env_75_geology (1=Alluvial_&_Lake_deposits, 2=Unconsolidated_Sands, 4=Mangrove_Swamp, 5=Tertiary_Limestones_+_Marls_&_Chalks, 6=Sandstones, 7=Mesozoic_Limestones_+_Marls_(inc._"Tsingy"), 9=Lavas_(including_Basalts_&_Gabbros), 10=Basement_Rocks_(Ign_&_Met), 11=Ultrabasics, 12=Quartzites, 13=Marble_(Cipolin))']

            >>> # Example to get bioclim layers only
            >>> bioclim_labels = [label for label in madaclim_info.get_layers_labels(as_descriptive_labels=True) if "bio" in label]
        """
        all_layers_df = self.all_layers.copy()
        layers_numbers = all_layers_df["layer_number"].to_list()
        possible_geoclim_types = all_layers_df["geoclim_type"].unique()


        # Convert layers_subset to a list of ints for all inputs
        
        if isinstance(layers_subset, list):
            try:
                layers_subset = [int(layer) for layer in layers_subset]
            except:
                raise TypeError("'layers_subet' list elements must be int or can be converted to int.")
            
        elif isinstance(layers_subset, (str, int)):
            # Extract layer num if layers_subset is a geoclim_type
            if layers_subset in possible_geoclim_types:
                layers_subset = all_layers_df.loc[all_layers_df["geoclim_type"] == layers_subset, "layer_number"].to_list()

            else:
                try:
                    layers_subset = [int(layers_subset)]
                except ValueError:
                    raise ValueError(f"'layers_subset' must be an int (or can be converted to int) or be one of {possible_geoclim_types}")
        
        elif layers_subset is None:
            layers_subset = layers_numbers

        else:
            raise TypeError("'layers_subset' must be a str, int or a list of str or ints.")

        # Validate layer number range for layer_numbers as in
        min_layer, max_layer = min(layers_numbers), max(layers_numbers)

        for layer_number in layers_subset:
            if not min_layer <= layer_number <= max_layer:
                raise ValueError(f"layer_number must fall between {min_layer} and {max_layer}. {layer_number=} is not valid.")
            
        # Validate unique entries
        if len(all_layers_df) != len(all_layers_df["layer_number"].unique()) != len(all_layers_df["layer_name"].unique()):
            raise ValueError("'layer_number' and 'layer_name' columns in the all_layers dataframe have non-unique entries.")
        
         # Fetch rows according to layers_subset
        select_df = all_layers_df[all_layers_df["layer_number"].isin(layers_subset)]
        
        # Fetch layer_<num> and descriptive labels
        sub_selection = list(zip(select_df["layer_number"], select_df["geoclim_type"], select_df["layer_name"], select_df["layer_description"]))
        layers_description = {
            f"layer_{num}": f"{geoclim}_{num}_{name} ({', '.join(desc) if type(desc) == list else desc})" for num, geoclim, name, desc in sub_selection
        }
        
        if as_descriptive_labels:
            unique_labels = list(layers_description.values())
        else:
            unique_labels = list(layers_description.keys())
        
        return unique_labels
        

    def fetch_specific_layers(self, layers_labels: Union[int, str, List[Union[int, str]]], as_descriptive_labels: bool=False, return_list: bool=False) -> Union[dict, pd.DataFrame, list]:
        """
        Fetches specific layers from the `all_layers` DataFrame based on the given input.

        Args:
            layers_labels (Union[int, str, List[Union[int, str]]]): The layer labels to fetch. Can be a single int or str value, or a list of int or str values.
                The input can also be in the format "layer_{num}" or "{geotype}_{num}_{name}_({description})" (output from `get_layers_labels(as_descriptive_labels=True)` method).
            as_descriptive_labels (bool): If True, only the layer descriptions are returned. Defaults to False.

        Returns:
            Union[dict, pd.DataFrame]: If `as_descriptive_labels` is True, returns a dictionary with the layer descriptions.
                Otherwise, returns a DataFrame with the specified layers.

        Raises:
            TypeError: If any value in layers_labels cannot be converted to an int or is not in the "layer_{num}" format.
            ValueError: If any layer_number does not fall between the minimum and maximum layer numbers.

        Example:
            >>> from coffeaphylogeo.madaclim_layers import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> madaclim_info.fetch_specific_layers([1, 15, 55, 71])    # Output is a pd.DataFrame
                layer_number                        geoclim_feature geoclim_type layer_name                                layer_description
            0              1  Monthly minimum temperature (°C x 10)         clim      tmin1  Monthly minimum temperature (°C x 10) - January
            14            15  Monthly maximum temperature (°C x 10)         clim      tmax3    Monthly maximum temperature (°C x 10) - March
            54            55        Bioclimatic variables (bioclim)         clim      bio19                 Precipitation of coldest quarter
            0             71                           Altitude (m)          env   altitude                                             None
            
            >>> # Using the output from get_layers_labels() method
            >>> bioclim_labels = [label for label in madaclim_info.get_layers_labels(as_descriptive_labels=True) if "bio" in label]
            >>> bio1_to_bio5_labels = bioclim_labels[0:6]
            >>> madaclim_info.fetch_specific_layers(bio1_to_bio5_labels)
                layer_number                  geoclim_feature geoclim_type layer_name                                  layer_description
            36            37  Bioclimatic variables (bioclim)         clim       bio1                            Annual mean temperature
            37            38  Bioclimatic variables (bioclim)         clim       bio2  Mean diurnal range (mean of monthly (max temp ...
            38            39  Bioclimatic variables (bioclim)         clim       bio3                  Isothermality (BIO2/BIO7) (x 100)
            39            40  Bioclimatic variables (bioclim)         clim       bio4  Temperature seasonality (standard deviation x ...
            40            41  Bioclimatic variables (bioclim)         clim       bio5                   Max temperature of warmest month
            
            >>> # Or from descriptive_labels as well
            >>> monthly_layers = [layer for layer in madaclim_info.get_layers_labels(as_descriptive_labels=True) if "Monthly" in layer]
            >>> monthly_layers[:3]
            ['clim_1_tmin1 (Monthly minimum temperature (°C x 10) - January)', 'clim_2_tmin2 (Monthly minimum temperature (°C x 10) - February)', 'clim_3_tmin3 (Monthly minimum temperature (°C x 10) - March)']
            >>> madaclim_info.fetch_specific_layers(monthly)[:5]
            layer_number  ...                                 layer_description
            0             1  ...   Monthly minimum temperature (°C x 10) - January
            1             2  ...  Monthly minimum temperature (°C x 10) - February
            2             3  ...     Monthly minimum temperature (°C x 10) - March
            3             4  ...     Monthly minimum temperature (°C x 10) - April
            4             5  ...       Monthly minimum temperature (°C x 10) - May

            
            >>> # Fetch description only with output as dict (instead of pd.DataFrame)
            >>> madaclim_info.fetch_specific_layers(bio1_to_bio5_labels, as_descriptive_labels=True)
            {'layer_37': 'clim_37_bio1 (Annual mean temperature)', 'layer_38': 'clim_38_bio2 (Mean diurnal range (mean of monthly (max temp - min temp)))', 'layer_39': 'clim_39_bio3 (Isothermality (BIO2/BIO7) (x 100))', 'layer_40': 'clim_40_bio4 (Temperature seasonality (standard deviation x 100))', 'layer_41': 'clim_41_bio5 (Max temperature of warmest month)'}

        """
        all_layers_df = self.all_layers.copy()    # Reference to all clim and env metadata df

        # Validate layers_labels
        possible_layers_num_format = [f"layer_{num}" for num in all_layers_df["layer_number"].to_list()]
        possible_layers_desc_format = self.get_layers_labels(as_descriptive_labels=True)

        if isinstance(layers_labels, list):
            # Check if all elements are in unique label format
            layers_num_format = all([layer_label in possible_layers_num_format for layer_label in layers_labels])
            layers_desc_format = all([layer_label in possible_layers_desc_format for layer_label in layers_labels])
            
            # Save as list of ints for later filtering
            if layers_num_format or layers_desc_format:
                layers_numbers = [int(layer_label.split("_")[1]) for layer_label in layers_labels]

            # layers_labels as list of ints
            else:
                try:
                    layers_numbers = [int(layer) for layer in layers_labels]
                except (ValueError, TypeError):
                    raise TypeError("layers_labels must be either a list of int values (or str that can be converted to int) or the output format from the 'get_layers_labels' method.")
        
        # Single layers_labels type check 
        else:
            if layers_labels in possible_layers_num_format:    # Check layer_<num> str format
                layers_numbers = [int(layers_labels.split("_")[1])]
            try:
                layers_numbers = [int(layers_labels)]
            except (ValueError, TypeError):
                raise TypeError("layers_labels must be either a single int value (or a str that can be converted to an int) or in the output format from the 'get_layers_labels' method.")
  
        # Validate layer number range for layer_numbers as in
        min_layer = min(all_layers_df["layer_number"])
        max_layer = max(all_layers_df["layer_number"])

        for layer_number in layers_numbers:
            if not min_layer <= layer_number <= max_layer:
                raise ValueError(f"layer_number must fall between {min_layer} and {max_layer}. {layer_number=} is not valid.")
            
        # Fetch rows according to layer selection
        select_df = all_layers_df[all_layers_df["layer_number"].isin(layers_numbers)]

        if as_descriptive_labels:
            # Generate dict with key as layer_<num> and values containing layer information
            sub_selection = list(zip(select_df["layer_number"], select_df["geoclim_type"], select_df["layer_name"], select_df["layer_description"]))
            layers_description = {
                f"layer_{num}": f"{geoclim}_{num}_{name} ({', '.join(desc) if type(desc) == list else desc})" for num, geoclim, name, desc in sub_selection
            }
            
            if return_list:
                return list(layers_description.values())
            else:
                return layers_description
        
        else:    # Save whole row of df
            return select_df

    def download_data(self, save_dir: Optional[pathlib.Path]=None):
        """Downloads climate and environment raster files from the Madaclim website.

        This method downloads the climate and environment raster data from the Madaclim website
        and saves them to the specified directory. If no directory is specified, the data is saved
        to the current working directory.

        Args:
            save_dir (Optional[pathlib.Path]): The directory where the data should be saved. If not specified,
                the data is saved to the current working directory.

        Raises:
            ValueError: If save_dir is not a directory.

        Examples:
            >>> madaclim_info.download_data()    # Defaults to current working directory

            ####   Trying get request to Madaclim website...   ####
            madaclim_current.tif is 21.8 MB
            Server response OK from madaclim.cirad.fr, starting to download madaclim_current.tif
            Progress for madaclim_current.tif : 100.00 % completed of 21.8 MB downloaded [ average speed of  3.3 MB/s ]
            Done downloading madaclim_current.tif in 6.65 seconds !

            ####   Trying get request to Madaclim website...   ####
            madaclim_enviro.tif is 5.5 MB
            Server response OK from madaclim.cirad.fr, starting to download madaclim_enviro.tif
            Progress for madaclim_enviro.tif : 100.00 % completed of 5.5 MB downloaded [ average speed of  2.7 MB/s ]
            Done downloading madaclim_enviro.tif in 2.07 seconds !
        """
        
        def download_single_file(url: str, dir_savepath: pathlib.Path, filename: str)-> pathlib.Path:
            """Downloads a single file from a URL and saves it to the specified directory.

            This function downloads a file from the specified URL and saves it to the specified directory
            with the specified filename. The download progress is printed to the console.

            Args:
                url (str): The URL of the file to download.
                dir_savepath (pathlib.Path): The directory where the file should be saved.
                filename (str): The name of the file to save.

            Returns:
                pathlib.Path: The path of the downloaded file.

            """
            try :
                response = requests.get(url, stream=True)
                print("\n####   Trying get request to Madaclim website...   ####")
                
                # Calculate file size in MB
                total_size = (float(response.headers['Content-Length']))/1000000
                print(f"{filename} is {total_size:.1f} MB")
                
                # Download in chunks of 1 MB and save to disk
                if response.status_code == 200 : 
                    print(f"Server response OK from {url.split('/')[2]}, starting to download {filename}")
                    with open(dir_savepath / filename, 'wb') as f :
                        chunksize = 1024 * 1000
                        start_time = time.time()
                        for n, chunk in enumerate(response.iter_content(chunk_size=chunksize)) :
                            percent = (n * chunksize / (total_size*1000000)) * 100
                            now_time = time.time()
                            current_speed = (n * chunksize) / (now_time - start_time) / 1000000
                            print(
                                f"Progress for {filename} : {percent:.2f} % completed of {total_size:.1f} MB downloaded [ current speed of  {current_speed:.1f} MB/s ]",
                                end="\r"
                            )
                            time.sleep(0.1)
                            f.write(chunk)
                        end_time = time.time()
                        print(
                            f"Progress for {filename} : 100.00 % completed of {total_size:.1f} MB downloaded [ average speed of  {total_size/(end_time-start_time):.1f} MB/s ]",
                            end="\r"
                        )
                        print()
                    print(f"Done downloading {filename} in {end_time-start_time:.2f} seconds !")
                            
            except :
                print(f"File {filename} cannot be downloaded. Status code : {response.status_code}")
            return dir_savepath / filename
        
        # Validate save target directory
        save_dir = Path.cwd() if save_dir is None else save_dir
        if not save_dir.is_dir():
            raise FileNotFoundError(f"Cannot find directory: {save_dir}")

        # Download rasters
        download_single_file(    # Climate raster
            url=Constants.MADACLIM_URLS["clim_raster"],
            dir_savepath=save_dir,
            filename=Constants.DEFAULT_CLIM_RASTER_FILENAME
        )

        download_single_file(
            url=Constants.MADACLIM_URLS["env_raster"],
            dir_savepath=save_dir,
            filename=Constants.DEFAULT_ENV_RASTER_FILENAME
        )

    def get_bandnums_from_layers(self, layers_labels : Union[int, str, List[Union[int, str]]]) -> List[int]:
        """
        Retrieves band numbers corresponding to the provided layers' labels.

        This method accepts labels for a subset of layers (specified as either layer numbers, "layer_<num>" format, or descriptive labels)
        and returns the corresponding band numbers from the `all_layers` dataframe. If the input is in the descriptive label format or "layer_<num>" format,
        it should match the output of the `get_layers_labels` method.

        Args:
            layers_labels (Union[int, str, List[Union[int, str]]]): A list of layer labels in various formats, or a single layer label.

        Raises:
            TypeError: If elements of `layers_labels` cannot be converted to int or if they do not match the format produced by the `get_layers_labels` method.
            ValueError: If the derived layer numbers do not fall within the valid range of layer numbers in the `all_layers` dataframe.

        Returns:
            List[int]: A list of band numbers corresponding to the provided layer labels.

        Example:
            >>> from coffeaphylogeo.madaclim_layers import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> last_20 = madaclim_info.get_layers_labels()[-20:]
            >>> band_nums = madaclim_info.get_bandnums_from_layers(last_20)
            >>> band_nums
            [60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        all_layers_df = self.all_layers.copy()    # Reference to all clim and env metadata df

        # Validate layers_labels
        possible_layers_num_format = [f"layer_{num}" for num in all_layers_df["layer_number"].to_list()]
        possible_layers_desc_format = self.get_layers_labels(as_descriptive_labels=True)

        if isinstance(layers_labels, list):
            # Check if all elements are in unique label format
            layers_num_format = all([layer_label in possible_layers_num_format for layer_label in layers_labels])
            layers_desc_format = all([layer_label in possible_layers_desc_format for layer_label in layers_labels])
            
            # Save as list of ints for later filtering
            if layers_num_format or layers_desc_format:
                layers_numbers = [int(layer_label.split("_")[1]) for layer_label in layers_labels]

            # layers_labels as list of ints
            else:
                try:
                    layers_numbers = [int(layer) for layer in layers_labels]
                except (ValueError, TypeError):
                    raise TypeError("layers_labels must be either a list of int values (or str that can be converted to int) or the output format from the 'get_layers_labels' method.")
        
        # Single layers_labels type check 
        else:
            if layers_labels in possible_layers_num_format:    # Check layer_<num> str format
                layers_numbers = [int(layers_labels.split("_")[1])]
            try:
                layers_numbers = [int(layers_labels)]
            except (ValueError, TypeError):
                raise TypeError("layers_labels must be either a single int value (or a str that can be converted to an int) or in the output format from the 'get_layers_labels' method.")
  
        # Validate layer number range for layer_numbers as in
        min_layer = min(all_layers_df["layer_number"])
        max_layer = max(all_layers_df["layer_number"])

        for layer_number in layers_numbers:
            if not min_layer <= layer_number <= max_layer:
                raise ValueError(f"layer_number must fall between {min_layer} and {max_layer}. {layer_number=} is not valid.")
            
        # Fetch rows and geoclim type(s) according to selected layers
        select_df = all_layers_df[all_layers_df["layer_number"].isin(layers_numbers)]
        geoclim_types = select_df["geoclim_type"].unique()
        
        # Get band num according to geoclim_type
        total_clim_layers = (all_layers_df["geoclim_type"] == "clim").sum()
        band_nums = []

        for geoclim_type in geoclim_types:
            self._validate_raster(f"{geoclim_type}_raster")    # Validate raster attr val, IO path and integrity
            if geoclim_type == "clim":
                band_nums += select_df.loc[select_df["geoclim_type"] == geoclim_type, "layer_number"].to_list()
            else:
                band_nums += (select_df.loc[select_df["geoclim_type"] == geoclim_type, "layer_number"] - total_clim_layers).to_list()
        
        return band_nums

    #! Deprecated: Unefficient method with raster.read I/O operation for each band
    def _get_band_from_layer_number(self, layer_number: Union[str, int], geoclim_type: str)->int:
        """Get the band number in a raster file corresponding to a given layer number and geoclim type.

        Args:
            layer_number (Union[str, int]): The layer number to find the corresponding band number for.
            geoclim_type (str): The geoclim type, either "clim" or "env".

        Raises:
            TypeError: If geoclim_type is not a string.
            ValueError: If geoclim_type is not one of the possible geoclim types ("clim" or "env").
            ValueError: If layer_number is out of range for the selected geoclim raster.
            ValueError: If the band number could not be retrieved for the given layer number.

        Returns:
            int: The band number in the raster file corresponding to the given layer number and geoclim type.
        """
        
        # Validate geoclim_type
        if not isinstance(geoclim_type, str):
            raise TypeError("geoclim_type must be a string.")
        
        possible_geoclim_types = ["clim", "env"]
        if geoclim_type not in possible_geoclim_types:
            raise ValueError(f"geoclim_type must be one of {possible_geoclim_types}")
        
        # Validate layer_number according to geoclim type
        select_geoclim_layers_range = self.select_geoclim_type_layers(geoclim_type)["layer_number"].to_list()

        if layer_number not in select_geoclim_layers_range:
            raise ValueError(f"{layer_number=} is out of range for the selected geoclim raster of {geoclim_type}. Choose between {min(select_geoclim_layers_range)} and {max(select_geoclim_layers_range)}")

        # Get the number of bands for the selected geoclim_type raster
        if geoclim_type == "clim":
            self._validate_raster(f"{geoclim_type}_raster")    # Validate raster before I/O
            with rasterio.open(self.clim_raster) as clim_raster:
                raster_bands = clim_raster.read().shape[0]
        
        if geoclim_type == "env":
            self._validate_raster(f"{geoclim_type}_raster")    # Validate raster before I/O
            with rasterio.open(self.env_raster) as env_raster:
                raster_bands = env_raster.read().shape[0]

        # Return the band number according to the layer number
        band_number = layer_number - select_geoclim_layers_range[0] + 1
        
        if band_number not in range(1, raster_bands + 1):
            raise ValueError(f"Could not retrieve {band_number=} for {layer_number=}")
        
        return band_number

    def _get_madaclim_layers(self) -> pd.DataFrame :
        """Private method that will generate the `all_layers` attributes based on the format and metada files from the Madaclim db by accessing their corresponding attributes.

        Args:
            climate_dir (Path): The directory path for climate-related data.
            enviro_dir (Path): The directory path for environment-related data.
            clim_dataformat_filename (str): The file name of the climate data format file.
            clim_meta_filename (str): The file name of the climate metadata file.
            env_dataformat_filename (str): The file name of the environment data format file.
            env_meta_filename (str): The file name of the environment metadata file.

        Returns:
            pd.DataFrame: A DataFrame containing a complete and formatted version of all Madaclim layers.
        """
        
        def split_layers(row, layers_col_name: str):
            """Split "Layers" column and create new rows

            Args:
                row (pd.Series): Row of data from a DataFrame
                layers_col_name (str): Name of the column to split

            Returns:
                pd.DataFrame: DataFrame with split rows

            Examples:
                >>> df = pd.DataFrame({"A": ["x1-3", "y4"], "B": ["foo", "bar"]})
                >>> result = df.apply(split_repeating_vars, axis=1, col_to_split="A", col_to_keep="B")
                >>> df_result = pd.concat(result.tolist(), axis=0).reset_index(drop=True)
                >>> print(df_result)
                    A    B
                0   x1  foo
                1   x2  foo
                2   x3  foo
                3   y4  bar
            """
            if "-" in row[layers_col_name]:
                start, end = map(int, row[layers_col_name].split("-"))
                index = range(start, end+1)
            else:
                index = [int(row[layers_col_name])]
            
            data = {
                "Climate variable": [row["Climate variable"]] * len(index),
                "data_type": [row["data_type"]] * len(index)
            }
            
            return pd.DataFrame(data, index=index)
        
        def split_repeating_vars(row: pd.Series, col_to_split: str, col_to_keep: str) -> pd.DataFrame:
            """Split a column containing repeating values into separate rows.

            This function takes a row of data from a DataFrame and splits the value in the column specified by `col_to_split` if it contains a hyphen. The function returns a new DataFrame with rows for each value in the specified range and with values from the column specified by `col_to_keep`. The month name is appended to the `col_to_keep` values.

            Args:
                row (pd.Series): A row of data from a DataFrame.
                col_to_split (str): The name of the column to split.
                col_to_keep (str): The name of the column to keep.

            Returns:
                pd.DataFrame: A DataFrame with split rows.
            """
            
            # Extract the range and layername to a new smaller df of len(range(start, end))
            if "-" in row[col_to_split]:
                start = int(re.search("\d+", row[col_to_split].split("-")[0]).group())
                end = int(row[col_to_split].split("-")[1])
                name = re.search("[a-z]*", row[col_to_split]).group()
                        
                # Create a DataFrame with the split values and the description column
                df = pd.DataFrame({col_to_split: [f"{name}{month}" for month in range(start, end+1)],
                                col_to_keep: [f"{row[col_to_keep]} - {month_name[month]}" for month in range(start, end+1)]})
                return df
            
            # When no changes to changes to row
            else:
                df = pd.DataFrame({col_to_split: [row[col_to_split]], col_to_keep: [row[col_to_keep]]})
                return df
            
        def add_layer_numbers_bio_monthly(df: pd.DataFrame) -> pd.DataFrame:
            """Adds layer numbers to bio_monthly dataframe.

            This function takes the bio_monthly ddf and adds a new column
            "layer_number" that assigns a layer number to each row based on the value
            of the "layer_name" column.

            Args:
                df (pd.DataFrame): The bio_monthly dataframe.

            Returns:
                pd.DataFrame: The input dataframe with an additional "layer_number" column.
            """
            
            # Get the geoclim feature name and sub df for each layer_name categories
            monthly_features = {}
            categories_layer_name = df["layer_name"].str.extract("(^[a-zA-Z]+)")[0].unique()

            # Add layer number to bio_monthly df
            current_start_layer = 1
            for layer_name in categories_layer_name:
                
                # Extract common base layer name
                category_feature_val = df[df["layer_name"].str.contains(layer_name)]["layer_description"].str.split(" - ").str[0].unique()[0]
                monthly_features[f"{layer_name}_cat_feature"] = category_feature_val
                
                # Get df associated with common category
                category_feature_df = pd.DataFrame(df[df["layer_name"].str.contains(layer_name)])
                monthly_features[f"{layer_name}_df"] = category_feature_df

                # Assign layer number according to current category
                df.loc[df["layer_name"].str.contains(layer_name), "layer_number"] = range(current_start_layer, current_start_layer + len(category_feature_df))
                current_start_layer += len(category_feature_df)

            df["layer_number"] = df["layer_number"].astype(int)

            return df
        
        def meta_merge_clim_df(clim_df: pd.DataFrame, meta_dfs: List[pd.DataFrame]) -> pd.DataFrame:
            """Merges a the original clim_df dataframe with multiple metadata dataframes.

            Args:
                clim_df (pd.DataFrame): A climate dataframe with a "layer_number" column.
                meta_dfs (List[pd.DataFrame]): A list of metadata dataframes with a "layer_number" column.

            Returns:
                pd.DataFrame: The concatenated merged dataframe.
            """
            
            merge_result_dfs = []
            for meta_df in meta_dfs:
                result_df = pd.merge(clim_df, meta_df, on="layer_number")
                merge_result_dfs.append(result_df)
            
            merged_df = pd.concat(merge_result_dfs).reset_index(drop=True)

            return merged_df
        
        
        # Extract climate data and format it using the metadata
        df_clim = pd.read_json(self._clim_dataformat["table_0"])
        df_clim["data_type"] = "clim"    # Tag for latter id in merge

        # Split the layers to get initial layer_num for all clim-related layers
        df_clim = pd.concat((df_clim.apply(split_layers, args=("Layers", ), axis=1)).to_list()).reset_index()
        df_clim.columns = ["layer_number", "geoclim_feature", "geoclim_type"]

        # Formatting + add layers to the monthly bioclim metadata
        bio_monthly_feats = pd.read_json(self._clim_metadata["table_0"])
        bio_monthly_feats.columns = ["layer_name", "layer_description"]
        
        bio_monthly_feats = pd.concat(    # Split according to range
            bio_monthly_feats.apply(split_repeating_vars, axis=1, args=("layer_name", "layer_description", )).to_list(),
            ignore_index=True
        )
        bio_monthly_feats = add_layer_numbers_bio_monthly(bio_monthly_feats)    # Append layer numbers

        # Formatting + add layers to the other bioclim metadata (non-monthly)
        bioclim_feats = pd.read_json(self._clim_metadata["table_1"])
        bioclim_feats.columns = ["layer_name", "layer_description"]

        current_start_layer = len(bio_monthly_feats) + 1    # Save the current state of the layer number for clim_df
        bioclim_feats["layer_number"] = range(current_start_layer, current_start_layer + len(bioclim_feats))

        # Formatting + add layers to monthly and annual evapotranspiration metadata
        evap_feats = pd.read_json(self._clim_metadata["table_2"])
        evap_feats.columns = ["layer_name", "layer_description"]
        
        evap_feats = pd.concat(    # Split monthly evapo data
            evap_feats.apply(split_repeating_vars, axis=1, args=("layer_name", "layer_description", )).to_list(),
            ignore_index=True
        )
        current_start_layer = max(bioclim_feats["layer_number"]) + 1    # Save the current state of the layer number for clim_df
        evap_feats["layer_number"] = range(current_start_layer, current_start_layer + len(evap_feats))

        # Formatting + add layers to the bioclim water-related metadata
        biowater_feats = pd.read_json(self._clim_metadata["table_3"])
        biowater_feats.columns = ["layer_name", "layer_description"]
        
        current_start_layer = max(evap_feats["layer_number"]) + 1    # Save the current state of the layer number for clim_df
        biowater_feats["layer_number"] = range(current_start_layer, current_start_layer + len(biowater_feats))

        # Merge meta_dfs with original clim_df for class attribute
        meta_dfs = [bio_monthly_feats, bioclim_feats, evap_feats, biowater_feats]
        df_clim = meta_merge_clim_df(df_clim, meta_dfs)


        # Extract environmental data and format it using its related metadata
        df_env = pd.read_json(self._env_dataformat["table_0"])
        df_env.columns = ["layer_number", "geoclim_feature"]
        df_env["geoclim_type"] = "env"    # Tag for latter id in merge

        current_start_layer = max(df_clim["layer_number"]) + 1    # Save the current state of the layer number according to the final clim_df
        df_env["layer_number"] = range(current_start_layer, current_start_layer + len(df_env))

        # Generate layer_name since absent from metadata
        df_env["layer_name"] = df_env["geoclim_feature"].str.split(" ").str[0].str.lower()
        df_env.loc[df_env["layer_number"] == 79, "layer_name"] = "forestcover"    # Fix first word with more informative info
        df_env["layer_description"] = None

        # Assign dummy var information for geology layer to layer_description
        geology_description = []
    
        env_meta_str = self._env_metadata["table_0"]
        env_meta_data = json.loads(env_meta_str)

        for i, val in env_meta_data["Raster value"].items():
            rock_type = env_meta_data["Rock type"][i]
            rock_type = "_".join(rock_type.split(" "))
            rock_type_categorical = f"{val}={rock_type}"
            geology_description.append(rock_type_categorical)

        df_env.at[(df_env["layer_name"] == "geology").idxmax(), "layer_description"] = geology_description
        
        # Concat both clim and env final dfs
        df = pd.concat([df_clim, df_env])
        df = df.reset_index().drop(columns="index")

        return df
    
    def _load_dataformat(self, filepath: pathlib.Path) -> dict:
        """Parse the json dataformat (either clim or env) file into a dictionary.

        Args:
            filepath (pathlib.Path): Path to the dataformat file.

        Raises:
            ValueError: If filepath is not a path to an existing file.
            ValueError: If json file cannot be parsed.

        Returns:
            dict: A dictionary containing the data format information contained in the dataformat file.
        """
        
        try:
            with open(filepath) as f:
                return json.load(f)
        except FileNotFoundError:
            raise ValueError(f".json dataformat file not found: {filepath}")
        except json.JSONDecodeError:
            raise ValueError(f"Failed to parse json dataformat file: {filepath}")
        
    def _load_metadata(self, filepath: pathlib.Path) -> dict:
        """Parse the json metadata (either clim or env) file into a dictionary.

        Args:
            filepath (pathlib.Path): Path to the metadata file.

        Raises:
            ValueError: If filepath is not a path to an existing file.
            ValueError: If json file cannot be parsed.

        Returns:
            dict: A dictionary containing the metadata contained in the metadata file.
        """
        
        try:
            with open(filepath) as f:
                return json.load(f)
        except FileNotFoundError:
            raise ValueError(f".json dataformat file not found: {filepath}")
        except json.JSONDecodeError:
            raise ValueError(f"Failed to parse json dataformat file: {filepath}")
    
    def _validate_raster(self, raster_attr_name: str):
        """
        Validates the specified raster attribute.

        This method checks whether the specified raster attribute exists, whether the corresponding raster file exists, 
        and whether the raster file can be opened without any IO errors.

        Args:
            raster_attr_name (str): Name of the raster attribute to be validated.

        Raises:
            ValueError: If the specified raster attribute is not defined.
            FileExistsError: If the corresponding raster file does not exist.
            IOError: If the raster file cannot be opened due to IO errors.
        """    
        raster_to_check = getattr(self, raster_attr_name)    # Fetch current value of raster attr

        if raster_to_check is None:
            raise ValueError(f"Undefined attribute: '{raster_attr_name}'. You need to assign a valid pathlib.Path to the related raster attribute first.")

        # Check if raster file exists
        if not raster_to_check.is_file():
            raise FileExistsError(f"Could not find '{raster_attr_name}' file: {raster_to_check}")
                
        # Catch any IO errors
        try:
            with rasterio.open(raster_to_check) as raster_file:
                return
        except rasterio.errors.RasterioIOError as e:
            raise IOError(f"Could not open '{raster_attr_name}': {raster_to_check}. Error: {e}")
