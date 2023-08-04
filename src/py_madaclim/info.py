import json
import re
import pathlib
import requests
import time
from calendar import month_name
from pathlib import Path
from typing import List, Union, Optional, Dict
from itertools import zip_longest

import pandas as pd
import numpy as np
import rasterio
import pyproj

from py_madaclim._constants import Constants

class MadaclimLayers:
    """A class that represents all of the information and data from the climate and environmental variable layers 
        that can be found from the rasters of the Madaclim database.

    The main metadata retrieval tool for the Madaclim database. Access all layers information with the `all_layers` attribute.
        Also provides methods to filter, generate unique labels from the `all_layers` attr.
        Access the crs and band number from the climate and environmental rasters when they are provided in the constructor.
        Categorical data can be explored in details with the `categorical_layers` attribute and the value:category pairs with the `get_categorical_combinations`.
    
    Attributes:
        clim_raster (pathlib.Path): The path to the Madaclim climate raster GeoTiff file. Defaults to None if not specified.
        env_raster (pathlib.Path): The path to the Madaclim environmental raster GeoTif file. Defaults to None if not specified.
        all_layers (pd.DataFrame): A DataFrame containing a complete and formatted version of all Madaclim layers.
        categorical_layers (pd.DataFrame):  A DataFrame containing the in depth information about the layers with categorical data.
    
    """
    def __init__(self, clim_raster: Optional[pathlib.Path]=None, env_raster: Optional[pathlib.Path]=None):
        """Initializes a new instance of the MadaclimLayers class.

        This constructor generates a DataFrame containing all Madaclim layers' information though the `all_layers` attribute.
            If given the current climate or environmental raster from the Madaclim db, it will also get the CRS for each raster.
            The instance can also access various methods to extract relevant information from the Madaclim database.

        
        Example:
            >>> from py_madaclim.info import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> # Access the all layers df
            >>> all_layers_df = madaclim_info.all_layers
            >>> all_layers_df.shape
            (79, 6)
            >>> all_layers_df
            geoclim_type  layer_number layer_name                                  layer_description  is_categorical                                              units
            0          clim             1      tmin1              Monthly minimum temperature - January           False                                            °C x 10
            1          clim             2      tmin2             Monthly minimum temperature - February           False                                            °C x 10
            2          clim             3      tmin3                Monthly minimum temperature - March           False                                            °C x 10
            3          clim             4      tmin4                Monthly minimum temperature - April           False                                            °C x 10
            4          clim             5      tmin5                  Monthly minimum temperature - May           False                                            °C x 10
            ..          ...           ...        ...                                                ...             ...                                                ...
            74          env            75        geo                                         Rock types            True  [1=Alluvial_&_Lake_deposits, 2=Unconsolidated_...
            75          env            76        soi                                         Soil types            True  [1=Bare_Rocks, 2=Raw_Lithic_Mineral_Soils, 3=P...
            76          env            77        veg                                   Vegetation types            True  [1=VegCat_1, 2=VegCat_2, 3=VegCat_3, 4=VegCat_...
            77          env            78        wat                                         Watersheds            True  [1=N-Bemarivo, 2=S-Bemarivo,_N-Mangoro, 3=S-Ma...
            78          env            79     forcov  Percentage of forest cover in 1 km by 1 km gri...           False                                                  %

            [79 rows x 6 columns]

            >>> # Categorical layers only df
            >>> madaclim_info.categorical_layers
            geoclim_type  layer_number layer_name layer_description value                              category
            0           env            75        geo        Rock types     1              Alluvial_&_Lake_deposits
            1           env            75        geo        Rock types     2                  Unconsolidated_Sands
            2           env            75        geo        Rock types     4                        Mangrove_Swamp
            3           env            75        geo        Rock types     5  Tertiary_Limestones_+_Marls_&_Chalks
            4           env            75        geo        Rock types     6                            Sandstones
            ..          ...           ...        ...               ...   ...                                   ...
            73          env            78        wat        Watersheds    20                  ret-disp_Tsiribihina
            74          env            78        wat        Watersheds    21                    ret-disp_Betsiboka
            75          env            78        wat        Watersheds    22            ret-disp_Maevarana_(1/2_N)
            76          env            78        wat        Watersheds    23                    ret-disp_Sambirano
            77          env            78        wat        Watersheds    24                     ret-disp_Mahavavy

            [78 rows x 6 columns]


            >>> # 'clim_raster' and 'env_raster' attributes are empty by default
            >>> print(madaclim_info.clim_raster)
            None
            >>> print(madaclim_info.env_raster)
            None
            
            >>> # Certain attributes and methods need a valid 'clim_raster' or 'env_raster' attribute a-priori
            >>> madaclim_info.clim_crs
            Traceback (most recent call last):
            File "<stdin>", line 1, in <module>
            File "/home/local/USHERBROOKE/lals2906/programming/python_projects/coffeaphylogeo/src/py_madaclim/madaclim_info.py", line 223, in clim_crs
                self._validate_raster("clim_raster")
            File "/home/local/USHERBROOKE/lals2906/programming/python_projects/coffeaphylogeo/src/py_madaclim/madaclim_info.py", line 1111, in _validate_raster
                raise ValueError(f"Undefined attribute: '{raster_attr_name}'. You need to assign a valid pathlib.Path to the related raster attribute first.")
            ValueError: Undefined attribute: 'clim_raster'. You need to assign a valid pathlib.Path to the related raster attribute first.
            
            >>> # You can download the rasters using the 'download_data' method
            >>> madaclim_info.download_data()    # Defaults to current working dir otherwise specify save_dir pathlib.Path
        
            >>> madaclim_info.clim_raster = "madaclim_current.tif"
            >>> madaclim_info.env_raster = "madaclim_enviro.tif"
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
        self.clim_raster = clim_raster
        self.env_raster = env_raster

        self._clim_dataformat = self._load_dataformat(Constants.CLIM_DATAFORMAT_FILE)
        self._clim_metadata = self._load_metadata(Constants.CLIM_METADATA_FILE)

        self._env_dataformat = self._load_dataformat(Constants.ENV_DATAFORMAT_FILE)
        self._env_metadata = self._load_metadata(Constants.ENV_METADATA_FILE)
        
        self._all_layers = self._get_madaclim_layers()
        self._categorical_layers = self._get_categorical_df()

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
    def clim_crs(self) -> pyproj.crs.crs.CRS:
        """
        Retrieves the Coordinate Reference System (CRS) from the Madaclim climate raster.

        This property first validates the `clim_raster` attribute, ensuring its integrity and existence. 
            It then opens the raster file and retrieves the CRS in EPSG format. The EPSG code is used to 
            create and return a pyproj CRS object.

        Returns:
            pyproj.crs.crs.CRS: The CRS object derived from the EPSG code of the climate raster.

        Example:
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
    def env_crs(self) -> pyproj.crs.crs.CRS:
        """
        Retrieves the Coordinate Reference System (CRS) from the Madaclim environmental raster.

        This property first validates the `env_raster` attribute, ensuring its integrity and existence. 
        It then opens the raster file and retrieves the CRS in EPSG format. The EPSG code is used to 
        create and return a pyproj CRS object.

        Returns:
            pyproj.crs.crs.CRS: The CRS object derived from the EPSG code of the environmental raster.
        
        Example:
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
    
    @property
    def all_layers(self) -> pd.DataFrame:
        """
        Retrieves the 'all_layers' Dataframe using the private '_get_madaclim_layers' method.

        Contains all information about all the raster layers in the Madaclim db.

        Returns:
            pd.DataFrame: A DataFrame containing a complete and formatted version of all Madaclim layers.
        """
        return self._all_layers
    
    @property
    def categorical_layers(self) -> pd.DataFrame:
        """
        Retrieves the 'categorical_layers' Dataframe using the private '_get_categorical_df' method.

        Contains detailed information about the categorical layers from the rasters in the Madaclim db.

        Returns:
            pd.DataFrame: A DataFrame containing information for each categorical value in each layer
        """
        return self._categorical_layers
    
    def __str__(self) -> str:
        """Prints the `MadaclimLayers` instance's attributes.

        Returns:
            str: All the object's attributes as attr_name for keys and attr_value for values.
        """
        categ_layer_nums = self._categorical_layers['layer_number'].unique().astype(str)
        
        custom_public_methods = []
        for attr in dir(self):
            attr_passed = False
            try:
                attr_val = getattr(self, attr)
                attr_passed = True
            except:
                attr_val = None
            if (attr_passed and
                callable(attr_val) and not
                attr.startswith("_")):
                custom_public_methods.append(attr)
        # Remove raster specific methods if not present in instance
        if not self.clim_raster and not self.env_raster:
            raster_specif_pub_meth = ["get_bandnums_from_layers"] 
            custom_public_methods = [
                method for method in custom_public_methods if method not in raster_specif_pub_meth
            ]
        
        clim_raster_info = (
            f"\tclim_raster = {self.clim_raster.name if self.clim_raster else None}\n"
            f"\tclim_crs = {self.clim_crs if self.clim_raster else None}\n"
        ) if self.clim_raster else ""
        env_raster_info = (
            f"\tenv_raster = {self.env_raster.name if self.env_raster else None}\n"
            f"\tenv_crs = {self.env_crs if self.env_raster else None}\n"
        ) if self.env_raster else ""

        base_info = (
            
            f"\tall_layers = {type(self._all_layers).__name__}({self._all_layers.shape[0]} rows x {self._all_layers.shape[1]} columns)\n"
            f"\tcategorical_layers = {type(self._categorical_layers).__name__}"
            f"(Layers {', '.join(categ_layer_nums)} "
            f"with a total of {len(self._categorical_layers)} categories\n"
            f"{clim_raster_info}"
            f"{env_raster_info}"
            f"\tpublic methods -> {', '.join(custom_public_methods[:3])}\n"
            f"\t\t\t {', '.join(custom_public_methods[3:])}\n"
        )

        full_info = "MadaclimLayers(\n" + base_info + ")"
        return full_info
    
    def __repr__(self) -> str:
        return self.__str__()

    def select_geoclim_type_layers(self, geoclim_type: str) -> pd.DataFrame:
        """Method that selects the desired geoclimatic type layers as a dataframe.

        Args:
            geoclim_type (str): The desired geoclimatic layers type to extract.

        Returns:
            pd.DataFrame: A slice of the `all_layers` dataframe containing the desired geoclimatic type layers.

        Raises:
            TypeError: If geoclim_type is not a string.
            ValueError: If geoclim_type does not corresponds to a valid geoclim type.
        
        Example:
            >>> from py_madaclim.madaclim_layers import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> clim_df = madaclim_info.select_geoclim_type_layers(geoclim_type="clim")
            >>> clim_df.head()
            geoclim_type  layer_number layer_name                       layer_description  is_categorical    units
            0         clim             1      tmin1   Monthly minimum temperature - January           False  °C x 10
            1         clim             2      tmin2  Monthly minimum temperature - February           False  °C x 10
            2         clim             3      tmin3     Monthly minimum temperature - March           False  °C x 10
            3         clim             4      tmin4     Monthly minimum temperature - April           False  °C x 10
            4         clim             5      tmin5       Monthly minimum temperature - May           False  °C x 10

        """
        all_layers_df = self._all_layers.copy()
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
            >>> from py_madaclim.madaclim_layers import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> all_layers = madaclim_info.get_layers_labels()
            >>> len(all_layers)
            79
            >>> # Basic format 'layer_<num>'
            >>> all_layers[:5]
            ['layer_1', 'layer_2', 'layer_3', 'layer_4', 'layer_5']

            >>> # Specify a geoclim subset
            >>> env_layers = madaclim_info.get_layers_labels(layers_subset="env")
            >>> env_layers
            ['layer_71', 'layer_72', 'layer_73', 'layer_74', 'layer_75', 'layer_76', 'layer_77', 'layer_78', 'layer_79']
            
            >>> # Extract more information
            >>> # Format is a list of 'type_num_uniqname_description (units)' elements
            >>> informative_labels = madaclim_info.get_layers_labels(as_descriptive_labels=True)
            >>> informative_labels[:5]
            >>> informative_labels[:2]
            ['clim_1_tmin1_Monthly minimum temperature - January (°C x 10)', 'clim_2_tmin2_Monthly minimum temperature - February (°C x 10)']
            >>> # Any output from `get_layers_labels` can be used as input for other methods in other classes such as `fetch_specific_layers` from `MadaclimLayers`

            >>> # Specify a single layer or a subset of layers
            >>> madaclim_info.get_layers_labels(37, as_descriptive_labels=True)
            ['clim_37_bio1_Annual mean temperature (degrees)']
            >>> madaclim_info.get_layers_labels([68, 75], as_descriptive_labels=True)
            ['clim_68_pet_Annual potential evapotranspiration from the Thornthwaite equation (mm)', 'env_75_geo_Rock types (categ_vals: 1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13)']

            >>> # Example to get bioclim layers only
            >>> bioclim_labels = [label for label in madaclim_info.get_layers_labels(as_descriptive_labels=True) if "bio" in label]
        """
        all_layers_df = self._all_layers.copy()
        layers_numbers = all_layers_df["layer_number"].to_list()
        possible_geoclim_types = all_layers_df["geoclim_type"].unique()


        # Convert layers_subset to a list of ints for all inputs
        
        if isinstance(layers_subset, list):
            try:
                layers_subset = [int(layer) for layer in layers_subset]
            except:
                raise TypeError("'layers_subet' list elements must be int or can be converted to int.")
            
        elif isinstance(layers_subset, (str, int, np.integer)):
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
        sub_selection = list(
            zip(
                select_df["layer_number"], 
                select_df["geoclim_type"], 
                select_df["layer_name"], 
                select_df["layer_description"], 
                select_df["units"]
            )
        )
        layers_description = {}
        for num, geoclim, name, desc, unit in sub_selection:
            # Simplify units from 'all_layers' for categorical
            if isinstance(unit, list):
                unit = "categ_vals: " + ", ".join([ele.split('=')[0] for ele in unit])
            layers_description[f"layer_{num}"] = f"{geoclim}_{num}_{name}_{desc} ({unit})"
        
        if as_descriptive_labels:
            unique_labels = list(layers_description.values())
        else:
            unique_labels = list(layers_description.keys())
        
        return unique_labels
        

    def fetch_specific_layers(self, layers_labels: Union[int, str, List[Union[int, str]]], *args: str) -> Union[dict, pd.DataFrame]:
        """
        Fetches specific layers from the `all_layers` DataFrame based on the given input and returns either the entire rows or 
        certain columns as a dictionary. 

        Args:
            layers_labels (Union[int, str, List[Union[int, str]]]): The layer labels to fetch. Can be a single int or str value,
                or a list of int or str values. The input can also be in the format "layer_{num}" or 
                "{geotype}_{num}_{name}_({description})" (output from `get_layers_labels(as_descriptive_labels=True)` method).
            *args (str): Optional. One or more column names in `all_layers` DataFrame. If specified, only these columns
            will be returned as a dictionary.

        Returns:
            Union[dict, pd.DataFrame]: If `args` is specified, returns a nested dictionary with the format:
                {
                    layer_<num>: {
                        <arg1>: <value>,
                        <arg2>: <value>,
                        ...
                    },
                    ...
                }
                Otherwise, returns a DataFrame with the specified layers.

        Raises:
            TypeError: If any value in layers_labels cannot be converted to an int or is not in the "layer_{num}" format.
            ValueError: If any layer_number does not fall between the minimum and maximum layer numbers.
            KeyError: If any value in args is not a column in `all_layers` DataFrame.

        Example:
            >>> from py_madaclim.madaclim_layers import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> madaclim_info.fetch_specific_layers([1, 15, 55, 71])
            geoclim_type  layer_number layer_name                      layer_description  is_categorical         units
            0          clim             1      tmin1  Monthly minimum temperature - January           False       °C x 10
            14         clim            15      tmax3    Monthly maximum temperature - March           False       °C x 10
            54         clim            55      bio19       Precipitation of coldest quarter           False  mm.3months-1
            70          env            71        alt                               Altitude           False        meters

            >>> # Using the output from `get_layers_labels` method
            >>> bioclim_labels = [label for label in madaclim_info.get_layers_labels(as_descriptive_labels=True) if "bio" in label]
            >>> bio1_to_bio5_labels = bioclim_labels[0:6]
            >>> madaclim_info.fetch_specific_layers(bio1_to_bio5_labels)
            geoclim_type  layer_number layer_name                 layer_description  is_categorical                                       units
            36         clim            37       bio1           Annual mean temperature           False                                     degrees
            37         clim            38       bio2                Mean diurnal range           False  mean of monthly max temp - monthy min temp
            38         clim            39       bio3   Isothermality = BIO2/BIO7 x 100           False                                    No units
            39         clim            40       bio4           Temperature seasonality           False                    standard deviation x 100
            40         clim            41       bio5  Max temperature of warmest month           False                                     degrees
            41         clim            42       bio6  Min temperature of coldest month           False                                     degrees

            >>> # Or from descriptive_labels as well
            >>> len(pet_layers)
            13
            >>> pet_layers[-1]
            'clim_68_pet_Annual potential evapotranspiration from the Thornthwaite equation (mm)'
            >>> madaclim_info.fetch_specific_layers(pet_layers)
            geoclim_type  layer_number layer_name                                  layer_description  is_categorical       units
            55         clim            56       pet1  Monthly potential evapotranspiration from the ...           False  mm.month-1
            56         clim            57       pet2  Monthly potential evapotranspiration from the ...           False  mm.month-1
            57         clim            58       pet3  Monthly potential evapotranspiration from the ...           False  mm.month-1
            58         clim            59       pet4  Monthly potential evapotranspiration from the ...           False  mm.month-1
            59         clim            60       pet5  Monthly potential evapotranspiration from the ...           False  mm.month-1
            60         clim            61       pet6  Monthly potential evapotranspiration from the ...           False  mm.month-1
            61         clim            62       pet7  Monthly potential evapotranspiration from the ...           False  mm.month-1
            62         clim            63       pet8  Monthly potential evapotranspiration from the ...           False  mm.month-1
            63         clim            64       pet9  Monthly potential evapotranspiration from the ...           False  mm.month-1
            64         clim            65      pet10  Monthly potential evapotranspiration from the ...           False  mm.month-1
            65         clim            66      pet11  Monthly potential evapotranspiration from the ...           False  mm.month-1
            66         clim            67      pet12  Monthly potential evapotranspiration from the ...           False  mm.month-1
            67         clim            68        pet  Annual potential evapotranspiration from the T...           False          mm

            >>> # Fetch as dict with keys as layer_<num> and vals of choice using 
            >>> madaclim_info.fetch_specific_layers([15, 55, 75], "geoclim_type", "layer_name", "is_categorical")
            {
                'layer_15': {
                    'geoclim_type': 'clim',
                    'layer_name': 'tmax3',
                    'is_categorical': False
                },
                'layer_55': {
                    'geoclim_type': 'clim',
                    'layer_name': 'bio19',
                    'is_categorical': False
                },
                'layer_75': {
                    'geoclim_type': 'env',
                    'layer_name': 'geo',
                    'is_categorical': True}
                }
            }
            >>> # Only col names will be accepted as additionnal args
            >>> bio1 = next((layer for layer in mada_info.get_layers_labels(as_descriptive_labels=True) if "bio1" in layer), None)
            >>> madaclim_info.fetch_specific_layers(bio1, "band_number")
            Traceback (most recent call last):
            File "<stdin>", line 1, in <module>
            File "/home/local/USHERBROOKE/lals2906/programming/python_projects/py_madaclim/src/py_madaclim/madaclim_info.py", line 604, in fetch_specific_layers
                if not min_layer <= layer_number <= max_layer:
                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            KeyError: "Invalid args: ['band_number']. Args must be one of a key of ['geoclim_type', 'layer_number', 'layer_name', 'layer_description', 'is_categorical', 'units'] or 'all'"
            >>> # Get all keys with the `all` argument
            >>> madaclim_info.fetch_specific_layers(bio1, "all")
            {
                'layer_37': {
                    'geoclim_type': 'clim',
                    'layer_number': 37,
                    'layer_name': 'bio1',
                    'layer_description': 'Annual mean temperature',
                    'is_categorical': False,
                    'units': 'degrees'
                }
            }

        """
        all_layers_df = self._all_layers.copy()    # Reference to all clim and env metadata df

        # Validate layers_labels
        possible_layers_num_format = self.get_layers_labels()
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
            if isinstance(layers_labels, str):
                if layers_labels in possible_layers_num_format or layers_labels in possible_layers_desc_format:    # Check layer_<num> or desc str format
                    layers_numbers = [int(layers_labels.split("_")[1])]
                else:
                    raise TypeError("layers_labels must be in the output format from the 'get_layers_labels' method if it is a string.")
            
            elif isinstance(layers_labels, (int, np.integer)):
                try:
                    layers_numbers = [layers_labels]
                except (ValueError, TypeError):
                    raise TypeError("layers_labels must be an int if not in string format.")
            
            else:
                raise TypeError("layers_labels must be either a single int value, a string that can be converted to an int, or in the output format from the 'get_layers_labels' method.")

        # Validate layer number range for layer_numbers as in
        min_layer = min(all_layers_df["layer_number"])
        max_layer = max(all_layers_df["layer_number"])

        for layer_number in layers_numbers:
            if not min_layer <= layer_number <= max_layer:
                raise ValueError(f"layer_number must fall between {min_layer} and {max_layer}. {layer_number=} is not valid.")
            
        # Fetch rows according to layer selection
        select_df = all_layers_df[all_layers_df["layer_number"].isin(layers_numbers)]

        # Validate args for presence in `all_layers` attr or for `all` possible keys
        if len(args) > 0:
            missing_args = list(set(args) - set(list(select_df.columns) + ["all"]))
            if missing_args:
                raise KeyError(f"Invalid args: {missing_args}. Args must be one of a key of {list(select_df.columns)} or 'all'")
            if "all" in args and len(args) > 1:
                raise ValueError("Cannot have additional arguments when 'all' is specified")
            
            cols = select_df.columns if "all" in args else args
            # Return nested dicts with values of specified arg for the fetched layers
            if not missing_args:
                select_layers_dict = {}
                for layer_num in layers_numbers:
                    layer_label = f"layer_{layer_num}"
                    select_layers_dict[layer_label] = {}
                    for col in cols:
                        select_layers_dict[layer_label][col] = select_df[select_df["layer_number"] == layer_num][col].values[0]    
                return select_layers_dict
        
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

        Example:
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
            >>> from py_madaclim.madaclim_layers import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> last_20 = madaclim_info.get_layers_labels()[-20:]
            >>> band_nums = madaclim_info.get_bandnums_from_layers(last_20)
            >>> band_nums
            [60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        all_layers_df = self._all_layers.copy()    # Reference to all clim and env metadata df

        # Validate layers_labels
        possible_layers_num_format = self.get_layers_labels()
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
    
    def get_categorical_combinations(
            self, 
            layers_labels: Optional[Union[int, str, List[Union[int, str]]]]=None, 
            as_descriptive_keys: bool=False
        ) -> Union[dict, Dict[str, Dict[int, str]]]:
        """
        Returns a dictionary representation of the specified categorical layers corresponding the the categorical value encoding.

        Args:
            layers_labels (Optional[Union[int, str, List[Union[int, str]]]]): The layer labels to fetch. Can be a single 
            integer or string value, or a list of integer or string values. The input can also be in the format 
            "layer_{num}" or "{geotype}_{num}_{name}_({unit})" (output from `get_layers_labels(as_descriptive_labels=True)` method).
            If `layers_labels` is `None`, all categorical layers are fetched.
            as_descriptive_keys(bool)

        Raises:
            TypeError: If `layers_labels` is not a list of integers or strings, a single integer or a string 
            that can be converted to an integer, or in the output format from the 'get_layers_labels' method.
            ValueError: If a layer number in `layers_labels` is not a valid categorical layer number.

        Returns:
            Union[dict, Dict[str, Dict[int, str]]]: A dictionary of the specified categorical layers. 
            If multiple layers were specified, the dictionary keys are 'layer_{num}', and the values are dictionaries 
            with layer values as keys and their corresponding categories as values.
            If a single layer was specified, the dictionary keys are the categorical values, and the values are the 
            categories themselves. 

        Example:
            >>> # If multiple layers specified, it returns:
            >>> madaclim_info = MadaclimLayers()
            >>> >>> madaclim_info.get_categorical_combinations([75, 76])
            {
                'layer_75': {
                    1: 'N-Bemarivo',
                    2: 'S-Bemarivo,_N-Mangoro',
                    ...
                },
                'layer_76': {
                    1: 'Bare_Rocks',
                    2: 'Raw_Lithic_Mineral_Soils',
                    ...
                },
                ...
            }
            >>> # If a single layer is specified, it returns:
            >>> madaclim_info.get_categorical_combinations("layer_76")
            {
                'layer_76: {
                    1: 'Bare_Rocks',
                    2: 'Raw_Lithic_Mineral_Soils',
                ...
                }
            }
            >>> # For more descriptive keys (same output from as_descriptive_labels)
            >>> madaclim_info.get_categorical_combinations("layer_76", as_descriptive_keys=True)
            {
                'env_76_soi_Soil types (categ_vals: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)': {
                    1: 'Bare_Rocks',
                    2: 'Raw_Lithic_Mineral_Soils',
                ...
                }
            }
        """
        cat_df = self._categorical_layers.copy()
        
        # Validate layers_labels
        possible_layers_num_format = self.get_layers_labels()
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
        elif isinstance(layers_labels, str):
            if layers_labels in possible_layers_num_format or layers_labels in possible_layers_desc_format:    # Check layer_<num> or desc str format
                layers_numbers = [int(layers_labels.split("_")[1])]
            else:
                raise TypeError("layers_labels must be in the output format from the 'get_layers_labels' method if it is a string.")
        
        elif isinstance(layers_labels, (int, np.integer)):
            try:
                layers_numbers = [layers_labels]
            except (ValueError, TypeError):
                raise TypeError("layers_labels must be an int if not in string format.")
            
        elif layers_labels is None:
            layers_numbers = list(cat_df["layer_number"].unique())

        else:
            raise ValueError("layers_labels must be either a single int value, a string that can be converted to an int, or in the output format from the 'get_layers_labels' method.")
            
        # Validate layer number range for layer_numbers as in
        all_cat_layers_num = cat_df["layer_number"].unique()

        for layer_number in layers_numbers:
            if layer_number not in all_cat_layers_num:
                raise ValueError(f"layer_number must be one of the categorical layers: {all_cat_layers_num}. {layer_number=} is not valid.")
            
        select_cat_df = cat_df.loc[cat_df["layer_number"].isin(layers_numbers)]    # Fetch validated layers

        categorical_dict = {}    # Nested container dict
        for layer_number in layers_numbers:
            if as_descriptive_keys:
                categorical_key = self.get_layers_labels(layer_number, as_descriptive_labels=True)[0]
            else:
                categorical_key = self.get_layers_labels(layer_number)[0]
            categorical_dict[categorical_key] = {
                int(row["raster_value"]): row["category"] 
                for _, row in select_cat_df.loc[select_cat_df["layer_number"] == layer_number].iterrows()
            }
        return categorical_dict
    
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
        
        def get_env_layer_name(row: pd.Series) -> str:
            """
            Generate an environment layer name based on the `geoclim_feature` field in the given row.
            
            Args:
                row (pd.Series): A row of the DataFrame.
            
            Returns:
                str: The generated environment layer name.
            """
            sourceless_unitless = row["geoclim_feature"].split(" (")[0]
            sless_uless_list = sourceless_unitless.split(" ")
            if len(sless_uless_list) > 2:
                name = f"{sless_uless_list[2][:3]}" + f"{sless_uless_list[3][:3]}"
            elif len(sless_uless_list) > 1:
                name = f"{sless_uless_list[0][:3]}" + f"{sless_uless_list[1][:3]}"
            else:
                name = sourceless_unitless[:3]
            return name.lower()
        
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
        
        # Define your function to extract units
        def extract_units(row: pd.Series) -> Union[str, list]:
            """
            Extract units from the `layer_description` field if the `units` field is None.
            Otherwise, return the value in the `units` field.
            
            Args:
                row (pd.Series): A row of the DataFrame.
                
            Returns:
                Union[str, list]: The extracted units or the original value in the `units` field.
            """
            if row["units"] is None: 
                match = re.search(r'\((.*?)\)', row["layer_description"])
                return match.group(1) if match else "No units"
            else:
                return row["units"] 

            
        def remove_units_from_desc(row: pd.Series) -> str:
            """
            Remove units from the `layer_description` field if `units` field is not a list.
            Otherwise, return the original `layer_description`.
            
            Args:
                row (pd.Series): A row of the DataFrame.
                
            Returns:
                str: The processed `layer_description` without units if `units` is not a list, 
                    otherwise the original `layer_description`.
            """
            if isinstance(row["units"], list):  # Apply only non-categorical (null) entries
                return row["layer_description"]
            else:
                return re.sub(r' \(.*?\)', "", row["layer_description"])
        
        
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
        df_clim["units"] = None


        # Extract environmental data and format it using its related metadata
        df_env = pd.read_json(self._env_dataformat["table_0"])
        df_env.columns = ["layer_number", "geoclim_feature"]
        df_env["geoclim_type"] = "env"    # Tag for latter id in merge

        current_start_layer = max(df_clim["layer_number"]) + 1    # Save the current state of the layer number according to the final clim_df
        df_env["layer_number"] = range(current_start_layer, current_start_layer + len(df_env))

        # Generate layer_name since absent from metadata
        df_env["layer_name"] = df_env.apply(get_env_layer_name, axis=1)

        # Assign dummy var int-val to information for all categorical data layers
        env_meta_geol = json.loads(self._env_metadata["table_0"])    # env-raster band num 5 == 'geology'
        geology_categorical = []

        for i, val in env_meta_geol["Raster value"].items():
            rock_type = env_meta_geol[list(env_meta_geol.keys())[1]][i]
            rock_type = "_".join(rock_type.split(" "))
            rock_vals_type = f"{val}={rock_type}"
            geology_categorical.append(rock_vals_type)

        env_meta_soil_vals = [v for k, v in self._env_metadata["table_1"].items() if k != "source"]
        env_meta_soil = list(zip(*env_meta_soil_vals))    # env-raster band num 6 == 'soil'
        soil_categorical = [f"{val}={'_'.join(soil.split(' '))}" for val, soil in env_meta_soil]

        env_meta_vege_vals = [v for k, v in self._env_metadata["table_2"].items() if k != "source"]
        env_meta_vegetation = list(zip(*env_meta_vege_vals))    # env-raster band num 7 == 'vegetation'
        vege_categorical = [f"{val}={'_'.join(vege.split(' '))}" for val, vege in env_meta_vegetation]
        
        env_meta_watersheds_vals = [v for k, v in self._env_metadata["table_3"].items() if k != "source"]
        env_meta_watersheds = list(zip(*env_meta_watersheds_vals))    # env-raster band num 7 == 'vegetation'
        watersheds_categorical = [f"{val}={'_'.join(watershed.split(' '))}" for val, watershed in env_meta_watersheds]
        
        # put placeholder nones for desc and units
        df_env["layer_description"] = None
        df_env["units"] = None

        # Add categorical values to units col
        df_env.at[(df_env["layer_name"] == "geo").idxmax(), "units"] = geology_categorical
        df_env.at[(df_env["layer_name"] == "soi").idxmax(), "units"] = soil_categorical
        df_env.at[(df_env["layer_name"] == "veg").idxmax(), "units"] = vege_categorical
        df_env.at[(df_env["layer_name"] == "wat").idxmax(), "units"] = watersheds_categorical
        
        # Add layer_description entries to 'null' categoricals
        df_env.at[(df_env["layer_name"] == "geo").idxmax(), "layer_description"] = list(env_meta_geol.keys())[1]
        df_env.at[(df_env["layer_name"] == "soi").idxmax(), "layer_description"] = list(self._env_metadata["table_1"].keys())[1]
        df_env.at[(df_env["layer_name"] == "veg").idxmax(), "layer_description"] = list(self._env_metadata["table_2"].keys())[1]
        df_env.at[(df_env["layer_name"] == "wat").idxmax(), "layer_description"] = list(self._env_metadata["table_3"].keys())[1]

        # Add description to None numerical environmental vars
        env_meta_others = list(zip(*self._env_metadata["table_4"].values()))
        for name, desc in env_meta_others:
            df_env.at[(df_env["layer_name"] == name).idxmax(), "layer_description"] = desc
    
        # Concat both clim and env final dfs
        df = pd.concat([df_clim, df_env])
        df = df.reset_index().drop(columns="index")

        # Extract the units to final col except unitless + categorical
        df["units"] = df.apply(extract_units, axis=1)
        df["units"] = df["units"].apply(lambda x: x.strip() if isinstance(x, str) else x)
        
        # Remove units from layer_description
        df["layer_description"] = df.apply(remove_units_from_desc, axis=1)
        df["layer_description"] = df["layer_description"].str.strip()

        # Add categorical col for easier id
        df["is_categorical"] = df["units"].apply(lambda x: isinstance(x, list))


        # Reorder cols
        df = df.loc[:, ["geoclim_type", "layer_number", "layer_name", "layer_description", "is_categorical", "units"]]

        return df
    
    def _get_categorical_df(self) -> pd.DataFrame:
        """
        Private method to extract the categorical data from the metadata of all layers of a raster database 
        related to climate and environment of Madagascar. It returns a DataFrame where each row represents a 
        distinct category within each categorical layer.

        Returns:
            pd.DataFrame: A DataFrame containing information for each categorical value in each layer. The DataFrame
                has the following columns:
                - 'geoclim_type': Type of geoclimatic data.
                - 'layer_number': Numeric identifier of the layer.
                - 'layer_name': Name of the layer.
                - 'layer_description': Description of the layer.
                - 'value': Numeric identifier of the category within the layer.
                - 'category': Description of the category within the layer.
        """
        # Extract categorical only from all layers
        all_layers_df = self._all_layers.copy()
        cat_df = all_layers_df.loc[all_layers_df["is_categorical"] == True]

        # Split units list into separate val, category columns
        cat_df = cat_df.explode("units").reset_index().drop(columns=["index", "is_categorical"])
        cat_df[["raster_value", "category"]] = cat_df["units"].str.split("=", expand=True)
        cat_df = cat_df.drop(columns="units")

        return cat_df

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
            raise AttributeError(f"Undefined attribute: '{raster_attr_name}'. You need to assign a valid pathlib.Path to the related raster attribute first.")

        # Check if raster file exists
        if not raster_to_check.is_file():
            raise FileExistsError(f"Could not find '{raster_attr_name}' file: {raster_to_check}")
                
        # Catch any IO errors
        try:
            with rasterio.open(raster_to_check) as raster_file:
                return
        except rasterio.errors.RasterioIOError as e:
            raise IOError(f"Could not open '{raster_attr_name}': {raster_to_check}. Error: {e}")
