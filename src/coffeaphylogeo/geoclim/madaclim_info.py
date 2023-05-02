import json
import re
import pathlib
import requests
import time
import inspect
from calendar import month_name
from pathlib import Path
from typing import List, Union, Optional, Tuple, Dict

import pandas as pd
import geopandas as gpd
import numpy as np
import rasterio
from shapely.geometry import Point

from coffeaphylogeo.definitions import Definitions

defs = Definitions()

class MadaclimLayers:
    """A class that represents all of the information and data from the climate and environmental variable layers that can be found from the rasters of the Madaclim database.
    
    Attributes:
        climate_dir (Path): The directory path for climate-related data.
        enviro_dir (Path): The directory path for environment-related data.
        clim_dataformat_filename (str): The file name of the climate data format file.
        clim_meta_filename (str): The file name of the climate metadata file.
        env_dataformat_filename (str): The file name of the environment data format file.
        env_meta_filename (str): The file name of the environment metadata file.
        all_layers (pd.DataFrame): A DataFrame containing a complete and formatted version of all Madaclim layers.
    
    """
    def __init__(self):
        """Initializes a new instance of the MadaclimLayers class.

        This constructor sets the directory paths and file names for climate and environment data,
        and generates a DataFrame containing all Madaclim layers.
        
        Examples:
            >>> from coffeaphylogeo.geoclim.madaclim_layers import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> # Save the all layers df
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
        """
        self.climate_dir = defs.get_geoclim_path("climate_data")    # Path dir for clim-related data
        self.enviro_dir = defs.get_geoclim_path("environment_data")    # Path dir for env-related data
        
        self.clim_dataformat_filename = defs.geoclim_files["clim_data_format"]
        self.clim_meta_filename = defs.geoclim_files["clim_metadata"]
        self.env_dataformat_filename = defs.geoclim_files["env_data_format"]
        self.env_meta_filename = defs.geoclim_files["env_metadata"]
        self.clim_raster_filename = defs.geoclim_files["madaclim_current"]
        self.env_raster_filename = defs.geoclim_files["madaclim_enviro"]
        
        self.all_layers = self._get_madaclim_layers(
            climate_dir=self.climate_dir,
            enviro_dir=self.enviro_dir,
            clim_dataformat_filename=self.clim_dataformat_filename,
            clim_meta_filename=self.clim_meta_filename,
            env_dataformat_filename=self.env_dataformat_filename,
            env_meta_filename=self.env_meta_filename
        )


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
    def clim_dataformat_filename(self):
        """str: The file name of the climate data file."""
        return self._clim_dataformat_filename
    
    @clim_dataformat_filename.setter
    def clim_dataformat_filename(self, value):
        """Sets the file name of the climate data file.

        Args:
            value (str): The file name of the climate data file.

        Raises:
            TypeError: If the provided value is not a string.
            ValueError: If the provided file does not exist in the climate data directory.
        """
        # Validate type
        if not isinstance(value, str):
            raise TypeError("clim_dataformat_filename attribute nust be a string.")
        
        # Validate path and file
        valid_file = self.climate_dir / value
        if not valid_file.exists():
            raise ValueError(f"{valid_file} does not exists")
        
        else:
            self._clim_dataformat_filename = value

    @property
    def clim_meta_filename(self):
        """str: The file name of the climate metadata file."""
        return self._clim_meta_filename
    
    @clim_meta_filename.setter
    def clim_meta_filename(self, value):
        """Sets the file name of the climate metadata file.

        Args:
            value (str): The file name of the climate metadata file.

        Raises:
            TypeError: If the provided value is not a string.
            ValueError: If the provided file does not exist in the climate data directory.
        """
        # Validate type
        if not isinstance(value, str):
            raise TypeError("clim_meta_filename attribute nust be a string.")
        
        # Validate path and file
        valid_file = self.climate_dir / value
        if not valid_file.exists():
            raise ValueError(f"{valid_file} does not exists")
        
        else:
            self._clim_meta_filename = value
    
    @property
    def env_dataformat_filename(self):
        """str: The file name of the environmental data file."""
        return self._env_dataformat_filename
    
    @env_dataformat_filename.setter
    def env_dataformat_filename(self, value):
        """Sets the file name of the environmental data file.

        Args:
            value (str): The file name of the environmental data file.

        Raises:
            TypeError: If the provided value is not a string.
            ValueError: If the provided file does not exist in the environmental data directory.
        """
        # Validate type
        if not isinstance(value, str):
            raise TypeError("env_dataformat_filename attribute nust be a string.")
        
        # Validate path and file
        valid_file = self.enviro_dir / value
        if not valid_file.exists():
            raise ValueError(f"{valid_file} does not exists")
        
        else:
            self._env_dataformat_filename = value

    @property
    def env_meta_filename(self):
        """str: The file name of the environmental metadata file."""
        return self._env_meta_filename
    
    @env_meta_filename.setter
    def env_meta_filename(self, value):
        """Sets the file name of the environmental metadata file.

        Args:
            value (str): The file name of the environmental metadata file.

        Raises:
            TypeError: If the provided value is not a string.
            ValueError: If the provided file does not exist in the environmental data directory.
        """
        # Validate type
        if not isinstance(value, str):
            raise TypeError("env_meta_filename attribute nust be a string.")
        
        # Validate path and file
        valid_file = self.enviro_dir / value
        if not valid_file.exists():
            raise ValueError(f"{valid_file} does not exists")
        
        else:
            self._env_meta_filename = value

    @property
    def clim_raster_filename(self):
        """str: The file name of the current_climate raster file."""
        return self._clim_raster_filename
    
    @clim_raster_filename.setter
    def clim_raster_filename(self, value):
        """Sets the file name of the current_climate raster file.

        Args:
            value (str): The file name of the current_climate raster file.

        Raises:
            TypeError: If the provided value is not a string.
            ValueError: If the provided file does not exist in the climate data directory.
        """
        # Validate type
        if not isinstance(value, str):
            raise TypeError("env_meta_filename attribute nust be a string.")
        
        # Validate path and file
        test_raster_file = self.climate_dir / value
        if not test_raster_file.exists():
            raise ValueError(f"{test_raster_file} does not exists")
        
        else:
            self._clim_raster_filename = value

    @property
    def env_raster_filename(self):
        """str: The file name of the environmental raster file."""
        return self._env_raster_filename
    
    @env_raster_filename.setter
    def env_raster_filename(self, value):
        """Sets the file name of the environmental raster file.

        Args:
            value (str): The file name of the environmental raster file.

        Raises:
            TypeError: If the provided value is not a string.
            ValueError: If the provided file does not exist in the environmental data directory.
        """
        # Validate type
        if not isinstance(value, str):
            raise TypeError("env_meta_filename attribute nust be a string.")
        
        # Validate path and file
        test_raster_file = self.enviro_dir / value
        if not test_raster_file.exists():
            raise ValueError(f"{test_raster_file} does not exists")
        
        else:
            self._env_raster_filename = value

    def update_all_layers(self):
        """Updates the all_layers attribute with the current values of the instance attributes.

        This method calls the _get_madaclim_layers method with the current values of the climate_dir,
        enviro_dir, clim_dataformat_filename, clim_meta_filename, env_dataformat_filename, and env_meta_filename
        attributes and assigns its return value to the all_layers attribute.

        Args:
            self (MadaclimLayers): The instance of the MadaclimLayers class.

        Returns:
            None
        """
        self.all_layers = self._get_madaclim_layers(
            climate_dir=self.climate_dir,
            enviro_dir=self.enviro_dir,
            clim_dataformat_filename=self.clim_dataformat_filename,
            clim_meta_filename=self.clim_meta_filename,
            env_dataformat_filename=self.env_dataformat_filename,
            env_meta_filename=self.env_meta_filename
        )

    def _get_madaclim_layers(self, climate_dir, enviro_dir, clim_dataformat_filename, clim_meta_filename, env_dataformat_filename, env_meta_filename) -> pd.DataFrame :
        """Private method that will generate the all_layers attributes based on the climate/enviro dirs and all the data and metada files found by accessing their corresponding attriubtes.

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
        
        
        # Open data and metadata tables
        with open(climate_dir / clim_dataformat_filename, "r") as f:
            clim_format = json.load(f)
        with open(climate_dir /  clim_meta_filename,"r") as f:
            clim_meta = json.load(f)
        with open(enviro_dir / env_dataformat_filename, "r") as f:
            env_format = json.load(f)
        with open(enviro_dir / env_meta_filename, "r") as f:
            env_meta = json.load(f)

        # Extract climate data and format it using the metadata
        df_clim = pd.read_json(clim_format["table_0"])
        df_clim["data_type"] = "clim"    # Tag for latter id in merge

        # Split the layers to get initial layer_num for all clim-related layers
        df_clim = pd.concat((df_clim.apply(split_layers, args=("Layers", ), axis=1)).to_list()).reset_index()
        df_clim.columns = ["layer_number", "geoclim_feature", "geoclim_type"]

        # Formatting + add layers to the monthly bioclim metadata
        bio_monthly_feats = pd.read_json(clim_meta["table_0"])
        bio_monthly_feats.columns = ["layer_name", "layer_description"]
        
        bio_monthly_feats = pd.concat(    # Split according to range
            bio_monthly_feats.apply(split_repeating_vars, axis=1, args=("layer_name", "layer_description", )).to_list(),
            ignore_index=True
        )
        bio_monthly_feats = add_layer_numbers_bio_monthly(bio_monthly_feats)    # Append layer numbers

        # Formatting + add layers to the other bioclim metadata (non-monthly)
        bioclim_feats = pd.read_json(clim_meta["table_1"])
        bioclim_feats.columns = ["layer_name", "layer_description"]

        current_start_layer = len(bio_monthly_feats) + 1    # Save the current state of the layer number for clim_df
        bioclim_feats["layer_number"] = range(current_start_layer, current_start_layer + len(bioclim_feats))

        # Formatting + add layers to monthly and annual evapotranspiration metadata
        evap_feats = pd.read_json(clim_meta["table_2"])
        evap_feats.columns = ["layer_name", "layer_description"]
        
        evap_feats = pd.concat(    # Split monthly evapo data
            evap_feats.apply(split_repeating_vars, axis=1, args=("layer_name", "layer_description", )).to_list(),
            ignore_index=True
        )
        current_start_layer = max(bioclim_feats["layer_number"]) + 1    # Save the current state of the layer number for clim_df
        evap_feats["layer_number"] = range(current_start_layer, current_start_layer + len(evap_feats))

        # Formatting + add layers to the bioclim water-related metadata
        biowater_feats = pd.read_json(clim_meta["table_3"])
        biowater_feats.columns = ["layer_name", "layer_description"]
        
        current_start_layer = max(evap_feats["layer_number"]) + 1    # Save the current state of the layer number for clim_df
        biowater_feats["layer_number"] = range(current_start_layer, current_start_layer + len(biowater_feats))

        # Merge meta_dfs with original clim_df for class attribute
        meta_dfs = [bio_monthly_feats, bioclim_feats, evap_feats, biowater_feats]
        df_clim = meta_merge_clim_df(df_clim, meta_dfs)


        # Extract environmental data and format it using its related metadata
        df_env = pd.read_json(env_format["table_0"])
        df_env.columns = ["layer_number", "geoclim_feature"]
        df_env["geoclim_type"] = "env"    # Tag for latter id in merge

        current_start_layer = max(df_clim["layer_number"]) + 1    # Save the current state of the layer number according to the final clim_df
        df_env["layer_number"] = range(current_start_layer, current_start_layer + len(df_env))

        # Generate layer_name since absent from metadata
        df_env["layer_name"] = df_env["geoclim_feature"].str.split(" ").str[0].str.lower()
        df_env.loc[df_env["layer_number"] == 79, "layer_name"] = "forestcover"    # Fix first word with more informative info
        df_env["layer_description"] = None

        # Assign dummy var information for geology layer to layer_description
        geology_description = {}
    
        env_meta_str = env_meta["table_0"]
        env_meta_data = json.loads(env_meta_str)

        for i, val in env_meta_data["Raster value"].items():
            rock_type = env_meta_data["Rock type"][i]
            geology_description[val] = rock_type

        df_env.loc[df_env["layer_name"] == "geology", "layer_description"] = [geology_description]    # Single-item list since dict unsupported

        # Concat both clim and env final dfs
        df = pd.concat([df_clim, df_env])

        return df
    
    def select_geoclim_type_layers(self, geoclim_type: str) -> pd.DataFrame:
        """Method that selects the desired geoclimatic type layers as a dataframe.

        Args:
            geoclim_type (str): The desired geoclimatic layers type to extract.

        Returns:
            pd.DataFrame: A slice of the all_layers dataframe containing the desired geoclimatic type layers.

        Raises:
            TypeError: If geoclim_type is not a string.
            ValueError: If geoclim_type does not corresponds to a valid geoclim type.
        
        Examples:
        >>> from coffeaphylogeo.geoclim.madaclim_layers import MadaclimLayers
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
        # Validate geoclim_type
        if not isinstance(geoclim_type, str):
            raise TypeError("geoclim_type must be a string.")
        
        possible_geoclim_types = ["clim", "env"]
        if geoclim_type not in possible_geoclim_types:
            raise ValueError(f"geoclim_type must be one of {possible_geoclim_types}")
        
        select_df = self.all_layers[self.all_layers["geoclim_type"] == geoclim_type]

        return select_df
    
    def unique_labels_layers(self, geoclim_type: str="all") -> dict:
        """Extract unique combinations from all available layers in the Madaclim db as a dictionary of unique keys and values for each layer.

        This method takes an optional `geoclim_type` argument that specifies the type of geoclimatic layers to include in the result. If `geoclim_type` is not provided, all layers are included by default.

        Args:
            geoclim_type (str, optional): The desired geoclimatic layers type to extract Defaults to "all".

        Returns:
            dict: A dictionary containing unique labels for each layer of the Madaclim db.
        
        Raises:
            TypeError: If geoclim_type is not a string.
            ValueError: If geoclim_type does not corresponds to a valid geoclim type.
            ValueError: If values for "layer_number", "layer_name" are not unique for each entry of the all_layers dataframe.

        Example:
            >>> from coffeaphylogeo.geoclim.madaclim_layers import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            #TODO REDO OUTPUT EXAMPLE
            >>> >>> madaclim_info.unique_labels_layers(geoclim_type="env")
            {'layer_71': 'env_altitude', 'layer_72': 'env_slope', 'layer_73': 'env_aspect', 'layer_74': 'env_solar', 'layer_75': 'env_geology', 'layer_76': 'env_soil', 'layer_77': 'env_vegetation', 'layer_78': 'env_watersheds', 'layer_79': 'env_forestcover'}
        """
        # Validate geoclim_type
        if not isinstance(geoclim_type, str):
            raise TypeError("geoclim_type must be a string.")
        
        possible_geoclim_types = ["clim", "env", "all"]
        if geoclim_type not in possible_geoclim_types:
            raise ValueError(f"geoclim_type must be one of {possible_geoclim_types}")
        
        # Validate unique entries
        if len(self.all_layers) != len(self.all_layers["layer_number"].unique()) != len(self.all_layers["layer_name"].unique()):
            raise ValueError("'layer_number' and 'layer_name' columns in the all_layers dataframe have non-unique entries.")
        
        # Get dict for unique labels according to geoclim_type selection
        if geoclim_type != "all":
            select_df = self.all_layers[self.all_layers["geoclim_type"] == geoclim_type]
        else: 
            select_df = self.all_layers
            
        unique_labels = {f"layer_{num}": f"{geoclim}_{name}" for num, geoclim, name in list(zip(select_df["layer_number"], select_df["geoclim_type"], select_df["layer_name"]))}
        return unique_labels
        

    def fetch_specific_layers(self, layer_numbers: Union[int, str, List[Union[int, str]]], from_unique_labels: bool=False, as_descriptive_labels: bool=False) -> Union[dict, pd.DataFrame]:
        """Fetches specific layers from the all_layers DataFrame.

        Args:
            layer_numbers (Union[int, str, List[Union[int, str]]]): The layer number(s) to fetch. Can be a single int or str value or a list of int or str values.
            from_unique_labels (bool): If True, layer_numbers are treated as keys from the unique_labels method output. Defaults to False.
            as_descriptive_labels (bool): If True, only the layer description is returned. Defaults to False.

        Returns:
            Union[dict, pd.DataFrame]: If as_descriptive_labels is True, returns a dictionary with the layer descriptions. Otherwise, returns a DataFrame with the specified layers.

        Raises:
            ValueError: If from_unique_labels is True and any value in layer_numbers is not a valid key from the unique_labels method output.
            TypeError: If any value in layer_numbers cannot be converted to an int.

        Examples:
        >>> # From a single list of ints
        >>> from coffeaphylogeo.geoclim.madaclim_layers import MadaclimLayers
        >>> madaclim_info = MadaclimLayers()
        >>> madaclim_info.fetch_specific_layers([1, 15, 55])
            layer_number                        geoclim_feature geoclim_type layer_name                                layer_description
        0              1  Monthly minimum temperature (°C x 10)         clim      tmin1  Monthly minimum temperature (°C x 10) - January
        14            15  Monthly maximum temperature (°C x 10)         clim      tmax3    Monthly maximum temperature (°C x 10) - March
        54            55        Bioclimatic variables (bioclim)         clim      bio19                 Precipitation of coldest quarter
        
        >>> # From a unique label set
        >>> unique_env_labels = list(madaclim_info.unique_labels_layers(geoclim_type="env").keys())
        >>> madaclim_info.fetch_specific_layers(unique_env_labels, from_unique_labels=True)
        layer_number                                    geoclim_feature  ...   layer_name                                  layer_description
        0            71                                       Altitude (m)  ...     altitude                                               None
        1            72                                  Slope (in degree)  ...        slope                                               None
        2            73           Aspect (clockwise from North, in degree)  ...       aspect                                               None
        3            74                     Solar radiation (Wh.m-2.day-1)  ...        solar                                               None
        4            75               Geology (Kew Botanical Garden, 1997)  ...      geology  {1: 'Alluvial & Lake deposits', 2: 'Unconsolid...
        5            76                             Soil (Pelletier, 1981)  ...         soil                                               None
        6            77            Vegetation (Kew Botanical Garden, 2007)  ...   vegetation                                               None
        7            78                         Watersheds (Pearson, 2009)  ...   watersheds                                               None
        8            79  Percentage of forest cover for the year 2010 (%).  ...  forestcover                                               None

        [9 rows x 5 columns]
        #TODO REDO EXAMPLE FOR AS_DESCRIPTIVE
        
        """
        # Validate layer_numbers
        if isinstance(layer_numbers, list):
            if from_unique_labels:    # layer_numbers same as keys from unique_labels method output
                possible_unique_labels = [f"layer_{num}" for num in self.all_layers["layer_number"].to_list()]
                for layer_number in layer_numbers:
                    if layer_number not in possible_unique_labels:
                        raise ValueError(f"{layer_number} not one of {possible_unique_labels[:3]}...{possible_unique_labels[:-3]}")
                # Save as list of ints after check
                layer_numbers = [int(layer.split("_")[1]) for layer in layer_numbers]

            # layer_numbers as list of ints
            else:
                try:
                    layer_numbers = [int(layer) for layer in layer_numbers]
                except (ValueError, TypeError):
                    raise TypeError("layer_numbers must be either a single int value or a string that can be converted to an int, or a list of int values or strings that can be converted to int values")
        # Single layer_numbers type check 
        else:
            try:
                layer_numbers = [int(layer_numbers)]
            except (ValueError, TypeError):
                raise TypeError("layer_numbers must be either a single int value or a string that can be converted to an int, or a list of int values or strings that can be converted to int values")
            
        # Validate layer number range for layer_numbers as int(s)
        min_layer = min(self.all_layers["layer_number"])
        max_layer = max(self.all_layers["layer_number"])

        for layer_number in layer_numbers:
            if not min_layer <= layer_number <= max_layer:
                raise ValueError(f"layer_number must fall between {min_layer} and {max_layer} (You entered {layer_number=}).")
            
        # Fetch rows according to layer selection
        if as_descriptive_labels:
            # Save subset df with selected layers
            select_df = self.all_layers[self.all_layers["layer_number"].isin(layer_numbers)]
            sub_selection = list(zip(select_df["layer_number"], select_df["geoclim_type"], select_df["layer_name"], select_df["layer_description"]))

            # Generate dict with key as layer_num and values containing layer information
            description = {f"layer_{num}": f"{geotype}_{num}_{name} ({desc})" for num, geotype, name, desc in sub_selection}
            
            return description
        
        else:    # Save whole row of df
            select_df = self.all_layers[self.all_layers["layer_number"].isin(layer_numbers)]
            return select_df

    def download_data(self, save_dir: Optional[pathlib.Path]=None):
        """Downloads climate and environment data from the Madaclim website.

        This method downloads the climate and environment raster data from the Madaclim website
        and saves them to the specified directory. If no directory is specified, the data is saved
        to the default directory.

        Args:
            save_dir (Optional[pathlib.Path]): The directory where the data should be saved. If not specified,
                the data is saved to the default directory.

        Raises:
            ValueError: If save_dir is not a directory.
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
        
        # Validate save_dir when specified
        if save_dir is not None:
            if not save_dir.is_dir():
                raise ValueError("save_dir is not a directory.")
        
        # Download rasters with default save_dir
        if save_dir is None:
            download_single_file(    # Climate raster
                url=defs.urls["madaclim_current_raster"],
                dir_savepath=self.climate_dir,
                filename=self.clim_raster_filename
            )

            download_single_file(
                url=defs.urls["environment_raster"],
                dir_savepath=self.enviro_dir,
                filename=self.env_raster_filename
            )
        # Download to specified path
        else:
            download_single_file(    # Climate raster
                url=defs.urls["madaclim_current_raster"],
                dir_savepath=save_dir,
                filename=self.clim_raster_filename
            )

            download_single_file(
                url=defs.urls["environment_raster"],
                dir_savepath=save_dir,
                filename=self.env_raster_filename
            )

    def sample_rasters_from_gdf(self, gdf: gpd.GeoDataFrame, geometry_col_name: str="geometry", as_descriptive_labels: bool=True)->Tuple[gpd.GeoDataFrame, Dict[str, List[np.ndarray]]]:
        """Samples raster values from a GeoDataFrame containing Point geometries.

        Args:
            gdf (gpd.GeoDataFrame): The GeoDataFrame containing the Point geometries to sample from.
            geometry_col_name (str, optional): The name of the column in the GeoDataFrame that contains the Point geometries. Defaults to "geometry".
            as_descriptive_labels (bool, optional): Whether to return descriptive labels for the raster values. Defaults to True.

        Raises:
            ValueError: If the geometry column in the GeoDataFrame is not a Point geometry or if it contains empty Point objects.
            ValueError: If the raster projections are on different coordinate reference systems (CRS).

        Returns:
            Tuple[gpd.GeoDataFrame, Dict[str, List[np.ndarray]]]: A tuple containing a copy of the input GeoDataFrame with additional columns for the raster values and a dictionary containing the raster values for each filename.
        #TODO ADD EXAMPLES!
        """
        # Make a copy of the geodf
        gdf_copy = gdf.copy()

        # Get the geometry column from the GeoDataFrame
        geom_data = gdf_copy[geometry_col_name]

        # Check if the geometry column is a Point geometry
        if geom_data.dtype != 'geometry' or geom_data.geom_type.unique()[0] != 'Point':
            raise ValueError(f"The '{geometry_col_name}' column in the GeoDataFrame is not a Point geometry.")
        
        # Check if geodataframe contains empty Point objects
        empty_points = (geom_data.is_empty).sum()
        if empty_points > 0:
            raise ValueError("Empty Point objects cannot be sampled. To remove them use:\ngpd.GeoDataFrame(gdf.loc[gdf[geometry].is_empty == False])")

        # Define rasters path
        raster_clim = rasterio.open(self.climate_dir / self.clim_raster_filename)
        raster_env = rasterio.open(self.enviro_dir / self.env_raster_filename)

        # Sanity check for rasters
        if raster_clim.crs != raster_env.crs:
            raise ValueError("Raster projections are on different coordinate reference system (CRS).")

        # Validate the CRS from the geodataframe to both rasters
        if gdf_copy.crs != raster_clim.crs or gdf_copy.crs!= raster_env.crs:
            try:
                print(f"Reprojecting gdf[geometry] Point objects from {gdf_copy.crs} to the rasters' {raster_clim.crs}...\n")
                gdf_copy = gdf_copy.to_crs(raster_clim.crs)
            except Exception as e:
                raise(f"Error {e}: Could not convert GeoDataFrame to the raster's CRS.")

        # Extract nodata value for each rasters
        nodata_clim = raster_clim.nodata
        nodata_env = raster_env.nodata

        # Get geometry points coordinates
        coord_list = [(x,y) for x,y in zip(gdf_copy[geometry_col_name].x, gdf_copy[geometry_col_name].y)]

        # Sample rasters from Point objects
        raster_samples = {}
        raster_samples[self.clim_raster_filename] = [x for x in raster_clim.sample(coord_list)]    # Climate raster sampling        
        print(f"Extracted {len(raster_samples[self.clim_raster_filename])} specimens from {self.clim_raster_filename}")

        raster_samples[self.env_raster_filename] = [x for x in raster_env.sample(coord_list)]    # Environmental raster sampling
        print(f"Extracted {len(raster_samples[self.env_raster_filename])} specimens from {self.env_raster_filename}\n")

        # Catch specimens with nodata values
        if nodata_clim is None:    # Climate nodata samples 
            print(f"Raster {self.clim_raster_filename} does not contain nodata values")
        else:
            num_nodata_clim = sum([nodata_clim in val for val in raster_samples[self.clim_raster_filename]])
            if num_nodata_clim > 0 :
                print(f"BEWARE! {num_nodata_clim} specimens contains 'nodata' values when extracting from {self.clim_raster_filename}:")
                print(f"Index of specimen(s) with at least 1 'nodata' point: {[i for i, val in enumerate(raster_samples[self.clim_raster_filename]) if nodata_clim in val]}\n")
            else:
                print(f"0 specimen containing nodata values in any of the {len(raster_clim.indexes)} variables from {self.clim_raster_filename}\n")
        
        if nodata_env is None:    # Environmental nodata samples
            print(f"Raster {self.env_raster_filename} does not contain nodata values")
        else:
            num_nodata_env = sum([nodata_env in val for val in raster_samples[self.env_raster_filename]])    
            if num_nodata_env > 0 :
                print(f"BEWARE! {num_nodata_env} specimens contains 'nodata' values when extracting from {self.env_raster_filename}:")
                print(f"Index of specimen(s) with at least 1 'nodata' point: {[i for i, val in enumerate(raster_samples[self.env_raster_filename]) if nodata_env in val]}\n")
            else:
                print(f"0 specimen containing nodata values in any of the {len(raster_env.indexes)} variables from {self.env_raster_filename}\n")
        
        # Assign a new column with the list arrays sampled
        gdf_copy[f"raster_samples[{self.clim_raster_filename}]"] = raster_samples[self.clim_raster_filename]
        gdf_copy[f"raster_samples[{self.env_raster_filename}]"] = raster_samples[self.env_raster_filename]

        # Extract layer_info to use as col name when transposing raster samples to geodataframe
        clim_number_layers = self.select_geoclim_type_layers("clim")["layer_number"].to_list()
        clim_layers_info = self.fetch_specific_layers(layer_numbers=clim_number_layers, as_descriptive_labels=True)

        env_number_layers = self.select_geoclim_type_layers("env")["layer_number"].to_list()
        env_layers_info = self.fetch_specific_layers(layer_numbers=env_number_layers, as_descriptive_labels=True)

        # Append sampled data to geodataframe
        for index, (k, v) in enumerate(clim_layers_info.items()):
            for val, specimen_index in zip(raster_samples[self.clim_raster_filename], gdf_copy.index.to_list()):
                col_layer_name = v if as_descriptive_labels else k    # layer num as col names or more descriptive
                gdf_copy.loc[specimen_index, col_layer_name] = val[index]
        
        for index, (k, v) in enumerate(env_layers_info.items()):
            for val, specimen_index in zip(raster_samples[self.env_raster_filename], gdf_copy.index.to_list()):
                col_layer_name = v if as_descriptive_labels else k    # layer num as col names or more descriptive
                gdf_copy.loc[specimen_index, col_layer_name] = val[index]

        # Close rasters and return geodataframe + sampled raster values
        raster_clim.close()
        raster_env.close()
        return gdf_copy, raster_samples
    
    #TODO DEF SAMPLE_FROM_SINGLE_POINT()


    def __str__(self) -> str:
        # Get instance attributes
        attribute_keys = list(self.__dict__.keys())
        
        # Extract custom methods
        extract_methods = inspect.getmembers(self, predicate=inspect.ismethod)
        custom_methods = [method[0] for method in extract_methods if not method[0] == "__init__"]

        return f"Instance attributes:\n{attribute_keys}\n\nCustom methods:\n{custom_methods}\n"
    