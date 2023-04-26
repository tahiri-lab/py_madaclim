import json
import re
import pathlib
import pandas as pd

from calendar import month_name
from pathlib import Path
from typing import List

from coffeaphylogeo.definitions import Definitions

defs = Definitions()

class MadaclimLayers:
    """
    A class that represents all of the information and data from the climate and environmental variable layers that can be found 
    from the rasters of the Madaclim database.
    
    Attributes:
        climate_dir (Path): The directory path for climate-related data.
        enviro_dir (Path): The directory path for environment-related data.
        clim_data_file (str): The file name of the climate data format file.
        clim_meta_file (str): The file name of the climate metadata file.
        env_data_file (str): The file name of the environment data format file.
        env_meta_file (str): The file name of the environment metadata file.
        all_layers (pd.DataFrame): A DataFrame containing a complete and formatted version of all Madaclim layers.
    
    """
    def __init__(self):
        """Initializes a new instance of the MadaclimLayers class.

        This constructor sets the directory paths and file names for climate and environment data,
        and generates a DataFrame containing all Madaclim layers.
        """
        self.climate_dir = defs.get_geoclim_path("climate_data")    # Path dir for clim-related data
        self.enviro_dir = defs.get_geoclim_path("environment_data")    # Path dir for env-related data
        
        self.clim_data_file = defs.geoclim_files["clim_data_format"]
        self.clim_meta_file = defs.geoclim_files["clim_metadata"]
        self.env_data_file = defs.geoclim_files["env_data_format"]
        self.env_meta_file = defs.geoclim_files["env_metadata"]
        
        self.all_layers = self._get_madaclim_layers(
            climate_dir=self.climate_dir,
            enviro_dir=self.enviro_dir,
            clim_data_file=self.clim_data_file,
            clim_meta_file=self.clim_meta_file,
            env_data_file=self.env_data_file,
            env_meta_file=self.env_meta_file
        )

        # self.geology_raster_vars = {}

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

    def _get_madaclim_layers(self, climate_dir, enviro_dir, clim_data_file, clim_meta_file, env_data_file, env_meta_file) -> pd.DataFrame :
        """
        Private method that will generate the all_layers attributes based on the climate/enviro dirs and all the data and metada
        files found by accessing their corresponding attriubtes.

        Args:
            climate_dir (Path): The directory path for climate-related data.
            enviro_dir (Path): The directory path for environment-related data.
            clim_data_file (str): The file name of the climate data format file.
            clim_meta_file (str): The file name of the climate metadata file.
            env_data_file (str): The file name of the environment data format file.
            env_meta_file (str): The file name of the environment metadata file.

        Returns:
            pd.DataFrame: A DataFrame containing a complete and formatted version of all Madaclim layers.
        """
        
        def split_layers(row, layers_col_name: str):
            """
            Split "Layers" column and create new rows

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

    