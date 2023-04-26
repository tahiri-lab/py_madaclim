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
            """
            Split a column containing repeating values into separate rows.

            This function takes a row of data from a DataFrame and splits the value in the column specified by `col_to_split` if it contains a hyphen. The function returns a new DataFrame with rows for each value in the specified range and with values from the column specified by `col_to_keep`.

            Args:
                row (pd.Series): A row of data from a DataFrame.
                col_to_split (str): The name of the column to split.
                col_to_keep (str): The name of the column to keep.

            Returns:
                pd.DataFrame: A DataFrame with split rows.

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
            
            # Extract the range and layername to a new smaller df of len(range(start, end))
            if "-" in row[col_to_split]:
                start = int(re.search("\d+", row[col_to_split].split("-")[0]).group())
                end = int(row[col_to_split].split("-")[1])
                name = re.search("[a-z]*", row[col_to_split]).group()
                        
                # Create a DataFrame with the split values and the description column
                df = pd.DataFrame({col_to_split: [f"{name}{month}" for month in range(start, end+1)],
                                col_to_keep: row[col_to_keep]})
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
                category_feature_val = df[df["layer_name"].str.contains(layer_name)]["layer_description"].unique()[0]
                monthly_features[f"{layer_name}_cat_feature"] = category_feature_val
                
                # Get df associated with common category
                category_feature_df = pd.DataFrame(df[df["layer_name"].str.contains(layer_name)])
                monthly_features[f"{layer_name}_df"] = category_feature_df

                # Assign layer number according to current category
                df.loc[df["layer_description"] == category_feature_val, "layer_number"] = range(current_start_layer, current_start_layer + len(category_feature_df))
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
    