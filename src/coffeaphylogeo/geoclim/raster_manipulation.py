import csv
import pathlib
import inspect
from pathlib import Path
from typing import Optional, Union, List, Optional, Tuple, Dict
import time
from tqdm import tqdm

from coffeaphylogeo.definitions import Definitions
from coffeaphylogeo.geoclim.madaclim_info import MadaclimLayers

import rasterio
import pyproj
import shapely
from pyproj import Transformer
from shapely import Point
import pandas as pd


# Default dir and path for rasters
defs = Definitions()

climate_dir = defs.get_geoclim_path("climate_data")
clim_raster_filename = defs.geoclim_files["madaclim_current"]
default_clim_raster_path = climate_dir / clim_raster_filename

enviro_dir = defs.get_geoclim_path("environment_data") 
env_raster_filename = defs.geoclim_files["madaclim_enviro"]
default_env_raster_path = enviro_dir / env_raster_filename


class MadaclimPoint:
    
    """
    A class representing a specimen as a geographic point with a specific coordinate reference system (CRS)
    and additional attributes. The class provides methods for validating the point's coordinates
    and CRS, as well as sampling values from climate and environmental rasters of the Madaclim database.
    
    Attributes:
        specimen_id (str): An identifier for the point.
        source_crs (pyproj.crs.crs.CRS): The coordinate reference system of the point.
        latitude (float): The latitude of the point.
        longitude (float): The longitude of the point.
        mada_geom_point (shapely.geometry.point.Point): A Shapely Point object representing the point projected in the Madaclim rasters' CRS.
        ___base_attr (dict): A dictionary containing the base attributes names as keys and their values as values.
        
    """

    def __init__(self, specimen_id: str, latitude: float, longitude: float, source_crs: pyproj.crs.crs.CRS=pyproj.CRS.from_epsg(4326), **kwargs) -> None:
        """
        Initialize a MadaclimPoint object with the given specimen_id, latitude, longitude, and source_crs. The coordinates provided should respect the nature of the given source CRS native units' (i.e. degrees WGS84 or meters for EPSG:3857)
        Optionally, provide additional keyword arguments to store as instance attributes.
        
        Args:
            specimen_id (str): An identifier for the point.
            latitude (float): The latitude of the point.
            longitude (float): The longitude of the point.
            source_crs (pyproj.crs.crs.CRS, optional): The coordinate reference system of the point. Defaults to WGS84 (EPSG:4326).
            **kwargs: Additional keyword arguments to store as instance attributes.
        Examples:
            >>> from coffeaphylogeo.geoclim.raster_manipulation import MadaclimPoint
            >>> specimen_1 = MadaclimPoint(specimen_id="spe1_aren", latitude=-18.9333, longitude=48.2)    # Default CRS of EPSG:4326
            
            >>> # Creates a shapely point object according to the Madaclim CRS' projection when instantiating a new instance. 
            >>> specimen_1.mada_geom_point
            <POINT (837072.915 7903496.321)>

            >>> # Also accepts any other kwargs and saves them as attributes with specific typing
            >>> specimen_1 = MadaclimPoint(specimen_id="spe1_aren", latitude=-18.9333, longitude=48.2, genus="Coffea", species="arenesiana", has_sequencing=True)
            >>> specimen_1.species
            'arenesiana'
            >>> specimen_1.has_sequencing
            True
            >>> print(specimen_1)
            MadaclimPoint(
                genus = Coffea,
                mada_geom_point = POINT (837072.9150244407 7903496.320897499),
                source_crs = 4326,
                specimen_id = spe1_aren,
                species = arenesiana,
                has_sequencing = True,
                longitude = 48.2,
                latitude = -18.9333
            )

        """
        
        self.specimen_id = specimen_id
        self.source_crs = self.validate_crs(source_crs)
        self.latitude = self.validate_lat(latitude, crs=self.source_crs)
        self.longitude = self.validate_lon(longitude, crs=self.source_crs)
        self.mada_geom_point = self._construct_point(
            latitude=self.latitude, 
            longitude=self.longitude,
            source_crs=self.source_crs
        )
        
        # Store base attributes names/vals
        self.__base_attr = {k:v for k,v in self.__dict__.items()}

        # Store any additional keyword arguments as instance attributes
        base_args = self.get_args_names()[0] + [self.get_args_names()[1]]
        additional_args = [key for key in kwargs if key not in base_args]
        for key in additional_args:
            setattr(self, key, kwargs[key])
    
    @property
    def specimen_id(self) -> str:
        """
        Get the specimen_id attribute.

        Returns:
            str: The identifier for the MadaclimPoint.
        """
        return self._specimen_id
    @specimen_id.setter
    def specimen_id(self, value: str):
        """
        Set the specimen_id attribute.

        Args:
            value (str): The new identifier for the MadaclimPoint.
        """
        self._specimen_id = value
    
    @property
    def latitude(self) -> float:
        """
        Get the latitude attribute.

        Returns:
            float: The latitude of the point.
        """
        return self._latitude
    
    @latitude.setter
    def latitude(self, value: float):
        """
        Set the latitude attribute after validating the input value.

        Args:
            value (float): The latitude value for the point.
        """
        value = self.validate_lat(value, crs=self.source_crs)
        self._latitude = value
        # Update mada_geom_point when latitude is updated
        self._update_mada_geom_point()
    
    @property
    def longitude(self) -> float:
        """
        Get the longitude attribute.

        Returns:
            float: The longitude of the point.
        """
        return self._longitude
    
    @longitude.setter
    def longitude(self, value: float):
        """
        Set the longitude attribute after validating the input value.

        Args:
            value (float): The longitude value for the point.
        """
        value = self.validate_lon(value, crs=self.source_crs)
        self._longitude = value
        # Update mada_geom_point when longitude is updated
        self._update_mada_geom_point()
    
    @property
    def source_crs(self) -> pyproj.crs.CRS:
        """
        Get the source_crs attribute.

        Returns:
            pyproj.crs.CRS: The coordinate reference system of the point.
        """
        return self._source_crs
    
    @source_crs.setter
    def source_crs(self, value: pyproj.crs.CRS):
        """
        Set the source_crs attribute after validating the input value.

        Args:
            value (pyproj.crs.CRS): The coordinate reference system for the point.
        """
        value = self.validate_crs(value)
        self._source_crs = value
        
        # Update mada_geom_point when crs is updated
        self._update_mada_geom_point()

    @property
    def base_attr(self) -> dict:
        """Get the base attributes when constructing the instance

        Returns:
            dict: A dictionary containing the base attributes names as keys and their values as values.
        """
        return self.__base_attr

    def __str__(self) -> str:
        
        # Fetch base and additional attributes at construction
        base_attr = self.base_attr
        add_attr = self._get_additional_attributes()
        all_attr = {**base_attr, **add_attr}    # Append attr dictionaries
        
        # Pretty format
        all_attr_list = [f"{k.lstrip('_')} = {v.to_epsg() if k == '_source_crs' else v}" for k, v in all_attr.items()]
        all_attr_str = ",\n\t".join(all_attr_list)
        madapoint_obj = f"MadaclimPoint(\n\t{all_attr_str}\n)"
        
        return madapoint_obj

    def __repr__(self) -> str:
        return self.__str__()
    
    @staticmethod
    def get_args_names() -> Tuple[list, list]:
        """Gets the names of the required and default arguments of the MadaclimPoint constructor.

        This method uses the inspect module to introspect the MadaclimPoint constructor and extract the names of its arguments.
        It then separates these into required arguments (those that don't have default values) and default arguments (those that do).

        Returns:
            Tuple[list, list]: A tuple containing two lists:
                - The first list contains the names of the required arguments.
                - The second list contains the names of the default arguments.
                
        Note:
            - 'self' is excluded from the returned lists.
        """
        argspec = inspect.getfullargspec(MadaclimPoint)
        
        # Get required args names
        num_defaults = len(argspec.defaults) if argspec.defaults else 0
        num_required = len(argspec.args) - num_defaults
        required_args = argspec.args[1: num_required]    # Exclude self

        # Get defaults args names
        default_args = argspec.args[-num_defaults] if argspec.defaults else []
        return required_args, default_args
    
    @staticmethod
    def get_default_source_crs(as_epsg: bool=True) -> Union[pyproj.crs.crs.CRS, int]:
        """Extracts the default value of the source_crs attribute. By default, it will return the crs as the EPSG code.

        Args:
            as_epsg (bool, optional): The EPSG code of the source_CRS. Defaults to True.

        Returns:
            Union[pyproj.crs.crs.CRS, int]: The default value for the source_crs attribute. If true, source_crs is returned as the EPSG code of the crs.
        """
        source_crs_val = inspect.getfullargspec(MadaclimPoint).defaults[0]    # only 1 default args
        if as_epsg:
            source_crs_val = source_crs_val.to_epsg()
        return source_crs_val
    
    def validate_crs(self, crs):
        """
        Validate the input CRS and return a valid CRS object.

        Args:
            crs (pyproj.crs.CRS): The input CRS to validate.

        Returns:
            pyproj.crs.CRS: A valid CRS object.

        Raises:
            ValueError: If the input CRS is invalid.
        """
        try:
            valid_crs = pyproj.crs.CRS(crs)
        except pyproj.exceptions.CRSError:
            crs_error = (
                f"Invalid CRS: {crs}. For simplicy, validate your crs if using EPSG codes by using:\n"
                f">>> your_crs in [int(code) for code in pyproj.get_codes('EPSG', 'CRS')]\n"
                f"Or check the PyProj documentation: https://pyproj4.github.io/pyproj/"
            )
            raise ValueError(crs_error)
        return valid_crs
    
    def validate_lat(self, latitude, crs):
        """
        Validate the input latitude value and return a valid latitude.

        Args:
            latitude (float): The input latitude value.
            crs (pyproj.crs.CRS): The coordinate reference system of the point.

        Returns:
            float: A valid latitude value.

        Raises:
            TypeError: If the input latitude value cannot be converted to a float.
            ValueError: If the input latitude value is out of bounds for the given CRS.
        """
        # Validate float type
        try:
            latitude = float(latitude)
        except:
            raise TypeError(f"Could not convert {latitude} to float. Latitude must be a float.")
        
        # Validate max lat according to crs bounds
        if crs.is_geographic:
            bounds = crs.area_of_use.bounds
        else:
            # Extract bounds in native units of projection
            transformer = Transformer.from_crs(crs.geodetic_crs, crs, always_xy=True)
            bounds = transformer.transform_bounds(*crs.area_of_use.bounds)
        
        min_lat, max_lat = bounds[1], bounds[3]
        if not min_lat <= latitude <= max_lat:
            raise ValueError(f"{latitude=} is out of bounds for the crs=EPSG:{crs.to_epsg()}. Latitude must be between: {min_lat} and {max_lat}")
        
        return latitude
    
    def validate_lon(self, longitude, crs):
        """
        Validate the input longitude value and return a valid longitude.

        Args:
            longitude (float): The input longitude value.
            crs (pyproj.crs.CRS): The coordinate reference system of the point.

        Returns:
            float: A valid longitude value.

        Raises:
            TypeError: If the input longitude value cannot be converted to a float.
            ValueError: If the input longitude value is out of bounds for the given CRS.
        """
        # Validate float type
        try:
            longitude = float(longitude)
        except:
            raise TypeError(f"Could not convert {longitude} to float. longitude must be a float.")
        
        # Validate max lon according to crs bounds
        if crs.is_geographic:
            bounds = crs.area_of_use.bounds
        else:
            # Extract bounds in native units of projection
            transformer = Transformer.from_crs(crs.geodetic_crs, crs, always_xy=True)
            bounds = transformer.transform_bounds(*crs.area_of_use.bounds)

        min_lon, max_lon = bounds[0], bounds[2]
        if not min_lon <= longitude <= max_lon:
            raise ValueError(f"{longitude=} is out of bounds for the crs=EPSG:{crs.to_epsg()}. Longitude must be between {min_lon} and {max_lon}")
        
        return longitude
    
    def sample_from_rasters(
            self,
            layers_to_sample: Union[int, str, List[Union[int, str]]]="all", 
            layer_info: bool=False,
            return_nodata_layers: bool=False,
            clim_raster_path: Optional[pathlib.Path]=None, 
            env_raster_path: Optional[pathlib.Path]=None
        ) -> Union[dict, list]:
        """
        Samples geoclimatic data from raster files for specified layers at the location of the instances's lat/lon coordinates from the mada_geom_point attribute.

        Args:
            layers_to_sample (Union[int, str, List[Union[int, str]]], optional): The layer number(s) to sample from the raster files.
                Can be a single int, a single string in the format 'layer_<num>', or a list of ints or such strings. Defaults to 'all'.
            layer_info (bool, optional): Whether to use descriptive labels for the returned dictionary keys. Defaults to False.
            return_nodata_layers (bool, optional): Whether to return a list of layers with nodata values at the specimen location.
                Defaults to False.
            clim_raster_path (Optional[pathlib.Path], optional): Path to the climate raster file. Defaults to None.
            env_raster_path (Optional[pathlib.Path], optional): Path to the environment raster file. Defaults to None.

        Raises:
            TypeError: If the layers_to_sample is not valid, or if the mada_geom_point attribute is not a Point object.
            ValueError: If the layer_number is out of range or if the mada_geom_point object is empty.

        Returns:
            Union[dict, list]: A dictionary containing the sampled data, with keys being layer names or numbers depending
                on the layer_info parameter. If return_nodata_layers is True, also returns a list of layers with nodata values
                at the specimen location.
        Exmaples:
            >>> # Fetching bioclim layers from the MadaclimLayers class
            >>> from coffeaphylogeo.geoclim.madaclim_info import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> bioclim_labels = [label for label in madaclim_info.get_layers_labels(as_descriptive_labels=True) if "bio" in label]

            >>> # Sampling the bioclim layers
            >>> specimen_1 = MadaclimPoint(specimen_id="spe1_aren", latitude=-18.9333, longitude=48.2, genus="Coffea", species="arenesiana", has_sequencing=True)
            >>> spe1_arenesia_bioclim = specimen_1.sample_from_rasters(bioclim_layers)
            >>> {layer: val for layer, val in list(spe1_arenesia_bioclim.items())[:2]}
            {'clim_37_bio1 (Annual mean temperature)': 196, 'clim_38_bio2 (Mean diurnal range (mean of monthly (max temp - min temp)))': 112}

            >>> # Sample from both clim and env rasters with just a layer number
            >>> specimen_1.sample_from_rasters([68, 71])

            ######################################## Extracting data for: spe1_aren ########################################

            Sampling 1 layer(s) from madaclim_current.tif (geoclim_type=clim)...
            Extracting layer 68: pet (Annual potential evapotranspiration from the Thornthwaite equation (mm)):  1

            Sampling 1 layer(s) from madaclim_enviro.tif (geoclim_type=env)...
            Extracting layer 71: altitude (None):  100%|███████████████████████| layer 1/1 [Time remaining: 00:00]

            Finished raster sampling operation in 0.44 seconds.

            {'clim_68_pet (Annual potential evapotranspiration from the Thornthwaite equation (mm))': 891, 'env_71_altitude (None)': 899}

            >>> # Sample all layers with less descriptive layer names
            >>> specimen_2 = MadaclimPoint(specimen_id="spe2_humb", latitude=-12.716667	, longitude=45.066667, source_crs=4326)
            >>> spe2_all_layers = specimen_2.sample_from_rasters(layer_info=False)

            ######################################## Extracting data for: spe2 ########################################

            Sampling 70 layer(s) from madaclim_current.tif (geoclim_type=clim)...
            Extracting layer 70: ndm (Number of dry months in the year):  100%|█████████████████████| layer 70/70 [Time remaining: 00:00]

            Sampling 9 layer(s) from madaclim_enviro.tif (geoclim_type=env)...
            Extracting layer 79: forestcover (None):  100%|███████████████████████████████████████████| layer 9/9 [Time remaining: 00:00]
            BEWARE! 5 layer(s) contain a nodata value at the specimen location

            Finished raster sampling operation in 27.58 seconds.
            >>> {layer: val for layer, val in list(spe2_all_layers.items())[:5]}
            {'layer_1': 197, 'layer_2': 218, 'layer_3': 218, 'layer_4': 216, 'layer_5': 207}
            >>> spe2_all_layers["layer_37"]
            238

            >>> # We can easily access the nodata layers (still sampled with the method regardless)
            >>> spe2_all_layers, spe2_nodata_layers = specimen_2.sample_from_rasters("all", return_nodata_layers=True)
            >>> len(spe2_nodata_layers)
            5
            >>> spe2_nodata_layers
            ['env_75_geology (1=Alluvial_&_Lake_deposits, 2=Unconsolidated_Sands, 4=Mangrove_Swamp, 5=Tertiary_Limestones_+_Marls_&_Chalks, 6=Sandstones, 7=Mesozoic_Limestones_+_Marls_(inc._"Tsingy"), 9=Lavas_(including_Basalts_&_Gabbros), 10=Basement_Rocks_(Ign_&_Met), 11=Ultrabasics, 12=Quartzites, 13=Marble_(Cipolin))', 'env_76_soil (None)', 'env_77_vegetation (None)', 'env_78_watersheds (None)', 'env_79_forestcover (None)']
        """
        
        # Create a MadaclimLayers instance to get layers labels and validate layers to sample
        madaclim_info = MadaclimLayers()
        all_layers_df = madaclim_info.all_layers
        
        # Validate layers to sample
        possible_layers_num_format = [f"layer_{num}" for num in all_layers_df["layer_number"].to_list()]
        possible_layers_desc_format = madaclim_info.get_layers_labels(as_descriptive_labels=True)


        if isinstance(layers_to_sample, list):
            # Check if all elements are in layer_<num> format
            layers_num_format = all([layer_label in possible_layers_num_format for layer_label in layers_to_sample])
            layers_desc_format = all([layer_label in possible_layers_desc_format for layer_label in layers_to_sample])

            
            if layers_num_format or layers_desc_format:
                # Save as list of ints after check
                layers_numbers = [int(layer_label.split("_")[1]) for layer_label in layers_to_sample]

            # layers_to_sample as list of ints
            else:
                try:
                    layers_numbers = [int(layer) for layer in layers_to_sample]
                except (ValueError, TypeError):
                    raise TypeError("layers_to_sample must be either a single int value or a string that can be converted to an int, or a list of int values or strings that can be converted to int values")
        
        # Single layers_to_sample type check 
        else:
            if layers_to_sample == "all":    # Get all layers as default
                layers_numbers = all_layers_df["layer_number"].to_list()
            
            elif layers_to_sample in possible_layers_num_format:    # Check layer_<num> str format
                layers_numbers = [int(layers_to_sample.split("_")[1])]
            else:
                try:
                    layers_numbers = [int(layers_to_sample)]    # As single item int list
                except (ValueError, TypeError):
                    raise TypeError("layers_to_sample must be either a single int value or a string that can be converted to an int")
  
        # Validate layer number range for layer_numbers as in
        min_layer = min(all_layers_df["layer_number"])
        max_layer = max(all_layers_df["layer_number"])

        for layer_number in layers_numbers:
            if not min_layer <= layer_number <= max_layer:
                raise ValueError(f"layer_number must fall between {min_layer} and {max_layer}. {layer_number=} is not valid.")
            
        # Get possible layer numbers for each raster
        geoclim_types = ["clim", "env"]
        geoclim_layer_ranges = {geoclim_type: madaclim_info.select_geoclim_type_layers(geoclim_type)["layer_number"].to_list() for geoclim_type in geoclim_types}

        clim_raster_layers_to_sample = [layer_num for layer_num in layers_numbers if layer_num in geoclim_layer_ranges["clim"]]
        env_raster_layers_to_sample = [layer_num for layer_num in layers_numbers if layer_num in geoclim_layer_ranges["env"]]

        # # Validate if mada_geom_point attribute is of Point geom and not empty
        if not isinstance(self.mada_geom_point, shapely.geometry.point.Point):
            raise TypeError("The 'mada_geom_point' attribute must be a shapely.geometry.point.Point object.")
        
        if self.mada_geom_point.is_empty:
            raise ValueError("The 'mada_geom_point' object cannot be empty.")
        
        # Sample climate and env raster on demand
        sampled_data = {}
        nodata_layers = []
        
        print("\n" + "#" * 40 + f" \033[1mExtracting data for: {self.specimen_id}\033[0m " +"#" * 40)
        start_time = time.perf_counter()
        
        if clim_raster_layers_to_sample:
            total_clim_layers = len(clim_raster_layers_to_sample)
            
            with rasterio.open(clim_raster_path or default_clim_raster_path) as clim_raster:
                # Initialize reference and container to check for layers with nodata values
                nodata_clim = clim_raster.nodata
                
                # Status bar to display when sampling the raster
                print(f"\nSampling {total_clim_layers} layer(s) from {clim_raster.name.split('/')[-1]} (geoclim_type={geoclim_types[0]})...")
                with tqdm(
                    total=total_clim_layers, 
                    unit="layer",
                    bar_format="{desc} {percentage:.0f}%|{bar}| layer {n_fmt}/{total_fmt} [Time remaining: {remaining}]",
                ) as pbar :
                    
                    # Sample selected layers according to the coordinate
                    for layer_num in clim_raster_layers_to_sample:
                        # Get layer info for pbar display and label key for sampled_data
                        layer_name = madaclim_info.fetch_specific_layers(layers_labels=layer_num, as_descriptive_labels=True, return_list=True)[0]
                        layer_description_display = layer_name.split('_')[-1]
                        pbar.set_description(f"Extracting layer {layer_num}: {layer_description_display}")
                        pbar.update()
                        
                        # Sample using the self.mada_geom_point attributes coordinates for the current layer
                        band = madaclim_info.get_band_from_layer_number(layer_num, geoclim_types[0])
                        data = list(clim_raster.sample([(self.mada_geom_point.x, self.mada_geom_point.y)], indexes=band))[0]
                        
                        # Save extracted data with specified layer info/name
                        if not layer_info:
                            layer_name = f"layer_{layer_num}"
                        sampled_data[layer_name] = data[0]

                        if data[0] == nodata_clim:    # Save layers where nodata at specimen location
                            nodata_layers.append(layer_name)
                                    
        if env_raster_layers_to_sample:
            total_env_layers = len(env_raster_layers_to_sample)
            
            with rasterio.open(env_raster_path or default_env_raster_path) as env_raster:
                # Initialize reference and container to check for layers with nodata values
                nodata_env = env_raster.nodata
                
                # Status bar to display when sampling the raster
                print(f"\nSampling {total_env_layers} layer(s) from {env_raster.name.split('/')[-1]} (geoclim_type={geoclim_types[1]})...")
                with tqdm(
                    total=total_env_layers, 
                    unit="layer",
                    bar_format="{desc} {percentage:.0f}%|{bar}| layer {n_fmt}/{total_fmt} [Time remaining: {remaining}]",
                ) as pbar :
                    
                    # Sample selected layers according to the coordinate
                    for layer_num in env_raster_layers_to_sample:
                        # Get layer info for pbar display and label key for sampled_data
                        layer_name = madaclim_info.fetch_specific_layers(layers_labels=layer_num, as_descriptive_labels=True, return_list=True)[0]
                        layer_description_display = layer_name.split('_')[-1]
                        pbar.set_description(f"Extracting layer {layer_num}: {layer_description_display}")
                        pbar.update()
                        
                        # Sample using the self.mada_geom_point attributes coordinates for the current layer
                        band = madaclim_info.get_band_from_layer_number(layer_num, geoclim_types[1])
                        data = list(env_raster.sample([(self.mada_geom_point.x, self.mada_geom_point.y)], indexes=band))[0]
                        
                        # Save extracted data with specified layer info/name
                        if not layer_info:
                            layer_name = f"layer_{layer_num}"
                        sampled_data[layer_name] = data[0]

                        if data[0] == nodata_clim:    # Save layers where nodata at specimen location
                            nodata_layers.append(layer_name)
                                    

        if len(nodata_layers) > 0:    # No raising exception, just warning print
            print(f"BEWARE! {len(nodata_layers)} layer(s) contain a nodata value at the specimen location")

        end_time = time.perf_counter()
        elapsed_time = end_time - start_time    # Total raster sampling time
        print(f"\nFinished raster sampling operation in {elapsed_time:.2f} seconds.\n")
        
        if return_nodata_layers:
            return sampled_data, nodata_layers
        
        return sampled_data
    
    def _get_additional_attributes(self) -> dict:
        """
        Get additional attributes of the instance.

        This method retrieves any attributes that were added to the instance 
        after initialization, i.e., attributes that are not part of the base attributes.

        Note: This is a private method, indicated by the underscore prefix. 
        It's intended for internal use within the class, not for use by external code.

        Returns:
            dict: A dictionary containing the names and values of the additional attributes.
                The dictionary keys are attribute names and the dictionary values are attribute values.
                If there are no additional attributes, an empty dictionary is returned.
        """
        current_attributes = {k:v for k,v in self.__dict__.items() if k != "_MadaclimPoint__base_attr"}    # Remove the recursiveness of __base_attr
        additional_attributes = {k:v for k,v in current_attributes.items() if k not in self.__base_attr.keys()}
        return additional_attributes

    
    def _construct_point(self, latitude: float, longitude: float, source_crs: pyproj.crs.crs.CRS)-> shapely.geometry.point.Point:
        """
        Construct a Shapely Point object in the given source_crs, reprojecting the point if necessary.

        Args:
            latitude (float): The latitude of the point.
            longitude (float): The longitude of the point.
            source_crs (pyproj.crs.crs.CRS): The coordinate reference system of the point.

        Returns:
            shapely.geometry.point.Point: A Shapely Point object representing the point in the appropriate CRS.

        Raises:
            ValueError: If the climate and environmental rasters have different CRS, which is unexpected.
        """
        # Get the crs' from both rasters of the Madaclim db
        madaclim_info = MadaclimLayers()
        madaclim_crs = madaclim_info.clim_crs if madaclim_info.clim_crs == madaclim_info.env_crs else None

        # Sanity check for crs
        if madaclim_crs is None:
            raise ValueError("Beware, the clim and env rasters have different projections/CRS which is unexpected.")
        
        
        # Create a Point object in the source CRS
        point = Point(longitude, latitude)
        if source_crs == madaclim_crs:
            return point

        # Create a transformer object for converting coordinates between the two CRS
        else:
            transformer = Transformer.from_crs(source_crs, madaclim_crs, always_xy=True)
            x, y = transformer.transform(point.x, point.y)
            reprojected_point = Point(x, y)
            return reprojected_point
    
    def _update_mada_geom_point(self):
        """
        Update the mada_geom_point attribute by reconstructing the point with the current latitude, longitude, and source_crs.
        """
        if hasattr(self, '_latitude') and hasattr(self, '_longitude'):
            self.mada_geom_point = self._construct_point(
                latitude=self.latitude,
                longitude=self.longitude,
                source_crs=self.source_crs
            )
    
    
class MadaclimCollection:
    #TODO DOCSTRINGS CLS
    def __init__(self, madaclim_points: Optional[Union[MadaclimPoint, List[MadaclimPoint]]]=None) -> None:
        """Instantiate a collection of MadaclimPoint objects. By default, the MadaclimCollection is empty.
        It will populate the collection with a single instance or a list of MadaclimPoint instances by calling the add_points method with the given madaclim_points.

        Args:
            madaclim_points (Optional[Union[MadaclimPoint, List[MadaclimPoint]]], optional): A single MadaclimPoint object or a list of MadaclimPoint objects to be added to the MadaclimCollection. Initialize an empty MadaclimCollection by default (None).
        Examples:
            >>> from coffeaphylogeo.geoclim.raster_manipulation import MadaclimPoint, MadaclimCollection
            >>> specimen_1 = MadaclimPoint(specimen_id="spe1", latitude=-23.574583, longitude=46.419806, source_crs="epsg:4326")
            >>> specimen_2 = MadaclimPoint(specimen_id="spe2", latitude=-2622095.832726487, longitude=5048512.906023483, source_crs=3857)
            
            >>> # Add points to the collection when constructing
            >>> collection = MadaclimCollection([specimen_1, specimen_2])
            >>> print(collection)
            MadaclimCollection = [
                MadaclimPoint(specimen_id=spe1, mada_geom_point=POINT (644890.8921103649 7392153.658976035)),
                MadaclimPoint(specimen_id=spe2, mada_geom_point=POINT (536050.6239664567 7465523.013290589))
            ]
            >>> # You can also initiliaze an empty collection
            >>> collection = MadaclimCollection()
            >>> print(collection)
            No MadaclimPoint inside the collection.

            >>> # Add a single MadaclimPoint
            >>> collection.add_points(specimen_1)
            >>> print(collection)
            MadaclimCollection = [
                MadaclimPoint(specimen_id=spe1, mada_geom_point=POINT (644890.8921103649 7392153.658976035))
            ]
            
            >>> # Populate the collection from a csv or df
            >>> collection = MadaclimCollection.populate_from_csv("some_samples.csv")
            Warning! No source_crs column in the csv. Using the default value of EPSG:4326...
            Creating MadaclimPoint(specimen_id=sample_A...)
            \Creating MadaclimPoint(specimen_id=sample_B...)
            Creating MadaclimPoint(specimen_id=sample_C...)
            Creating MadaclimPoint(specimen_id=sample_D...)
            Creating MadaclimPoint(specimen_id=sample_E...)
            Created new MadaclimCollection with 5 samples.
            >>> print(collection)
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_A, mada_geom_point=POINT (837072.9150244407 7903496.320897499)),
                MadaclimPoint(specimen_id=sample_B, mada_geom_point=POINT (695186.2170220022 8197477.647690434)),
                MadaclimPoint(specimen_id=sample_C, mada_geom_point=POINT (761613.8281386737 7651088.106452912)),
                MadaclimPoint(specimen_id=sample_D, mada_geom_point=POINT (955230.600222457 8005985.896187438)),
                MadaclimPoint(specimen_id=sample_E, mada_geom_point=POINT (757247.2175273325 7618631.528869408))
            ]

            >>> # MadaclimPoints are stored in the .all_points attributes in a list
            >>> collection.all_points[0]
            MadaclimPoint(
                specimen_id = sample_A,
                source_crs = 4326,
                latitude = -18.9333,
                longitude = 48.2,
                mada_geom_point = POINT (837072.9150244407 7903496.320897499)
            )
        """
        self.__all_points = []
        if madaclim_points:
            self.add_points(madaclim_points)
        self.__sampled_rasters_data = None
        self.__nodata_layers = None

    @property
    def all_points(self) -> list:
        """Get the all_points attribute.

        Returns:
            list: A list of all the MadaclimPoint objects in the MadaclimCollection.
        Examples:
            >>> # MadaclimPoints are stored in the .all_points attributes in a list
            >>> sample_A = MadaclimPoint()
            >>> collection = MadaclimCollection(sample_A)
            
            >>> collection.all_points[0]
            MadaclimPoint(
                longitude = 48.2,
                specimen_id = sample_A,
                latitude = -18.9333,
                source_crs = 4326,
                mada_geom_point = POINT (837072.9150244407 7903496.320897499)
            )
        """
        return self.__all_points
    
    @property
    def sampled_rasters_data(self) -> Dict[str, Dict[str, float]]:
        """Get the sampled_rasters_data attribute.

        This attribute is a nested dictionary. The outer dictionary uses the MadaclimPoint.specimen_id as keys. 
        The corresponding value for each key is another dictionary, which uses layer_names as keys and sampled values from rasters as values.

        Returns:
            Dict[str, Dict[str, float]]: A dictionary with MadaclimPoint.specimen_id as keys and a dictionary of layer_names (str) and sampled values (float) as values or "nodata_layers" (str) and names of the layers with nodata values (list) as values. 
        """
        return self.__sampled_rasters_data

    @property
    def nodata_layers(self) -> Dict[str, Union[str, List[str]]]:
        """Get the nodata_layers attribute.

        This attribute is a dictionary that contains the MadaclimPoint.specimen_id as keys and the values as the 'nodata_layers' as str or list of str.
        
        Returns:
            Dict[str, Union[str, List[str]]]: A dictionary with MadaclimPoint.specimen_id as keys and values of str or list of str of the layers_name with nodata values.
        """
        return self.__nodata_layers
    
    def __str__(self) -> str:
        if len(self.__all_points) == 0:
            return "No MadaclimPoint inside the collection."
        else:
            all_points_short = [
                f"MadaclimPoint(specimen_id={point.specimen_id}, mada_geom_point={point.mada_geom_point})" for point in self.all_points
            ]
            return "MadaclimCollection = [\n" + "\t" + ",\n\t".join(all_points_short) + "\n]"
        
    def __repr__(self) -> str:
        return self.__str__()
    
    @classmethod
    def populate_from_csv(cls, csv_file: Union[str, pathlib.Path]) -> "MadaclimCollection":
        """Creates a new MadaclimCollection from a CSV file.
        
        Each row of the CSV file should represent a MadaclimPoint. The CSV file
        must have columns that correspond to the arguments of the MadaclimPoint
        constructor. If a 'source_crs' column is not provided, the method uses
        the default CRS value.

        Args:
            csv_file (Union[str, pathlib.Path]): The path to the CSV file.

        Returns:
            MadaclimCollection: A new MadaclimCollection instance with MadaclimPoint
            objects created from the rows of the CSV file.

        Raises:
            TypeError: If 'csv_file' is not a str or pathlib.Path object.
            FileNotFoundError: If the file specified by 'csv_file' does not exist.
            ValueError: If the CSV file headers are missing required arguments for
            constructing MadaclimPoint objects.
        Examples:
            # csv headers must contain the 3 required positional args for MadaclimPoint (specimen_id, lat, lon)
            # When no source_crs header is found, defaults to EPSG:4326

            # CSV example #1: "some_samples.csv"
            specimen_id,latitude,longitude
            sample_A,-18.9333,48.2
            sample_B,-16.295741,46.826763
            sample_C,-21.223,47.5204
            sample_D,-17.9869,49.2966
            sample_E,-21.5166,47.4833

            >>> collection = MadaclimCollection.populate_from_csv("some_samples.csv")
            Warning! No source_crs column in the csv. Using the default value of EPSG:4326...
            Creating MadaclimPoint(specimen_id=sample_A...)
            \Creating MadaclimPoint(specimen_id=sample_B...)
            Creating MadaclimPoint(specimen_id=sample_C...)
            Creating MadaclimPoint(specimen_id=sample_D...)
            Creating MadaclimPoint(specimen_id=sample_E...)
            Created new MadaclimCollection with 5 samples.
            >>> print(collection)
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_A, mada_geom_point=POINT (837072.9150244407 7903496.320897499)),
                MadaclimPoint(specimen_id=sample_B, mada_geom_point=POINT (695186.2170220022 8197477.647690434)),
                MadaclimPoint(specimen_id=sample_C, mada_geom_point=POINT (761613.8281386737 7651088.106452912)),
                MadaclimPoint(specimen_id=sample_D, mada_geom_point=POINT (955230.600222457 8005985.896187438)),
                MadaclimPoint(specimen_id=sample_E, mada_geom_point=POINT (757247.2175273325 7618631.528869408))
            ]

            # Can accept other non-required data for MadaclimPoint instantiation
            
            # CSV example #2: "other_samples.csv"
            specimen_id,latitude,longitude,source_crs,has_sequencing,specie
            sample_F,-19.9333,47.2,4326,True,bojeri
            sample_G,-18.295741,45.826763,4326,False,periwinkle
            sample_H,-21.223,44.5204,4326,False,spectabilis
            
            >>> other_collection = MadaclimCollection.populate_from_csv("other_samples.csv")
            Creating MadaclimPoint(specimen_id=sample_F...)
            Creating MadaclimPoint(specimen_id=sample_G...)
            Creating MadaclimPoint(specimen_id=sample_H...)
            Created new MadaclimCollection with 3 samples.
            >>> for point in other_collection.all_points:
            ...     print(point)
            ... 
            MadaclimPoint(
                specimen_id = sample_F,
                source_crs = 4326,
                latitude = -19.9333,
                longitude = 47.2,
                mada_geom_point = POINT (730272.0056458472 7794391.966030249),
                has_sequencing = True,
                specie = bojeri
            )
            MadaclimPoint(
                specimen_id = sample_G,
                source_crs = 4326,
                latitude = -18.295741,
                longitude = 45.826763,
                mada_geom_point = POINT (587378.6907481698 7976896.406900212),
                has_sequencing = False,
                specie = periwinkle
            )
            MadaclimPoint(
                specimen_id = sample_H,
                source_crs = 4326,
                latitude = -21.223,
                longitude = 44.5204,
                mada_geom_point = POINT (450229.7195355138 7653096.609718417),
                has_sequencing = False,
                specie = spectabilis
            )

        """
        # Convert str to pathlib.Path
        if isinstance(csv_file, str):
            csv_file = Path(csv_file)
        
        # Type and IO validation
        if not isinstance(csv_file, pathlib.Path):
            raise TypeError("'csv_file' must be a valid pathlib.Path object.")

        if not csv_file.is_file():
            raise FileNotFoundError(f"Could not find {csv_file}")
        
        # Get the required + default args to use to construct the MadaclimPoint objects
        madapoint_required_args, madapoint_default_crs_arg = MadaclimPoint.get_args_names()
        madapoint_default_crs_val = MadaclimPoint.get_default_source_crs()

        # Populate the collection from the input csv
        with open(csv_file, newline="") as f:
            csv_data = csv.DictReader(f)
            
            # Check if all required_args in csv
            col_names = csv_data.fieldnames
            missing_args = [req_arg for req_arg in madapoint_required_args if req_arg not in col_names]
            if len(missing_args) > 0:
                raise ValueError(f"csv file headers are missing the following required args to construct the MadaclimPoint objects:\n{missing_args}")
            
            # Warn for default EPSG if not provided
            if madapoint_default_crs_arg not in col_names:
                print(f"Warning! No {madapoint_default_crs_arg} column in the csv. Using the default value of EPSG:{madapoint_default_crs_val}...")

            # Initialize MadaclimPoint instances container to fill from csv
            points = []

            for row in csv_data:
                # If source_crs not in row use default val
                if madapoint_default_crs_arg not in row or not row[madapoint_default_crs_arg]:
                    row[madapoint_default_crs_arg] = madapoint_default_crs_val
                
                # Create a MadaclimPoint using the values of the row
                print(f"Creating MadaclimPoint(specimen_id={row['specimen_id']}...)")
                point = MadaclimPoint(**row)
                points.append(point)
        
        new_collection = cls(points)
        print(f"Created new MadaclimCollection with {len(points)} samples.")
        return new_collection
    
    @classmethod
    def populate_from_df(cls, df: pd.DataFrame) -> "MadaclimCollection":
        """
        Class method to populate a MadaclimCollection from a pandas DataFrame.

        This method takes a DataFrame where each row represents a MadaclimPoint and its 
        attributes. If the 'source_crs' column is not provided in the DataFrame, the default 
        CRS will be used. 

        Args:
            df (pd.DataFrame): DataFrame where each row represents a MadaclimPoint. Expected
                            columns are the same as the required arguments for the 
                            MadaclimPoint constructor.

        Returns:
            MadaclimCollection: A new MadaclimCollection instance populated with MadaclimPoints
                                created from the DataFrame.

        Raises:
            TypeError: If 'df' is not a pd.DataFrame.
            ValueError: If the DataFrame is missing any of the required arguments to construct 
                        a MadaclimPoint.
        Examples:
            >>> import pandas as pd

            >>> sample_df
            specimen_id    latitude  longitude
            0    sample_W  -16.295741  46.826763
            1    sample_X    -17.9869    49.2966
            2    sample_Y    -18.9333    48.2166
            3    sample_Z      -13.28      49.95

            >>> collection = MadaclimCollection.populate_from_df(sample_df)
            Warning! No source_crs column in the df. Using the default value of EPSG:4326...
            Creating MadaclimPoint(specimen_id=sample_W...)
            Creating MadaclimPoint(specimen_id=sample_X...)
            Creating MadaclimPoint(specimen_id=sample_Y...)
            Creating MadaclimPoint(specimen_id=sample_Z...)
            Created new MadaclimCollection with 4 samples.
            >>> collection.all_points[0]
            MadaclimPoint(
                specimen_id = sample_W,
                source_crs = 4326,
                latitude = -16.295741,
                longitude = 46.826763,
                mada_geom_point = POINT (695186.2170220022 8197477.647690434)
            )

        """
        if not isinstance(df, pd.DataFrame):
            raise TypeError("'df' is not a pd.DataFrame.")
        
        # Get the required + default args to use to construct the MadaclimPoint objects
        madapoint_required_args, madapoint_default_crs_arg = MadaclimPoint.get_args_names()
        madapoint_default_crs_val = MadaclimPoint.get_default_source_crs()

        # Check for required args present in df
        missing_args = [req_arg for req_arg in madapoint_required_args if req_arg not in df.columns]
        if len(missing_args) > 0:
            raise ValueError(f"df is missing the following required args to construct the MadaclimPoint objects:\n{missing_args}")
        
        # Warn for default EPSG if not provided
        if madapoint_default_crs_arg not in df.columns:
            print(f"Warning! No {madapoint_default_crs_arg} column in the df. Using the default value of EPSG:{madapoint_default_crs_val}...")

        # Initialize MadaclimPoint instances container to fill from df
        points = []

        for _, row in df.iterrows():
            # If source_crs not in row use default val
            if madapoint_default_crs_arg not in row or not row[madapoint_default_crs_arg]:
                row[madapoint_default_crs_arg] = madapoint_default_crs_val

            print(f"Creating MadaclimPoint(specimen_id={row['specimen_id']}...)")
            point = MadaclimPoint(**row)
            points.append(point)
        
        new_collection = cls(points)
        print(f"Created new MadaclimCollection with {len(points)} samples.")
        return new_collection

    def add_points(self, madaclim_points: Union[MadaclimPoint, List[MadaclimPoint]]) -> None:
        """
        Adds one or more MadaclimPoint objects to the MadaclimCollection.

        Args:
            madaclim_points (Union[MadaclimPoint, List[MadaclimPoint]]): A single MadaclimPoint object or a list of MadaclimPoint objects to be added to the MadaclimCollection.

        Raises:
            TypeError: If the input is not a MadaclimPoint object or a list of MadaclimPoint objects.
            ValueError: If the input MadaclimPoint(s) is/are already in the MadaclimCollection or if their specimen_id(s) are not unique.
        Examples:
            >>> from coffeaphylogeo.geoclim.raster_manipulation import MadaclimPoint, MadaclimCollection

            >>> specimen_1 = MadaclimPoint(specimen_id="spe1", latitude=-23.574583, longitude=46.419806, source_crs="epsg:4326")
            >>> collection = MadaclimCollection()
            >>> print(collection)
            No MadaclimPoint inside the collection.

            >>> # Add a single point
            >>> collection.add_points(specimen_1)
            >>> print(collection)
            MadaclimCollection = [
                MadaclimPoint(specimen_id=spe1, mada_geom_point=POINT (644890.8921103649 7392153.658976035))
            ]

            >>> # Add multiple points
            >>> specimen_2 = MadaclimPoint(specimen_id="spe2", latitude=-20.138470, longitude=46.054688, family="Rubiaceae", has_sequencing=True, num_samples=1)            
            >>> other_collection.add_points([specimen_1, specimen_2])
            >>> print(other_collection)
            MadaclimCollection = [
                MadaclimPoint(specimen_id=spe1, mada_geom_point=POINT (644890.8921103649 7392153.658976035)),
                MadaclimPoint(specimen_id=spe2, mada_geom_point=POINT (610233.867750987 7772846.143786541))
            ]

            # You cannot add duplicates, each point must have a unique specimen_id attribute.
            >>> other_collection.add_points(specimen_1)
            Traceback (most recent call last):
            File "<stdin>", line 1, in <module>
            File ".../coffeaPhyloGeo/src/coffeaphylogeo/geoclim/raster_manipulation.py", line 1013, in add_points
                MadaclimPoint(specimen_id=spe2, mada_geom_point=POINT (610233.867750987 7772846.143786541))
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ValueError: MadaclimPoint(
                specimen_id = spe1,
                source_crs = 4326,
                latitude = -23.574583,
                longitude = 46.419806,
                mada_geom_point = POINT (644890.8921103649 7392153.658976035)
            ) is already in the current MadaclimCollection instance.

        """
        # Add multiple MadaclimPoint objects
        if isinstance(madaclim_points, list):
            for point in madaclim_points:
                
                if not isinstance(point, MadaclimPoint):
                    raise TypeError(f"{point} is not a MadaclimPoint object. Accepted types are a single MadaclimPoint and a list of MadaclimPoint objects.")
                
                if point in self.__all_points:
                    raise ValueError(f"{point} is already in the current MadaclimCollection instance.")
                
                if point.specimen_id in [mada_pt.specimen_id for mada_pt in self.__all_points]:    # specimen_id unique id validation
                    raise ValueError(f"specimen_id={point.specimen_id} already exists inside current MadaclimCollection. Every MadaclimPoint must have a unique id.")
                
                self.__all_points.append(point) 
        else:
            # Add single MadaclimPoint object
            if not isinstance(madaclim_points, MadaclimPoint):
                raise TypeError("The madaclim_point to add is not a MadaclimPoint object.")
            
            if madaclim_points in self.__all_points:
                raise ValueError(f"{madaclim_points} is already in the current MadaclimCollection instance.")
            
            if madaclim_points.specimen_id in [mada_pt.specimen_id for mada_pt in self.__all_points]:    # specimen_id unique id validation
                    raise ValueError(f"specimen_id={madaclim_points.specimen_id} already exists inside current MadaclimCollection. Every MadaclimPoint must have a unique id.")
                
            
            self.__all_points.append(madaclim_points)

    def remove_points(self, *, madaclim_points: Optional[Union[MadaclimPoint, List[MadaclimPoint]]]=None, indices: Optional[Union[int, List[int]]]=None, clear:bool=False) -> None:
        """Removes MadaclimPoint objects from the MadaclimCollection based on specified criteria.

        This method allows removing MadaclimPoint objects from the collection by providing either
        MadaclimPoint instance(s), index/indices, or by clearing the whole collection.
        
        Args:
            madaclim_points (Optional[Union[MadaclimPoint, List[MadaclimPoint]]], optional): A single MadaclimPoint
                object or a list of MadaclimPoint objects to be removed from the collection. Defaults to None.
            indices (Optional[Union[int, List[int]]], optional): A single index or a list of indices of the MadaclimPoint
                objects to be removed from the collection. Defaults to None.
            clear (bool, optional): If set to True, removes all MadaclimPoint objects from the collection. When using
                this option, 'madaclim_points' and 'indices' must not be provided. Defaults to False.

        Raises:
            ValueError: If the MadaclimCollection is empty or if none of the input options are provided.
            ValueError: If 'madaclim_points' and 'indices' are both provided.
            ValueError: If 'clear' is set to True and either 'madaclim_points' or 'indices' are provided.
            TypeError: If an invalid type is provided for 'madaclim_points' or 'indices'.
            ValueError: If a provided MadaclimPoint object is not in the collection or if an index is out of bounds.
            IndexError: If an index is out of range.

        Examples:
            >>> import pandas as pd
            >>> sample_df
              specimen_id   latitude  longitude source_crs
            0    sample_W -16.295741  46.826763  EPSG:4326
            1    sample_X -17.986900  49.296600  EPSG:4326
            2    sample_Y -18.933300  48.216600  EPSG:4326
            3    sample_Z -13.280000  49.950000  EPSG:4326

            >>> collection = MadaclimCollection.populate_from_df(sample_df)
            Creating MadaclimPoint(specimen_id=sample_W...)
            Creating MadaclimPoint(specimen_id=sample_X...)
            Creating MadaclimPoint(specimen_id=sample_Y...)
            Creating MadaclimPoint(specimen_id=sample_Z...)
            Created new MadaclimCollection with 4 samples.

            # Remove points by passing in the MadaclimPoint objects, the index or the specimen_id to remove from the Collection.
            
            >>> # From MadaclimPoint instances
            >>> sample_W = collection.all_points[0]
            >>> sample_W
            MadaclimPoint(
                specimen_id = sample_W,
                source_crs = 4326,
                latitude = -16.295741,
                longitude = 46.826763,
                mada_geom_point = POINT (695186.2170220022 8197477.647690434)
            )

            >>> # You must specify madaclim_points keyword arg
            >>> collection.remove_points(sample_W)
            Traceback (most recent call last):
            File "<stdin>", line 1, in <module>
            TypeError: MadaclimCollection.remove_points() takes 1 positional argument but 2 were given
            >>> collection.remove_points(madaclim_points=sample_W)
            >>> >>> print(collection)
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_X, mada_geom_point=POINT (955230.600222457 8005985.896187438)),
                MadaclimPoint(specimen_id=sample_Y, mada_geom_point=POINT (838822.9378705097 7903464.491492193)),
                MadaclimPoint(specimen_id=sample_Z, mada_geom_point=POINT (1036778.1471182993 8526563.721996231))
            ]

            >>> # Using the position index of the instance
            >>> collection = MadaclimCollection.populate_from_df(sample_df)
            Creating MadaclimPoint(specimen_id=sample_W...)
            Creating MadaclimPoint(specimen_id=sample_X...)
            Creating MadaclimPoint(specimen_id=sample_Y...)
            Creating MadaclimPoint(specimen_id=sample_Z...)
            Created new MadaclimCollection with 4 samples.
            >>> print(collection)
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_W, mada_geom_point=POINT (695186.2170220022 8197477.647690434)),
                MadaclimPoint(specimen_id=sample_X, mada_geom_point=POINT (955230.600222457 8005985.896187438)),
                MadaclimPoint(specimen_id=sample_Y, mada_geom_point=POINT (838822.9378705097 7903464.491492193)),
                MadaclimPoint(specimen_id=sample_Z, mada_geom_point=POINT (1036778.1471182993 8526563.721996231))
            ]

            >>> collection.remove_points(indices=-1)    # Removes last point of the collection
            >>> print(collection)
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_W, mada_geom_point=POINT (695186.2170220022 8197477.647690434)),
                MadaclimPoint(specimen_id=sample_X, mada_geom_point=POINT (955230.600222457 8005985.896187438)),
                MadaclimPoint(specimen_id=sample_Y, mada_geom_point=POINT (838822.9378705097 7903464.491492193))
            ]

            >>> # Using the specimen.id attribute
            >>> collection.remove_points(madaclim_points="sample_Y")
            >>> print(collection)
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_W, mada_geom_point=POINT (695186.2170220022 8197477.647690434)),
                MadaclimPoint(specimen_id=sample_X, mada_geom_point=POINT (955230.600222457 8005985.896187438)),
                MadaclimPoint(specimen_id=sample_Z, mada_geom_point=POINT (1036778.1471182993 8526563.721996231))
            ]

            # We can also use a list to remove multiple points at the same time.
            # A list of str or MadaclimPoint or mixed types are accepted for the madaclim_points argument.
            >>> sample_w = collection.all_points[0]
            >>> to_remove = [sample_w, "sample_X"]
            >>> collection.remove_points(madaclim_points=to_remove)
            >>> print(collection)
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_Y, mada_geom_point=POINT (838822.9378705097 7903464.491492193)),
                MadaclimPoint(specimen_id=sample_Z, mada_geom_point=POINT (1036778.1471182993 8526563.721996231))
            ]

            # Or pass in a list of indices to the indices argument.
            >>> collection.remove_points(indices=[0, -1])    # Remove first and last point
            >>> print(collection)
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_X, mada_geom_point=POINT (955230.600222457 8005985.896187438)),
                MadaclimPoint(specimen_id=sample_Y, mada_geom_point=POINT (838822.9378705097 7903464.491492193))
            ]

            # Finaly we can clear the collection of all instances.
            >>> collection.remove_points(clear=True)
            No MadaclimPoint inside the collection.
        """
        # Handle empty MadaclimCollection
        if not len(self.__all_points) > 0:
            raise ValueError("No points to delete since the MadaclimCollection is empty.")
        
        # Drop all points
        if clear:
            if madaclim_points is not None or indices is not None:
                raise ValueError("When using 'clear', do not provide 'madaclim_points' or 'indices'.")
            else:
                self.__all_points.clear()
                return
        if madaclim_points is not None and indices is not None:
            raise ValueError("Either provide 'madaclim_points' or 'indices', not both.")
        
        if madaclim_points is None and indices is None:
            raise ValueError("At least one of madaclim_points or indices must be provided.")
        
        # Remove single/multiple points by MadaclimPoint instance(s)
        if madaclim_points is not None:
            # Containers for objects (MadaclimPoint) or specimen_id (str) to remove
            madaclim_points_objects = []
            madaclim_points_ids = []
            
            # Multiple objects removal (list arg)
            if isinstance(madaclim_points, list):
                for point_to_rm in madaclim_points:
                    if not isinstance(point_to_rm, (MadaclimPoint, str)):
                        raise TypeError(f"{point_to_rm} is not a MadaclimPoint object or str. Accepted types are a single MadaclimPoint, a str, a list of MadaclimPoint objects or a list of str.")
                    
                    if isinstance(point_to_rm, MadaclimPoint):
                        if point_to_rm not in self.__all_points:
                            raise ValueError(f"{point_to_rm} does not exists in the current MadaclimCollection instance.")
                        madaclim_points_objects.append(point_to_rm)
                    else:
                        if point_to_rm not in [point.specimen_id for point in self.__all_points]:
                            raise ValueError(f"{point_to_rm} is not a valid MadaclimPoint.specimen_id for all points of the collection.")
                        madaclim_points_ids.append(point_to_rm)

                self.__all_points = [point for point in self.__all_points if point not in madaclim_points_objects and point.specimen_id not in madaclim_points_ids]

            # Single object removal
            else:    
                if not isinstance(madaclim_points, (MadaclimPoint, str)):
                    raise TypeError(f"{madaclim_points} is not a MadaclimPoint object or str. Accepted types are a single MadaclimPoint, a str, a list of MadaclimPoint objects or a list of str.")
                
                if isinstance(madaclim_points, MadaclimPoint):
                    if madaclim_points not in self.__all_points:
                        raise ValueError(f"{madaclim_points} does not exists in the current MadaclimCollection instance.")
                    madaclim_points_objects.append(madaclim_points)
                else:
                    if madaclim_points not in [point.specimen_id for point in self.__all_points]:
                        raise ValueError(f"{madaclim_points} is not a valid MadaclimPoint.specimen_id for all points of the collection.")
                    madaclim_points_ids.append(madaclim_points)

                self.__all_points = [point for point in self.__all_points if point not in madaclim_points_objects and point.specimen_id not in madaclim_points_ids]

        
        
        # Remove multiple points from an indices list
        if indices is not None:
            if isinstance(indices, list):
                try:
                    indices = [int(index) if index >= 0 else len(self.__all_points) + int(index) for index in indices]
                except:
                    raise TypeError("Indices should be integers.")
                for index in indices:
                    if index < 0 or index >= len(self.__all_points):
                        raise ValueError(f"Index {index} is out of bounds.")
                self.__all_points = [point for i, point in enumerate(self.__all_points) if i not in indices]

            else:    # Remove single point
                try:
                    index = int(indices) if indices >= 0 else len(self.__all_points) + int(indices)
                except:
                    raise TypeError("Single index must be an integer.")
                
                if index < 0 or index >= len(self.__all_points):
                    raise IndexError("Index out of range.")
                
                else:
                    self.__all_points.pop(index)


    def sample_from_rasters(
            self,
            layers_to_sample: Union[int, str, List[Union[int, str]]]="all", 
            layer_info: bool=False,
            return_nodata_layers: bool=False,
            clim_raster_path: Optional[pathlib.Path]=None, 
            env_raster_path: Optional[pathlib.Path]=None
        ) -> Union[Dict[str, Dict[str, float]], Tuple[Dict[str, Dict[str, float]], Optional[Dict[str, List[str]]]]]:
        """
        Sample the given raster layers for all points in the MadaclimCollection.

        Args:
            layers_to_sample (Union[int, str, List[Union[int, str]]], optional): The raster layers to sample. Can be an integer (layer index), a string (layer name), or a list of integers or strings. If 'all', all layers are sampled. Defaults to 'all'.
            layer_info (bool, optional): If True, return detailed information about each layer. Defaults to False.
            return_nodata_layers (bool, optional): If True, return layers where the sampled value is nodata. Defaults to False.
            clim_raster_path (Optional[pathlib.Path], optional): Path to the climate raster file. If not provided, the default path is used. Defaults to None.
            env_raster_path (Optional[pathlib.Path], optional): Path to the environmental raster file. If not provided, the default path is used. Defaults to None.

        Returns:
            Union[Dict[str, Dict[str, float]], Tuple[Dict[str, Dict[str, float]], Optional[Dict[str, List[str]]]]]: 
                If return_nodata_layers is False: Returns a dictionary where the keys are specimen_ids and the values are dictionaries of sampled data.
                If return_nodata_layers is True: Returns a tuple where the first element is the same dictionary as when return_nodata_layers is False, and the second element is a dictionary where the keys are specimen_ids and the values are lists of layers where the sampled value is nodata.
       
            
        Raises:
            ValueError: If the MadaclimCollection doesn't contain any MadaclimPoints.

        Notes:
            This method also updates the 'sampled_rasters_data' and 'nodata_layers' attributes of the MadaclimCollection instance.

        Examples:


        """
        
        if not len(self.__all_points) > 0:
            raise ValueError("Not MadaclimPoint to sample from in the Collection.")
        
        # # Initialize containers for sampled data
        sampled_data = {}
        nodata_layers = {}

        # Sample rasters for whole collection
        for point in self.__all_points:
            point_data = point.sample_from_rasters(
                layers_to_sample=layers_to_sample,
                layer_info=layer_info,
                return_nodata_layers=return_nodata_layers,
                clim_raster_path=clim_raster_path,
                env_raster_path=env_raster_path
            )

            if return_nodata_layers:
                sampled_data_point, nodata_layers_point = point_data    # Tuple from single sample_from_rasters output
                
                sampled_data[point.specimen_id] = sampled_data_point    # Save raster data
    
                # Save layers name only when val is nodata
                if len(nodata_layers_point) > 0:
                    nodata_layers[point.specimen_id] = nodata_layers_point
            else:
                sampled_data[point.specimen_id] = point_data
        
        # Update instance attributes
        self.__sampled_rasters_data = sampled_data
        self.__nodata_layers = nodata_layers if len(nodata_layers) > 0 else None

        if return_nodata_layers:
            return sampled_data, nodata_layers if len(nodata_layers) > 0 else None
        else:
            return sampled_data
        