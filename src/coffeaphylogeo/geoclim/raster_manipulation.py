import pathlib
from typing import Optional, Union, List
import time
from tqdm import tqdm

from coffeaphylogeo.definitions import Definitions
from coffeaphylogeo.geoclim.madaclim_info import MadaclimLayers

import rasterio
import pyproj
import shapely
from pyproj import Transformer
from shapely import Point

# Default dir and path for rasters
defs = Definitions()

climate_dir = defs.get_geoclim_path("climate_data")
clim_raster_filename = defs.geoclim_files["madaclim_current"]
default_clim_raster_path = climate_dir / clim_raster_filename

enviro_dir = defs.get_geoclim_path("environment_data") 
env_raster_filename = defs.geoclim_files["madaclim_enviro"]
default_env_raster_path = enviro_dir / env_raster_filename


class MadaclimPoint:
    
    #TODO DOCSTRING
    def __init__(self, specimen_id: str, latitude: float, longitude: float, source_crs: pyproj.crs.crs.CRS=pyproj.CRS.from_epsg(4326), **kwargs) -> None:
        self.specimen_id = specimen_id
        self.source_crs = self.validate_crs(source_crs)
        self.latitude = self.validate_lat(latitude, crs=self.source_crs)
        self.longitude = self.validate_lon(longitude, crs=self.source_crs)
        self.mada_geom_point = self._construct_point(
            latitude=self.latitude, 
            longitude=self.longitude,
            source_crs=self.source_crs
        )
        
        # Store any additional keyword arguments as instance attributes
        additional_args = [key for key in kwargs if key not in ["id", "epsg", "x", "y"]]
        for key in additional_args:
            setattr(self, key, kwargs[key])
    
    @property
    def specimen_id(self):
        return self._specimen_id
    @specimen_id.setter
    def specimen_id(self, value):
        self._specimen_id = value
    
    @property
    def latitude(self):
        return self._latitude
    
    @latitude.setter
    def latitude(self, value):
        value = self.validate_lat(value, crs=self.source_crs)
        self._latitude = value
        # Update mada_geom_point when latitude is updated
        self._update_mada_geom_point()
    
    @property
    def longitude(self):
        return self._longitude
    
    @longitude.setter
    def longitude(self, value):
        value = self.validate_lon(value, crs=self.source_crs)
        self._longitude = value
        # Update mada_geom_point when longitude is updated
        self._update_mada_geom_point()
    
    @property
    def source_crs(self):
        return self._source_crs
    
    @source_crs.setter
    def source_crs(self, value):
        value = self.validate_crs(value)
        self._source_crs = value
        
        # Update mada_geom_point when crs is updated
        self._update_mada_geom_point()
    
    def validate_crs(self, crs):
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
        # Validate float type
        try:
            latitude = float(latitude)
        except:
            raise TypeError(f"Could not convert {latitude} to float. Latitude must be a float.")
        
        # Validate bounds according to crs
        bounds = crs.area_of_use.bounds
        min_lat, max_lat = bounds[0], bounds[2]
        if not min_lat <= latitude <= max_lat:
            raise ValueError(f"{latitude=} is out of bounds for the crs=EPSG:{crs.to_epsg()}. Latitude must be between: {min_lat} and {max_lat}")
        
        return latitude
    
    def validate_lon(self, longitude, crs):
        # Validate float type
        try:
            longitude = float(longitude)
        except:
            raise TypeError(f"Could not convert {longitude} to float. longitude must be a float.")
        
        # Validate bounds according to crs
        bounds = crs.area_of_use.bounds
        min_lon, max_lon = bounds[1], bounds[3]
        if not min_lon <= longitude <= max_lon:
            raise ValueError(f"{longitude=} is out of bounds for the crs=EPSG:{crs.to_epsg()}. Longitude must be between {min_lon} and {max_lon}")
        
        return longitude
    
    def _construct_point(self, latitude: float, longitude: float, source_crs: pyproj.crs.crs.CRS)-> shapely.geometry.point.Point:

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
        if hasattr(self, '_latitude') and hasattr(self, '_longitude'):
            self.mada_geom_point = self._construct_point(
                latitude=self.latitude,
                longitude=self.longitude,
                source_crs=self.source_crs
            )

    
    def sample_from_rasters(
            self,
            layers_to_sample: Union[int, str, List[Union[int, str]]]="all", 
            layer_info: bool=True,
            return_nodata_layers: bool=False,
            clim_raster_path: Optional[pathlib.Path]=None, 
            env_raster_path: Optional[pathlib.Path]=None
        ) -> Union[dict, list]:
        """
        Samples geoclimatic data from raster files for specified layers at the location of the instances's lat/lon coordinates from the mada_geom_point attribute.

        Args:
            layers_to_sample (Union[int, str, List[Union[int, str]]], optional): The layer number(s) to sample from the raster files.
                Can be a single int, a single string in the format 'layer_<num>', or a list of ints or such strings. Defaults to 'all'.
            layer_info (bool, optional): Whether to use descriptive labels for the returned dictionary keys. Defaults to True.
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
        """
        
        # Create a MadaclimLayers instance to get layers labels and validate layers to sample
        madaclim_info = MadaclimLayers()
        all_layers_df = madaclim_info.all_layers
        
        # Validate layers to sample
        possible_layers_num_format = [f"layer_{num}" for num in all_layers_df["layer_number"].to_list()]

        if isinstance(layers_to_sample, list):
            # Check if all elements are in layer_<num> format
            layers_num_format = all([layer_label in possible_layers_num_format for layer_label in layers_to_sample])
            
            if layers_num_format:
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

        if clim_raster_layers_to_sample:
            start_time = time.perf_counter()
            total_clim_layers = len(clim_raster_layers_to_sample)
            
            with rasterio.open(clim_raster_path or default_clim_raster_path) as clim_raster:
                # Initialize reference and container to check for layers with nodata values
                nodata_clim = clim_raster.nodata
                
                # Status bar to display when sampling the raster
                print(f"Sampling {total_clim_layers} layer(s) from {clim_raster.name.split('/')[-1]}...")
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
                                    
            end_time = time.perf_counter()
            elapsed_time = end_time - start_time
            
            print(f"Finished {geoclim_types[0]} raster sampling operation in {elapsed_time:.2f} seconds\n")
            
        if env_raster_layers_to_sample:
            start_time = time.perf_counter()
            total_env_layers = len(env_raster_layers_to_sample)
            
            with rasterio.open(env_raster_path or default_env_raster_path) as env_raster:
                # Initialize reference and container to check for layers with nodata values
                nodata_env = env_raster.nodata
                
                # Status bar to display when sampling the raster
                print(f"Sampling {total_env_layers} layer(s) from {env_raster.name.split('/')[-1]}...")
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
                                    
            end_time = time.perf_counter()
            elapsed_time = end_time - start_time
            
            print(f"Finished {geoclim_types[1]} raster sampling operation in {elapsed_time:.2f} seconds\n")

        if len(nodata_layers) > 0:
                print(f"BEWARE! {len(nodata_layers)} layer(s) contain a nodata value at the specimen location")
        else:
            print("No sampled layers with nodata values.")
        
        if return_nodata_layers:
            return sampled_data, nodata_layers
        
        return sampled_data
    
    
    
    def __str__(self) -> str:
        madapoint_obj = (
            f"MadaclimPoint(\n\tspecimen_id = '{self.specimen_id}',\n\tsource_crs = EPSG:{self.source_crs.to_epsg()},\n\t"
            f"latitude = {self.latitude},\n\tlongitude = {self.longitude},\n\tmada_geom_point = {self.mada_geom_point}\n)"
        )
        return madapoint_obj

    def __repr__(self) -> str:
        madapoint_obj = (
            f"MadaclimPoint(\n\tspecimen_id = '{self.specimen_id}',\n\tsource_crs = EPSG:{self.source_crs.to_epsg()},\n\t"
            f"latitude = {self.latitude},\n\tlongitude = {self.longitude},\n\tmada_geom_point = {self.mada_geom_point}\n)"
        )
        return madapoint_obj