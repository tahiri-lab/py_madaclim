import pathlib
from typing import Optional, Union, List

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
clim_raster_path = climate_dir / clim_raster_filename

enviro_dir = defs.get_geoclim_path("environment_data") 
env_raster_filename = defs.geoclim_files["madaclim_enviro"]
env_raster_path = enviro_dir / env_raster_filename


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
        point = Point(latitude, longitude)
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