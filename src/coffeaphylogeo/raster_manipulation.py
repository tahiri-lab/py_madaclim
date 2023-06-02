import csv
import pathlib
import inspect
from pathlib import Path
from typing import Optional, Union, List, Optional, Tuple, Dict
import time
from tqdm import tqdm

from coffeaphylogeo._constants import Constants
from coffeaphylogeo.madaclim_info import MadaclimLayers

import rasterio
import pyproj
import shapely
from pyproj import Transformer
from shapely import Point
import geopandas as gpd

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

class MadaclimRaster:
    """Handles operations on Madaclim climate and environmental raster files.

    Attributes:
        clim_raster (pathlib.Path): Path to the climate raster file.
        env_raster (pathlib.Path): Path to the environmental raster file.
    """

    def __init__(self, clim_raster: pathlib.Path, env_raster: pathlib.Path) -> None:
        """
        Initializes a MadaclimRaster object with climate and environmental raster files.

        Args:
            clim_raster (pathlib.Path): Path to the climate raster file.
            env_raster (pathlib.Path): Path to the environmental raster file.
        
        Examples:
            >>> from coffeaphylogeo.raster_manipulation import MadaclimRaster
            >>> mada_rasters = MadaclimRaster(clim_raster="madaclim_current.tif", env_raster="madaclim_enviro.tif")
            >>> mada_rasters.clim_crs
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
    
    @property
    def clim_raster(self) -> pathlib.Path:
        """
        Retrieves or sets the climate raster file path.

        Args:
            value (pathlib.Path): The new climate raster file path.

        Raises:
            TypeError: If the input value is not a pathlib.Path or str.
            ValueError: If a pathlib.Path object could not be created from the input value.
            FileExistsError: If the file does not exist.
            IOError: If the raster file could not be opened.

        Returns:
            pathlib.Path: The climate raster file path.
        """
        return self._clim_raster
    
    @clim_raster.setter
    def clim_raster(self, value: pathlib.Path) -> None:
        # Validate type
        if not isinstance(value, (pathlib.Path, str)):
            raise TypeError("'clim_raster' must be a pathlib.Path object or str.")
        
        # Validate path and file
        try:
            value = Path(value)
        except:
            raise ValueError(f"Could not create a pathlib.Path object from {value}")
            
        # Check if raster file exists
        if not value.exists():
            raise FileExistsError(f"Could not find 'clim_raster' file: {value}")
              
        # Catch any IO errors
        try:
            with rasterio.open(value) as raster_file:
                self._clim_raster = value
        except rasterio.errors.RasterioIOError as e:
            raise IOError(f"Could not open 'clim_raster' file: {value}. Error: {e}")

    @property
    def env_raster(self) -> pathlib.Path:
        """
        Retrieves or sets the environmental raster file path.

        Args:
            value (pathlib.Path): The new environmental raster file path.

        Raises:
            TypeError: If the input value is not a pathlib.Path or str.
            ValueError: If a pathlib.Path object could not be created from the input value.
            FileExistsError: If the file does not exist.
            IOError: If the raster file could not be opened.

        Returns:
            pathlib.Path: The environmental raster file path.
        """
        return self._env_raster
    
    @env_raster.setter
    def env_raster(self, value: pathlib.Path) -> None:
        # Validate type
        if not isinstance(value, (pathlib.Path, str)):
            raise TypeError("'env_raster' must be a pathlib.Path object or str.")
        
        # Validate path and file
        try:
            value = Path(value)
        except:
            raise ValueError(f"Could not create a pathlib.Path object from {value}")
            
        # Check if raster file exists
        if not value.exists():
            raise FileExistsError(f"Could not find 'env_raster' file: {value}")
              
        # Catch any IO errors
        try:
            with rasterio.open(value) as raster_file:
                self._env_raster = value
        except rasterio.errors.RasterioIOError as e:
            raise IOError(f"Could not open 'env_raster' file: {value}. Error: {e}")
        
    @property
    def clim_crs(self) -> pyproj.crs.crs.CRS:
        """
        Retrieves the Coordinate Reference System (CRS) from the Madaclim climate raster.

        This property opens the raster file and retrieves the CRS in EPSG format. The EPSG code is used to 
        create and return a pyproj CRS object.

        Returns:
            pyproj.crs.
        """
        # Get epsg from clim_raster
        with rasterio.open(self.clim_raster) as clim_raster:
            clim_epsg = clim_raster.crs.to_epsg()  # Get the EPSG code of the CRS
            clim_crs = pyproj.CRS.from_epsg(clim_epsg)  # Create a pyproj CRS object
        return clim_crs
    @property
    def env_crs(self) -> pyproj.crs.crs.CRS:
        """
        Retrieves the Coordinate Reference System (CRS) from the Madaclim environmental raster.

        This property opens the raster file and retrieves the CRS in EPSG format. The EPSG code is used to 
        create and return a pyproj CRS object.

        Returns:
            pyproj.crs.
        """
        # Get epsg from clim_raster
        with rasterio.open(self.env_raster) as env_raster:
            env_epsg = env_raster.crs.to_epsg()  # Get the EPSG code of the CRS
            env_crs = pyproj.CRS.from_epsg(env_epsg)  # Create a pyproj CRS object
        return env_crs
        
    def __str__(self) -> str:
        info = "\n ".join({f"{attr} = {val}" for attr, val in vars(self).items()})
        return info
    
    def __repr__(self) -> str:
        return self.__str__()
    
    def visualize_layer(self, layer: Union[str, int], figsize: Optional[Tuple[int, int]]=None, **kwargs) -> None:
        """
        Method to visualize a specific layer from the Madagascan climate/environmental raster datasets. The layer 
        is displayed as a raster map and its distribution is plotted in a histogram. 

        It accepts layer labels in the following formats: `layer_<num>` (e.g. "layer_1") and `<descriptive_layer_label>` 
        (e.g. "annual_mean_temperature"). Alternatively, the layer number can be supplied directly as an integer.
        
        Depending on whether the layer is categorical or continuous, the visualization will be different. For categorical 
        layers, it will display a map using different colors for each category and a legend mapping categories to colors. 
        For continuous layers, it will display a color gradient map with a color bar.
        
        This function does not return anything; it directly generates and displays the plots.

        Args:
            layer (Union[str, int]): Layer to visualize. Accepts layer numbers as integers, or layer labels in 
                                    descriptive or `layer_<num>` format.
            figsize (Optional[Tuple[int, int]]): Tuple representing the size of the figure to be displayed. If not provided,
                                                a default size is used based on data type.
            **kwargs: Additional arguments to customize the imshow, colorbar and histplot. 
                    Use "imshow_<arg>", "cax_<arg>", and "histplot_<arg>" formats to customize corresponding 
                    matplotlib/sns arguments.

        Raises:
            TypeError: If 'layer' is not a str or an int, or if 'figsize' is not a tuple of ints.
            ValueError: If 'layer' is not found within the range of layers, or if 'figsize' is not a tuple of 2 elements.
            ValueError: If the keys from `get_categorical_combinations` do not match raster values.

        #TODO FINISH EXAMPLES
        Example:
            >>> # Extract environmental layers labels
            >>> from coffeaphylogeo.madaclim_info import MadaclimLayers
            >>> mada_info = MadaclimLayers(clim_raster="madaclim_current.tif", env_raster="madaclim_enviro.tif")
            >>> env_layers_labels = mada_info.get_layers_labels("env", as_descriptive_labels=True)
            >>> >>> env_layers_labels[0]    # Using alttitude as our example
            'env_71_altitude (Altitude in meters)'
            
            >>> # Default visualization of the raster map
            >>> from coffeaphylogeo.raster_manipulation import MadaclimRaster
            >>> mada_rasters = MadaclimRaster(clim_raster=mada_info.clim_raster, env_raster=mada_info.env_raster)    # Using common attr btw the instances
            >>> mada_rasters.v
            mada_raster.visualize_layer('layer_2', figsize=(12, 6), imshow_cmap='hot', 
                            cax_size='7%', histplot_bins=30, histplot_color='blue')
        """
    
        # Fetch metadata and layer nums with a MadaclimLayers instance
        madaclim_info = MadaclimLayers(clim_raster=self.clim_raster, env_raster=self.env_raster)
        all_layers_df = madaclim_info.all_layers
        
        # Validate layers to sample
        possible_layers_num_format = madaclim_info.get_layers_labels()
        possible_layers_desc_format = madaclim_info.get_layers_labels(as_descriptive_labels=True)


        if not isinstance(layer, (str, int)):
            raise TypeError("'layer' must be a str or an int.")
        
        if layer in possible_layers_num_format or layer in possible_layers_desc_format:    # Check layer_<num> or descriptive format
            layer_num = int(layer.split("_")[1])
        else:
            try:
                layer_num = int(layer)    # As single item int list
            except (ValueError, TypeError):
                raise TypeError("layer must be either a single int value or a string that can be converted to an int")
  
        # Validate layer number range for layer_numbers as in
        min_layer = min(all_layers_df["layer_number"])
        max_layer = max(all_layers_df["layer_number"])

        if not min_layer <= layer_num <= max_layer:
            raise ValueError(f"layer_number must fall between {min_layer} and {max_layer}. {layer_num=} is not valid.")

        # Check if figsize is a tuple of size 2 with proper types
        if figsize:
            if not isinstance(figsize, tuple):
                raise TypeError("'figsize' must be a tuple.")
            if not len(figsize) == 2:
                raise ValueError("'figsize' must be a tuple of 2 elements.")
            for ele in figsize:
                if not isinstance(ele, int):
                    raise TypeError(f"{ele} is not an int. 'figsize' is a tuple of integers.")    
        
        # Get customizable matplotlib kwargs for both axes
        imshow_args = {k[7:]: v for k, v in kwargs.items() if k.startswith("imshow_")}
        cax_args = {k[4:]: v for k, v in kwargs.items() if k.startswith("cax_")}
        histplot_args = {k[9:]: v for k, v in kwargs.items() if k.startswith("histplot_")}

        # Set default for imshow
        imshow_cmap = imshow_args.pop("cmap", "inferno")

        # Set defaults cbar
        cax_position = cax_args.pop("position", "right")
        cax_size = cax_args.pop("size", "5%")
        cax_pad = cax_args.pop("pad", 0.10)

        # Set defaults for histplot
        histplot_color = histplot_args.pop("color", "grey")
        histplot_bins = histplot_args.pop("bins", "auto")
        histplot_kde = histplot_args.pop("kde", True)
        histplot_stat = histplot_args.pop("stat", "percent")
        histplot_line_kws = histplot_args.pop("line_kws", {"linestyle" : "--"})
        
        # Plotting the raster map with distribution data for the selected layer
        band_num = madaclim_info.get_bandnums_from_layers(layer_num)[0]
        layer_info = madaclim_info.fetch_specific_layers(layer_num, "layer_description", "is_categorical", "units")
        layer_key = list(layer_info.keys())[0]
        layer_info = layer_info[layer_key]    # Eliminate redundancy

        layer_units = "Categorical values" if layer_info["is_categorical"] else layer_info["units"]    # Simplify xlabel for categorical units
        layer_description = layer_info["layer_description"]

        # Get geoclim type from layer_num for raster IO selection            
        geoclim_types = ["clim", "env"]
        geoclim_type = [geotype for geotype in geoclim_types if layer_num in madaclim_info.select_geoclim_type_layers(geotype)["layer_number"].values][0]
        chosen_raster = self.clim_raster if geoclim_type == "clim" else self.env_raster

        with rasterio.open(chosen_raster) as raster:
            band_data = raster.read(band_num, masked=True)

            # Categorical data multi plots
            if layer_info["is_categorical"]:
            
                figsize = (20, 10) if figsize is None else figsize
                fig, axes = plt.subplots(1, 2, figsize=figsize)

                # Custom categorical palette for over 20 categories
                categ_combinations = dict(sorted(madaclim_info.get_categorical_combinations(layer_num).items(), key=lambda item: item[0]))   # Get categ dict {<numerical_val>:<category>}
                if len(categ_combinations) > 20:    
                    colors = [
                        "#ff0000", "#3cb44b", "#ffe432", "#4363d8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", 
                        "#bcf60c", "#fabebe", "#008080", "#bbbbbb", "#9a6324", "#fffac8", "#800000", "#aaffc3",
                        "#808000", "#ffd8b1", "#000075", "#e6beff", "#ffffff", "#000000", "#bda400", "#bd5447",
                        "#354555" 
                    ]
                else:
                    colors = sns.color_palette(palette="tab20", n_colors=len(categ_combinations))

                # Exclude masked values
                unmasked_values = [v for v in categ_combinations.keys() if v is not np.ma.masked]
                unmasked_colors = colors[:len(unmasked_values)]

                # Make sure categories are sorted and mapped to colors accordingly
                sorted_categ_combinations = dict(sorted(zip(unmasked_values, unmasked_values)))
                categ_colors = {category: color for category, color in zip(sorted_categ_combinations.values(), unmasked_colors)}
                cmap = ListedColormap(categ_colors.values())

                # Create boundaries based on your categories
                bounds = list(sorted_categ_combinations.keys())
                norm = BoundaryNorm(bounds, cmap.N)

                # Plot the raster map
                axes[0].imshow(band_data.squeeze(), cmap=cmap, norm=norm, **imshow_args)    # Draw raster map from masked array

                # Create the legend with categorical labels
                categ_labels = [categ_combinations[cat] for cat in sorted_categ_combinations.keys()]
                legend_elements = [Patch(facecolor=color, edgecolor=color, 
                                        label=category) for color, category in zip(categ_colors.values(), categ_labels)]
                axes[0].legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left")
                                
                axes[0].set_title(f"Madagascar {geoclim_type.capitalize()} Raster Map (band={band_num})", fontsize=10)    # Raster map interface cleanup
                axes[0].set_yticks([])
                axes[0].set_xticks([])
                axes[0].axis("off")
                
                # Compute categories counts, normalize and map colors
                categ, counts = np.unique(band_data.compressed(), return_counts=True)
                counts_percent = (counts / counts.sum()) * 100
                categ = categ.astype(int)
                sort_idx = np.argsort(categ)
                categ = categ[sort_idx]
                counts_percent = counts_percent[sort_idx]

                mapped_categ_combinations = [sorted_categ_combinations[c] if c in sorted_categ_combinations else None for c in categ]
                if None in mapped_categ_combinations:
                    raise ValueError("Keys from `get_categorical_combinations` does not match raster values")

                # Normalize counts to freq percent
                df_cat = pd.DataFrame({"category":categ, "mapped_category": mapped_categ_combinations, "counts_percent": counts_percent})

                # Categorical bar plot
                axes[1].bar(df_cat["category"], df_cat["counts_percent"], color=categ_colors.values())
                axes[1].set_title("Distribution of raster values at 1km resolution", fontsize=10)
                xticks = list(set(df_cat["category"].tolist() + [1]))
                axes[1].set_xticks(sorted(xticks))
                axes[1].set_xlabel(layer_units)
                axes[1].set_ylabel("Percent (%)")

            # Continous data for multi plots
            else:                  
                figsize = (12, 6) if figsize is None else figsize
                fig, axes = plt.subplots(1, 2, figsize=figsize)

                # Raster map with cbar
                imshow_cmap = imshow_args.pop("cmap", "inferno")
                im = axes[0].imshow(band_data.squeeze(), cmap=imshow_cmap, **imshow_args)    # Draw raster map from masked array
                
                # Colorbar customization
                divider = make_axes_locatable(axes[0])
                cax = divider.append_axes(position=cax_position, size=cax_size, pad=cax_pad, **cax_args)
                plt.colorbar(im, cax=cax)

                axes[0].set_title(f"Madagascar {geoclim_type.capitalize()} Raster Map (band={band_num})", fontsize=10)
                axes[0].set_yticks([])
                axes[0].set_xticks([])
                axes[0].axis("off")
                
                fig.text(    # Add units to raster map
                    x=(cax.get_position().x0 + cax.get_position().x1) / 2,
                    y=cax.get_position().y0 - 0.05, 
                    s=layer_units, 
                    ha="center", 
                    va="center"
                )
                
                # Plot histogram of raster vals
                sns.histplot(
                    data=band_data.compressed(),
                    ax=axes[1], 
                    color=histplot_color, 
                    bins=histplot_bins, 
                    kde=histplot_kde, 
                    stat=histplot_stat,
                    line_kws=histplot_line_kws,
                    **histplot_args
                )
                axes[1].lines[0].set_color("black")
                axes[1].set_title("Distribution of raster values at 1km resolution", fontsize=10)
                axes[1].set_xlabel(layer_units)

            fig.suptitle(
                f"Layer {layer_num}: {layer_description}",
                fontsize=16,
                weight="bold",
                ha="center"
            )

            fig.tight_layout()


class MadaclimPoint:
    
    """
    A class representing a specimen as a geographic point with a specific coordinate reference system (CRS)
    and additional attributes. The class provides methods for validating the point's coordinates
    and CRS, as well as sampling values from climate and environmental rasters of the Madaclim database.
    
    Attributes:
        specimen_id (str): An identifier for the point.
        latitude (float): The latitude of the point.
        longitude (float): The longitude of the point.
        #TODO CLIM/ENV_CRS + CLIM/ENV_RASTER ATTRS
        source_crs (pyproj.crs.crs.CRS): The coordinate reference system of the point.
        mada_geom_point (shapely.geometry.point.Point): A Shapely Point object representing the point projected in the Madaclim rasters' CRS.
        ___base_attr (dict): A dictionary containing the base attributes names as keys and their values as values.
        
    """

    def __init__(
            self, 
            specimen_id: str, 
            latitude: float, 
            longitude: float, 
            source_crs: pyproj.crs.crs.CRS=pyproj.CRS.from_epsg(4326), 
            **kwargs
        ) -> None:
        """
        Initialize a MadaclimPoint object with the given `specimen_id`, `latitude`, `longitude`, `clim_raster, `env_raster` and `source_crs`.
        The coordinates provided should respect the nature of the given source CRS native units' (i.e. degrees WGS84 or meters for EPSG:3857).
        The path to both climate and environmental rasters have to be provided to properly construc the `mada_geom_point` attribute and use other methods.
        Optionally, provide additional keyword arguments to store as instance attributes.
        
        Args:
            specimen_id (str): An identifier for the point.
            latitude (float): The latitude of the point.
            longitude (float): The longitude of the point.
            source_crs (pyproj.crs.crs.CRS, optional): The coordinate reference system of the point. Defaults to WGS84 (EPSG:4326).
            **kwargs: Additional keyword arguments to store as instance attributes.
        Examples:
            >>> from coffeaphylogeo.raster_manipulation import MadaclimPoint
            >>> specimen_1 = MadaclimPoint(specimen_id="spe1_aren", latitude=-18.9333, longitude=48.2)    # Default CRS of EPSG:4326
            
            >>> # Creates a shapely point object according to the Madaclim CRS' projection when instantiating a new instance. 
            >>> specimen_1.mada_geom_point
            <POINT (837072.915 7903496.321)>

            >>> # Uses default EPSG:4326 CRS when not specified
            >>> specimen_1.source_crs
            <Geographic 2D CRS: EPSG:4326>
            Name: WGS 84
            Axis Info [ellipsoidal]:
            - Lat[north]: Geodetic latitude (degree)
            - Lon[east]: Geodetic longitude (degree)
            Area of Use:
            - name: World.
            - bounds: (-180.0, -90.0, 180.0, 90.0)
            Datum: World Geodetic System 1984 ensemble
            - Ellipsoid: WGS 84
            - Prime Meridian: Greenwich

            >>> # Also accepts any other kwargs and saves them as attributes with specific typing
            >>> specimen_1 = MadaclimPoint(specimen_id="spe1_aren", latitude=-18.9333, longitude=48.2, genus="Coffea", species="arenesiana", has_sequencing=True)
            >>> specimen_1.species
            'arenesiana'
            >>> specimen_1.has_sequencing
            True
            >>> specimen_1
            MadaclimPoint(
                specimen_id = spe1_aren,
                source_crs = 4326,
                latitude = -18.9333,
                longitude = 48.2,
                mada_geom_point = POINT (837072.9150244407 7903496.320897499),
                sampled_data = None,
                nodata_layers = None,
                genus = Coffea,
                species = arenesiana,
                has_sequencing = True
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

        self._sampled_data = None
        self._nodata_layers = None

        # Store any additional keyword arguments as instance attributes
        base_args = self.get_args_names()[0] + [self.get_args_names()[1]]
        additional_args = [key for key in kwargs if key not in base_args]
        for key in additional_args:
            setattr(self, key, kwargs[key])

        # self.gdf = {
        #     "test": "geopandas_gdf"
        # }
    
    @property
    def specimen_id(self) -> str:
        """
        Get or sets the specimen_id attribute.

        Args:
            value (str): The new identifier for the MadaclimPoint.

        Returns:
            str: The identifier for the MadaclimPoint.
        """
        return self._specimen_id
    @specimen_id.setter
    def specimen_id(self, value: str):
        self._specimen_id = value
    
    @property
    def latitude(self) -> float:
        """
        Get or sets the latitude attribute.

        Args:
            value (float): The latitude value for the point.

        Returns:
            float: The latitude of the point.
        """
        return self._latitude
    
    @latitude.setter
    def latitude(self, value: float):
        value = self.validate_lat(value, crs=self.source_crs)
        self._latitude = value
        # Update mada_geom_point when latitude is updated
        self._update_mada_geom_point()
    
    @property
    def longitude(self) -> float:
        """
        Get or sets the longitude attribute.

        Args:
            value (float): The longitude value for the point.

        Returns:
            float: The longitude of the point.
        """
        return self._longitude
    
    @longitude.setter
    def longitude(self, value: float):
        value = self.validate_lon(value, crs=self.source_crs)
        self._longitude = value
        # Update mada_geom_point when longitude is updated
        self._update_mada_geom_point()
    
    @property
    def source_crs(self) -> pyproj.crs.CRS:
        """
        Get or sets the source_crs attribute.
       
         Args:
            value (pyproj.crs.CRS): The coordinate reference system for the point.
    
        Returns:
            pyproj.crs.CRS: The coordinate reference system of the point.
        """
        return self._source_crs
    
    @source_crs.setter
    def source_crs(self, value: pyproj.crs.CRS):
        value = self.validate_crs(value)
        self._source_crs = value
        
        # Update mada_geom_point when crs is updated
        self._update_mada_geom_point()

    @property
    def base_attr(self) -> dict:
        """Get the base attributes when constructing the instance.

        Returns:
            dict: A dictionary containing the base attributes names as keys and their values as values.
        """
        return self.__base_attr
    
    # @property
    # def gdf(self) -> gpd.GeoDataFrame:
    #     """Get the geodataframe using mada_geom_point as geometry.

    #     Returns:
    #         gpd.GeoDataFrame: A Geopandas GeoDataFrame generated from instance attributes and Point geometry.
    #     """
    #     return self.__gdf

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
    
    def validate_crs(self, crs) -> pyproj.crs.crs.CRS:
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
    
    def validate_lat(self, latitude, crs) -> float:
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
    
    def validate_lon(self, longitude, crs) -> float:
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
            clim_raster: pathlib.Path, 
            env_raster: pathlib.Path,
            layers_to_sample: Union[int, str, List[Union[int, str]]]="all", 
            layer_info: bool=False,
            return_nodata_layers: bool=False,
        ) -> Union[dict, list]:
        """
        Samples geoclimatic data from raster files for specified layers at the location of the instances's lat/lon coordinates from the mada_geom_point attribute.

        Args:
            clim_raster_path (pathlib.Path): Path to the climate raster file.
            env_raster_path (pathlib.Path): Path to the environment raster file.
            layers_to_sample (Union[int, str, List[Union[int, str]]], optional): The layer number(s) to sample from the raster files.
                Can be a single int, a single string in the format 'layer_<num>', or a list of ints or such strings. Defaults to 'all'.
            layer_info (bool, optional): Whether to use descriptive labels for the returned dictionary keys. Defaults to False.
            return_nodata_layers (bool, optional): Whether to return a list of layers with nodata values at the specimen location.
                Defaults to False.

        Raises:
            TypeError: If the layers_to_sample is not valid, or if the mada_geom_point attribute is not a Point object.
            ValueError: If the layer_number is out of range or if the mada_geom_point object is empty.

        Returns:
            Union[dict, list]: A dictionary containing the sampled data, with keys being layer names or numbers depending
                on the layer_info parameter. If return_nodata_layers is True, also returns a list of layers with nodata values
                at the specimen location.
        Exmaples:
            >>> # Fetching bioclim layers from the MadaclimLayers class
            >>> from coffeaphylogeo.madaclim_info import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> bioclim_labels = [label for label in madaclim_info.get_layers_labels(as_descriptive_labels=True) if "bio" in label]

            >>> # Sampling the bioclim layers
            >>> specimen_1 = MadaclimPoint(specimen_id="spe1_aren", latitude=-18.9333, longitude=48.2, genus="Coffea", species="arenesiana", has_sequencing=True)
            >>> spe1_bioclim = specimen_1.sample_from_rasters(
            ...     clim_raster="madaclim_current.tif",
            ...     env_raster="madaclim_enviro.tif",
            ...     layers_to_sample=bioclim_labels
            ... )
            >>> spe1_bioclim["layer_37"]
            196

            >>> # layer_info key as more descriptive and informative
            >>> spe1_bioclim = specimen_1.sample_from_rasters(
            ...     clim_raster="madaclim_current.tif",
            ...     env_raster="madaclim_enviro.tif",
            ...     layers_to_sample=bioclim_labels,
            ...     layer_info=True
            ... )

            >>> bio1_label = bioclim_labels[0]    # Example for first bioclim variable
            'clim_37_bio1 (Annual mean temperature)'
            >>> spe1_bioclim[bio1_label]
            196
            >>> {k:v for k, v in list(spe1_bioclim.items())[:3]}    # Print first 3 extracted items
            {'clim_37_bio1 (Annual mean temperature)': 196, 'clim_38_bio2 (Mean diurnal range (mean of monthly (max temp - min temp)))': 112, 'clim_39_bio3 (Isothermality (BIO2/BIO7) (x 100))': 64}

            
            >>> # Sample any rasters using layer numbers only
            >>> spe1_l68_l71 = specimen_1.sample_from_rasters(
            ...     clim_raster="madaclim_current.tif",
            ...     env_raster="madaclim_enviro.tif",
            ...     layers_to_sample=[68, 71],
            ...     layer_info=True
            ... )

            >>> spe1_l68_l71
            {'clim_68_pet (Annual potential evapotranspiration from the Thornthwaite equation (mm))': 891, 'env_71_altitude (Altitude in meters)': 899}

            >>> # Sample all layers with less descriptive layer names
            >>> specimen_2 = MadaclimPoint(specimen_id="spe2_humb", latitude=-12.716667, longitude=45.066667, source_crs=4326, genus="Coffea", species="humblotiana", has_sequencing=True)
            >>> spe2_all_layers = specimen_2.sample_from_rasters("madaclim_current.tif", "madaclim_enviro.tif")

            ######################################## Extracting data for: spe2_humb ########################################

            Sampling 70 layer(s) from madaclim_current.tif (geoclim_type=clim)...
            Extracting layer 70: ndm (Number of dry months in the year):  100%|██████████████████████████████████████████████████████████████████████████████████████| layer 70/70 [Time remaining: 00:00]

            Sampling 9 layer(s) from madaclim_enviro.tif (geoclim_type=env)...
            Extracting layer 79: forestcover (Percentage of forest cover in 1 km by 1 km grid cells):  100%|███████████████████████████████████████████████████████████| layer 9/9 [Time remaining: 00:00]
            BEWARE! 5 layer(s) contain a nodata value at the specimen location

            Finished raster sampling operation in 0.09 seconds.

            >>> spe2_all_layers["layer_68"]
            1213
            
            >>> # Note the BEWARE! message indicating that we have some NaN in the data extracted for that MadaclimPoint
            >>> # We can easily access the nodata layers (still sampled with the method regardless)
            >>> spe2_all_layers, spe2_nodata_layers = specimen_2.sample_from_rasters(
            ...     clim_raster="madaclim_current.tif",
            ...     env_raster="madaclim_enviro.tif",
            ...     layer_info=True, 
            ...     return_nodata_layers=True
            ... )
            >>> len(spe2_nodata_layers)
            5
            >>> spe2_nodata_layers[0]    # Example of a categorical feature description with raster-value/description associations
            'env_75_geology (1=Alluvial_&_Lake_deposits, 2=Unconsolidated_Sands, 4=Mangrove_Swamp, 5=Tertiary_Limestones_+_Marls_&_Chalks, 6=Sandstones, 7=Mesozoic_Limestones_+_Marls_(inc._"Tsingy"), 9=Lavas_(including_Basalts_&_Gabbros), 10=Basement_Rocks_(Ign_&_Met), 11=Ultrabasics, 12=Quartzites, 13=Marble_(Cipolin))'

            >>> # Calling the sample_from_rasters method also updates the 'sampled_data' and 'nodata_layers' attributes
            >>> specimen_2.sample_from_rasters(
            ...     clim_raster="madaclim_current.tif",
            ...     env_raster="madaclim_enviro.tif",
            ...     layers_to_sample=[37, 75]
            ... )
            >>> specimen_2
            MadaclimPoint(
                    specimen_id = spe2_humb,
                    source_crs = 4326,
                    latitude = -12.716667,
                    longitude = 45.066667,
                    mada_geom_point = POINT (507237.57495924993 8594195.741515966),
                    sampled_data = {'layer_37': 238, 'layer_75': -32768},
                    nodata_layers = ['layer_75'],
                    genus = Coffea,
                    species = humblotiana,
                    has_sequencing = True
            )
        """
        # Create a MadaclimRaster to validate both rasters
        mada_raster = MadaclimRaster(clim_raster=clim_raster, env_raster=env_raster)
        if mada_raster.clim_crs != Constants.MADACLIM_CRS:
            raise ValueError(f"The provided clim_raster's CRS does not corresponds to Madaclim db's expected crs: {Constants.MADACLIM_CRS}")

        if mada_raster.env_crs != Constants.MADACLIM_CRS:
            raise ValueError(f"The provided env_raster's CRS does not corresponds to Madaclim db's expected crs: {Constants.MADACLIM_CRS}")
        
        # Create a MadaclimLayers instance to get layers labels and validate layers to sample
        madaclim_info = MadaclimLayers(clim_raster=mada_raster.clim_raster, env_raster=mada_raster.env_raster)
        all_layers_df = madaclim_info.all_layers
        
        # Validate layers to sample
        possible_layers_num_format = madaclim_info.get_layers_labels()
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
            
            elif layers_to_sample in possible_layers_num_format or layers_to_sample in possible_layers_desc_format:    # Check layer_<num> or desc str format
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
            
        # Save layers and band numbers to sample for each raster
        clim_raster_sample_info, env_raster_sample_info = {}, {}

        geoclim_types = ["clim", "env"]
        geoclim_layer_ranges = {geoclim_type: madaclim_info.select_geoclim_type_layers(geoclim_type)["layer_number"].to_list() for geoclim_type in geoclim_types}

        clim_raster_sample_info["layers"] = [layer_num for layer_num in layers_numbers if layer_num in geoclim_layer_ranges["clim"]]
        clim_raster_sample_info["bands"] = madaclim_info.get_bandnums_from_layers(clim_raster_sample_info["layers"])

        env_raster_sample_info["layers"] = [layer_num for layer_num in layers_numbers if layer_num in geoclim_layer_ranges["env"]]
        env_raster_sample_info["bands"] = madaclim_info.get_bandnums_from_layers(env_raster_sample_info["layers"])

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
        
        if clim_raster_sample_info["layers"]:
            total_clim_layers = len(clim_raster_sample_info["layers"])
            
            with rasterio.open(mada_raster.clim_raster) as clim_raster:
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
                    for layer_num, band_num in zip(clim_raster_sample_info["layers"], clim_raster_sample_info["bands"]):
                        # Get layer info for pbar display and label key for sampled_data
                        layer_name = madaclim_info.fetch_specific_layers(layers_labels=layer_num, as_descriptive_labels=True, return_list=True)[0]
                        layer_description_display = layer_name.split('_')[-1]
                        pbar.set_description(f"Extracting layer {layer_num}: {layer_description_display}")
                        pbar.update()
                        
                        # Sample using the self.mada_geom_point attributes coordinates for the current layer
                        data = list(clim_raster.sample([(self.mada_geom_point.x, self.mada_geom_point.y)], indexes=band_num))[0]
                        
                        # Save extracted data with specified layer info/name
                        if not layer_info:
                            layer_name = f"layer_{layer_num}"
                        sampled_data[layer_name] = data[0]

                        if data[0] == nodata_clim:    # Save layers where nodata at specimen location
                            nodata_layers.append(layer_name)
        
        if env_raster_sample_info["layers"]:
            total_env_layers = len(env_raster_sample_info["layers"])
            
            with rasterio.open(mada_raster.env_raster) as env_raster:
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
                    for layer_num, band_num in zip(env_raster_sample_info["layers"], env_raster_sample_info["bands"]):
                        # Get layer info for pbar display and label key for sampled_data
                        layer_name = madaclim_info.fetch_specific_layers(layers_labels=layer_num, as_descriptive_labels=True, return_list=True)[0]
                        layer_description_display = layer_name.split('_')[-1]
                        pbar.set_description(f"Extracting layer {layer_num}: {layer_description_display}")
                        pbar.update()
                        
                        # Sample using the self.mada_geom_point attributes coordinates for the current layer
                        data = list(env_raster.sample([(self.mada_geom_point.x, self.mada_geom_point.y)], indexes=band_num))[0]
                        
                        # Save extracted data with specified layer info/name
                        if not layer_info:
                            layer_name = f"layer_{layer_num}"
                        sampled_data[layer_name] = data[0]

                        if data[0] == nodata_env:    # Save layers where nodata at specimen location
                            nodata_layers.append(layer_name)
                                    

        if len(nodata_layers) > 0:    # No raising exception, just warning print
            print(f"BEWARE! {len(nodata_layers)} layer(s) contain a nodata value at the specimen location")

        end_time = time.perf_counter()
        elapsed_time = end_time - start_time    # Total raster sampling time
        print(f"\nFinished raster sampling operation in {elapsed_time:.2f} seconds.\n")
        
        # Update instance attributes
        self._sampled_data = sampled_data
        self._nodata_layers = nodata_layers if len(nodata_layers) > 0 else None

        if return_nodata_layers:
            return sampled_data, nodata_layers
        
        return sampled_data
    
    def visualize_on_layer(layer: Union[int, str], **kwargs) -> None:
        pass
        # TODO: GET BAND NUMS FROM LAYER
        # TODO: IMPLEMENT MATPLOTLIB / GEOPANDAS PLOT
        # TODO: IMPLEMENT MATPLOTLIB CUSTOM
    
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
        madaclim_crs = Constants.MADACLIM_CRS
        
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
            >>> collection
            MadaclimCollection = [
                MadaclimPoint(specimen_id=spe1, mada_geom_point=POINT (644890.8921103649 7392153.658976035)),
                MadaclimPoint(specimen_id=spe2, mada_geom_point=POINT (536050.6239664567 7465523.013290589))
            ]
            >>> # You can also initiliaze an empty collection
            >>> collection = MadaclimCollection()
            >>> collection
            No MadaclimPoint inside the collection.

            >>> # Add a single MadaclimPoint
            >>> collection.add_points(specimen_1)
            >>> collection
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
            >>> collection
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
                mada_geom_point = POINT (837072.9150244407 7903496.320897499),
                sampled_data = None,
                nodata_layers = None
            )

        """
        self.__all_points = []
        if madaclim_points:
            self.add_points(madaclim_points)
        self.__sampled_raster_data = None
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
                specimen_id = sample_A,
                source_crs = 4326,
                latitude = -18.9333,
                longitude = 48.2,
                mada_geom_point = POINT (837072.9150244407 7903496.320897499),
                sampled_data = None,
                nodata_layers = None
            )
        """
        return self.__all_points
    
    @property
    def sampled_raster_data(self) -> Union[None, Dict[str, Dict[str, float]]]:
        """Get the sampled_raster_data attribute.

        This attribute is a nested dictionary. The outer dictionary uses the MadaclimPoint.specimen_id as keys. 
        The corresponding value for each key is another dictionary, which uses layer_names as keys and sampled values from rasters as values.

        Returns:
            Union[None, Dict[str, Dict[str, float]]]: A dictionary with MadaclimPoint.specimen_id as keys and a dictionary of layer_names (str) and sampled values (float) as values.
                None if Collection has not been sampled yet.
        """
        return self.__sampled_raster_data

    @property
    def nodata_layers(self) -> Union[None, Dict[str, Union[str, List[str]]]]:
        """Get the nodata_layers attribute.

        This attribute is a dictionary that contains the MadaclimPoint.specimen_id as keys and the values as the 'nodata_layers' as str or list of str.
        
        Returns:
            Dict[str, Union[str, List[str]]]: A dictionary with MadaclimPoint.specimen_id as keys and values of str or list of str of the layers_name with nodata values.
                None if Collection has not been sampled yet or all layers sampled contained valid data.
        """
        return self.__nodata_layers
    
    def __str__(self) -> str:
        if len(self.__all_points) == 0:
            return "No MadaclimPoint inside the collection."
        else:
            # Raster sampling state of each object of the collection.
            sampled = False
            if self.sampled_raster_data is not None:
                sampled = True

            all_points_short = [
                f"MadaclimPoint(specimen_id={point.specimen_id}, mada_geom_point={point.mada_geom_point}, {sampled=})" for point in self.all_points
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
            >>> # csv headers must contain the 3 required positional args for MadaclimPoint (specimen_id, lat, lon)
            >>> # When no source_crs header is found, defaults to EPSG:4326

            >>> # CSV example #1: "some_samples.csv"
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
            >>> collection
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_A, mada_geom_point=POINT (837072.9150244407 7903496.320897499), sampled=False),
                MadaclimPoint(specimen_id=sample_B, mada_geom_point=POINT (695186.2170220022 8197477.647690434), sampled=False),
                MadaclimPoint(specimen_id=sample_C, mada_geom_point=POINT (761613.8281386737 7651088.106452912), sampled=False),
                MadaclimPoint(specimen_id=sample_D, mada_geom_point=POINT (955230.600222457 8005985.896187438), sampled=False),
                MadaclimPoint(specimen_id=sample_E, mada_geom_point=POINT (757247.2175273325 7618631.528869408), sampled=False)
            ]

            >>> # Can accept other non-required data for MadaclimPoint instantiation
            
            >>> # CSV example #2: "other_samples.csv"
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
                sampled_data = None,
                nodata_layers = None,
                has_sequencing = True,
                specie = bojeri
            )
            MadaclimPoint(
                specimen_id = sample_G,
                source_crs = 4326,
                latitude = -18.295741,
                longitude = 45.826763,
                mada_geom_point = POINT (587378.6907481698 7976896.406900212),
                sampled_data = None,
                nodata_layers = None,
                has_sequencing = False,
                specie = periwinkle
            )
            MadaclimPoint(
                specimen_id = sample_H,
                source_crs = 4326,
                latitude = -21.223,
                longitude = 44.5204,
                mada_geom_point = POINT (450229.7195355138 7653096.609718417),
                sampled_data = None,
                nodata_layers = None,
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
                mada_geom_point = POINT (695186.2170220022 8197477.647690434),
                sampled_data = None,
                nodata_layers = None
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
            >>> collection
            No MadaclimPoint inside the collection.

            >>> # Add a single point
            >>> collection.add_points(specimen_1)
            >>> collection
            MadaclimCollection = [
                MadaclimPoint(specimen_id=spe1, mada_geom_point=POINT (644890.8921103649 7392153.658976035), sampled=False)
            ]

            >>> # Add multiple points
            >>> specimen_2 = MadaclimPoint(specimen_id="spe2", latitude=-20.138470, longitude=46.054688, family="Rubiaceae", has_sequencing=True, num_samples=1)            
            >>> other_collection.add_points([specimen_1, specimen_2])
            >>> print(other_collection)
            MadaclimCollection = [
                MadaclimPoint(specimen_id=spe1, mada_geom_point=POINT (644890.8921103649 7392153.658976035), sampled=False),
                MadaclimPoint(specimen_id=spe2, mada_geom_point=POINT (610233.867750987 7772846.143786541), sampled=False)
            ]

            >>> # You cannot add duplicates, each point must have a unique specimen_id attribute.
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
                mada_geom_point = POINT (644890.8921103649 7392153.658976035),
                sampled_data = None,
                nodata_layers = None
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

            >>> # Remove points by passing in the MadaclimPoint objects, the index or the specimen_id to remove from the Collection.
            
            >>> # From MadaclimPoint instances
            >>> sample_W = collection.all_points[0]
            >>> sample_W
            MadaclimPoint(
                specimen_id = sample_W,
                source_crs = 4326,
                latitude = -16.295741,
                longitude = 46.826763,
                mada_geom_point = POINT (695186.2170220022 8197477.647690434),
                sampled_data = None,
                nodata_layers = None
            )

            >>> # You must specify madaclim_points keyword arg
            >>> collection.remove_points(sample_W)
            Traceback (most recent call last):
            File "<stdin>", line 1, in <module>
            TypeError: MadaclimCollection.remove_points() takes 1 positional argument but 2 were given
            >>> collection.remove_points(madaclim_points=sample_W)
            >>> collection
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_X, mada_geom_point=POINT (955230.600222457 8005985.896187438), sampled=False),
                MadaclimPoint(specimen_id=sample_Y, mada_geom_point=POINT (838822.9378705097 7903464.491492193), sampled=False),
                MadaclimPoint(specimen_id=sample_Z, mada_geom_point=POINT (1036778.1471182993 8526563.721996231), sampled=False)
            ]

            >>> # Using the position index of the instance
            >>> collection = MadaclimCollection.populate_from_df(sample_df)
            Creating MadaclimPoint(specimen_id=sample_W...)
            Creating MadaclimPoint(specimen_id=sample_X...)
            Creating MadaclimPoint(specimen_id=sample_Y...)
            Creating MadaclimPoint(specimen_id=sample_Z...)
            Created new MadaclimCollection with 4 samples.
            >>> collection
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_W, mada_geom_point=POINT (695186.2170220022 8197477.647690434), sampled=False),
                MadaclimPoint(specimen_id=sample_X, mada_geom_point=POINT (955230.600222457 8005985.896187438), sampled=False),
                MadaclimPoint(specimen_id=sample_Y, mada_geom_point=POINT (838822.9378705097 7903464.491492193), sampled=False),
                MadaclimPoint(specimen_id=sample_Z, mada_geom_point=POINT (1036778.1471182993 8526563.721996231), sampled=False)
            ]

            >>> collection.remove_points(indices=-1)    # Removes last point of the collection
            >>> collection
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_W, mada_geom_point=POINT (695186.2170220022 8197477.647690434), sampled=False),
                MadaclimPoint(specimen_id=sample_X, mada_geom_point=POINT (955230.600222457 8005985.896187438), sampled=False),
                MadaclimPoint(specimen_id=sample_Y, mada_geom_point=POINT (838822.9378705097 7903464.491492193), sampled=False)
            ]

            >>> # Using the specimen.id attribute
            >>> collection.remove_points(madaclim_points="sample_Y")
            >>> collection
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_W, mada_geom_point=POINT (695186.2170220022 8197477.647690434), sampled=False),
                MadaclimPoint(specimen_id=sample_X, mada_geom_point=POINT (955230.600222457 8005985.896187438), sampled=False),
                MadaclimPoint(specimen_id=sample_Z, mada_geom_point=POINT (1036778.1471182993 8526563.721996231), sampled=False)
            ]

            >>> # We can also use a list to remove multiple points at the same time.
            >>> # A list of str or MadaclimPoint or mixed types are accepted for the madaclim_points argument.
            >>> sample_w = collection.all_points[0]
            >>> to_remove = [sample_w, "sample_X"]
            >>> collection.remove_points(madaclim_points=to_remove)
            >>> collection
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_Y, mada_geom_point=POINT (838822.9378705097 7903464.491492193), sampled=False),
                MadaclimPoint(specimen_id=sample_Z, mada_geom_point=POINT (1036778.1471182993 8526563.721996231), sampled=False)
            ]

            >>> # Or pass in a list of indices to the indices argument.
            >>> collection.remove_points(indices=[0, -1])    # Remove first and last point
            >>> collection
            MadaclimCollection = [
                MadaclimPoint(specimen_id=sample_X, mada_geom_point=POINT (955230.600222457 8005985.896187438), sampled=False),
                MadaclimPoint(specimen_id=sample_Y, mada_geom_point=POINT (838822.9378705097 7903464.491492193), sampled=False)
            ]

            >>> # Finaly we can clear the collection of all instances.
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
            This method also updates the 'sampled_raster_data' and 'nodata_layers' attributes of the MadaclimCollection instance.

        Examples:
            >>> # Start with a collection
            >>> from coffeaphylogeo.geoclim.raster_manipulation import MadaclimPoint, MadaclimCollection
            >>> specimen_1 = MadaclimPoint(specimen_id="spe1_aren", latitude=-18.9333, longitude=48.2, genus="Coffea", species="arenesiana", has_sequencing=True)
            >>> specimen_2 = MadaclimPoint(specimen_id="spe2_humb", latitude=-12.716667, longitude=45.066667, source_crs=4326, genus="Coffea", species="humblotiana", has_sequencing=True)
            >>> collection = MadaclimCollection()
            >>> collection.add_points([specimen_1, specimen_2])

            >>> # Fetch a specific set of layers to sample(using the MadaclimLayers class utilities)
            >>> from coffeaphylogeo.geoclim.madaclim_info import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> bioclim_labels = [label for label in madaclim_info.get_layers_labels(as_descriptive_labels=True) if "bio" in label]
            >>> bio1 = bioclim_labels[0]
            >>> bio1
            'clim_37_bio1 (Annual mean temperature)'

            >>> # Sample the current_climate raster
            >>> collection    # sampled status set to False
            MadaclimCollection = [
                    MadaclimPoint(specimen_id=spe1_aren, mada_geom_point=POINT (837072.9150244407 7903496.320897499),sampled=False),
                    MadaclimPoint(specimen_id=spe2_humb, mada_geom_point=POINT (507237.57495924993 8594195.741515966),sampled=False)
            ]
            >>> collection_bioclim_data = collection.sample_from_rasters(bioclim_labels)

            ######################################## Extracting data for: spe1_aren ########################################

            Sampling 19 layer(s) from madaclim_current.tif (geoclim_type=clim)...
            Extracting layer 55: bio19 (Precipitation of coldest quarter):  100%|████████████████████████████████████████████████████████████████████████████████████| layer 19/19 [Time remaining: 00:00]

            Finished raster sampling operation in 7.39 seconds.


            ######################################## Extracting data for: spe2_humb ########################################

            Sampling 19 layer(s) from madaclim_current.tif (geoclim_type=clim)...
            Extracting layer 55: bio19 (Precipitation of coldest quarter):  100%|████████████████████████████████████████████████████████████████████████████████████| layer 19/19 [Time remaining: 00:00]

            Finished raster sampling operation in 7.29 seconds.

            >>> collection    # sampled status updated
            MadaclimCollection = [
                    MadaclimPoint(specimen_id=spe1_aren, mada_geom_point=POINT (837072.9150244407 7903496.320897499),sampled=True),
                    MadaclimPoint(specimen_id=spe2_humb, mada_geom_point=POINT (507237.57495924993 8594195.741515966),sampled=True)
            ]
            >>> # Dictionary output with keys based on MadaclimPoint.specimen_id attribute
            >>> list(collection_bioclim_data.keys())
            ['spe1_aren', 'spe2_humb']
            >>> collection_bioclim_data["spe2_humb"]["layer_55"]
            66
            >>> # Results also stored in the .sampled_raster_data attribute
            >>> collection.sampled_raster_data["spe2_humb"]["layer_55"]
            66

            >>> # Sample all layers and examine nodata layers with more informative layers names
            >>> collection_all_layers, collection_nodata_layers = collection.sample_from_rasters(layer_info=True, return_nodata_layers=True)

            ######################################## Extracting data for: spe1_aren ########################################

            Sampling 70 layer(s) from madaclim_current.tif (geoclim_type=clim)...
            Extracting layer 70: ndm (Number of dry months in the year):  100%|██████████████████████████████████████████████████████████████████████████████████████| layer 70/70 [Time remaining: 00:00]

            Sampling 9 layer(s) from madaclim_enviro.tif (geoclim_type=env)...
            Extracting layer 79: forestcover (None):  100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████| layer 9/9 [Time remaining: 00:00]

            Finished raster sampling operation in 25.37 seconds.


            ######################################## Extracting data for: spe2_humb ########################################

            Sampling 70 layer(s) from madaclim_current.tif (geoclim_type=clim)...
            Extracting layer 70: ndm (Number of dry months in the year):  100%|██████████████████████████████████████████████████████████████████████████████████████| layer 70/70 [Time remaining: 00:00]

            Sampling 9 layer(s) from madaclim_enviro.tif (geoclim_type=env)...
            Extracting layer 79: forestcover (None):  100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████| layer 9/9 [Time remaining: 00:00]
            BEWARE! 5 layer(s) contain a nodata value at the specimen location

            Finished raster sampling operation in 25.22 seconds.

            >>> len(collection_nodata_layers)
            1
            >>> len(collection_nodata_layers["spe2_humb"])
            5
            >>> collection_nodata_layers["spe2_humb"][-1]
            'env_79_forestcover (None)'
            >>> # Also saved in the 'nodata_layers' attribute
            >>> collection.nodata_layers["spe2_humb"]
            ['env_75_geology (1=Alluvial_&_Lake_deposits, 2=Unconsolidated_Sands, 4=Mangrove_Swamp, 5=Tertiary_Limestones_+_Marls_&_Chalks, 6=Sandstones, 7=Mesozoic_Limestones_+_Marls_(inc._"Tsingy"), 9=Lavas_(including_Basalts_&_Gabbros), 10=Basement_Rocks_(Ign_&_Met), 11=Ultrabasics, 12=Quartzites, 13=Marble_(Cipolin))', 'env_76_soil (None)', 'env_77_vegetation (None)', 'env_78_watersheds (None)', 'env_79_forestcover (None)']


            >>> # layers_to_sample also accepts a single layer
            >>> collection.sample_from_rasters(37)

            {'spe1_aren': {'layer_37': 196}, 'spe2_humb': {'layer_37': 238}}
        """
        
        if not len(self.__all_points) > 0:
            raise ValueError("Not MadaclimPoint to sample from in the Collection.")
        
        # # Initialize containers for sampled data
        sampled_data = {}
        nodata_layers = {}

        # Sample rasters for whole collection
        for point in self.__all_points:
            sampled_data_point, nodata_layers_point = point.sample_from_rasters(
                layers_to_sample=layers_to_sample,
                layer_info=layer_info,
                return_nodata_layers=True,
                clim_raster_path=clim_raster_path,
                env_raster_path=env_raster_path
            )

            # Save sampled raster data (nested dicts)
            sampled_data[point.specimen_id] = sampled_data_point

            # Save layers name only when val is nodata
            if len(nodata_layers_point) > 0:
                nodata_layers[point.specimen_id] = nodata_layers_point
                
        # Update instance attributes
        self.__sampled_raster_data = sampled_data
        self.__nodata_layers = nodata_layers if len(nodata_layers) > 0 else None

        if return_nodata_layers:
            return sampled_data, nodata_layers if len(nodata_layers) > 0 else None
        else:
            return sampled_data
        