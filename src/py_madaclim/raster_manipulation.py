import csv
import pathlib
import inspect
from pathlib import Path
from typing import Optional, Union, List, Optional, Tuple, Dict
import time
from tqdm import tqdm
import re
import warnings

import py_madaclim
from py_madaclim._constants import Constants
from py_madaclim.info import MadaclimLayers

import rasterio
import rasterio.errors
import rasterio.plot
import pyproj
from pyproj import Transformer
import shapely
from shapely.geometry import Point
import geopandas as gpd

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1 import make_axes_locatable

class _PlotConfig:
    """
    A class that handles the additional keyword arguments to pass to the LayerPlotter instances.

    This class provides properties that filter a dictionary of arguments based on specified prefixes.
    The properties return a dictionary that consists of only the key-value pairs where the key starts with the specified prefix, 
    with the prefix removed from the keys.

    Attributes:
        POSSIBLE_PREFIXES (list[str]): The list of prefixes that can be used in the plot_element_args dictionary.
        _plot_element_args (dict): The original dictionary passed during the object initialization.

    Example:
        >>> arg_dict = {
            "imshow_data": 123,
            "imshow_plot": True,
            "cax_color": "red",
            "histplot_bins": 500,
        }
        >>> pltcfg = _PlotConfig(arg_dict)
        >>> imshow_args = pltcfg.imshow_args
        >>> print(imshow_args)
        {'data': 123, 'plot': True}
    """

    POSSIBLE_PREFIXES = ["imshow_", "cax_", "histplot_", "subplots_", "rasterpoint_", "vline_"]

    def __init__(self, plot_element_args: Optional[dict]=None) -> None:
        """
        Args:
            plot_element_args (dict, optional): The dictionary of plot element arguments. 
                The keys should start with one of the prefixes in POSSIBLE_PREFIXES. 
                If not provided, an empty dictionary will be used.
            imshow_args (dict): The dictionary of the imshow arguments stripped of the prefix.
            cax_args (dict): The dictionary of the cax arguments stripped of the prefix.
            histplot_args (dict): The dictionary of the histplot arguments stripped of the prefix.
            subplots_args (dict): The dictionary of the subplots arguments stripped of the prefix.
            rasterpoint_args (dict): The dictionary of the Geopandas.plot arguments stripped of the prefix.
            vline_args (dict): The dictionary of the axvline arguments stripped of the prefix.
        """
        self._plot_element_args = plot_element_args if plot_element_args else {}
        
        self._imshow_args = self._classify_args(prefix="imshow_", args=self._plot_element_args)
        self._cax_args = self._classify_args(prefix="cax_", args=self._plot_element_args)
        self._histplot_args = self._classify_args(prefix="histplot_", args=self._plot_element_args)
        self._subplots_args = self._classify_args(prefix="subplots_", args=self._plot_element_args)
        self._rasterpoint_args = self._classify_args(prefix="rasterpoint_", args=self._plot_element_args)
        self._vline_args = self._classify_args(prefix="vline_", args=self._plot_element_args)

    @property
    def imshow_args(self) -> dict:
        """
        Get the arguments for imshow, stripped of the 'imshow_' prefix.

        Returns:
            dict: A dictionary containing the imshow arguments with 'imshow_' prefix removed.
                Empty if no args with prefix.
        """
        return self._imshow_args
    
    @property
    def cax_args(self) -> dict:
        """
        Get the arguments for cax, stripped of the 'cax_' prefix.

        Returns:
            dict: A dictionary containing the cax arguments with 'cax_' prefix removed.
                Empty if no args with prefix.
        """
        return self._cax_args
    
    @property
    def histplot_args(self) -> dict:
        """
        Get the arguments for histplot, stripped of the 'histplot_' prefix.

        Returns:
            dict: A dictionary containing the histplot arguments with 'histplot_' prefix removed.
                Empty if no args with prefix.
        """
        return self._histplot_args
    
    @property
    def subplots_args(self) -> dict:
        """
        Get the arguments for histplot, stripped of the 'subplots_' prefix.

        Returns:
            dict: A dictionary containing the histplot arguments with 'subplots_' prefix removed.
                Empty if no args with prefix.
        """
        return self._subplots_args
    
    @property
    def rasterpoint_args(self) -> dict:
        """
        Get the arguments for the Geodataframe Point plotting on the raster representation.
        The arguments are stripped of the 'rasterpoint_' prefix.

        Returns:
            dict: A dictionary containing the Geodataframe Point plotting arguments for the raster 
                representation with the 'rasterpoint_' prefix removed. Empty if no args with prefix.
        """
        return self._rasterpoint_args
    
    @property
    def vline_args(self) -> dict:
        """
        Get the arguments for the vertical line on the distribution plot.
        The arguments are stripped of the 'vline_' prefix.

        Returns:
            dict: A dictionary containing the vertical line arguments for the histogram 
                representation with the 'vline_' prefix removed. Empty if no args with prefix.
        """
        return self._vline_args

    def _classify_args(self, prefix: str, args: dict):
        """
        Private helper method that classifies (filters and strips prefixes from) arguments.

        Args:
            prefix (str): The prefix to use for filtering and stripping.
            args (dict): The dictionary of arguments.

        Raises:
            TypeError: If the prefix is not a string or if args is not a dictionary.
            ValueError: If the prefix is not in POSSIBLE_PREFIXES or if a key in args doesn't start with a prefix in POSSIBLE_PREFIXES.

        Returns:
            dict: A dictionary containing the arguments with the specified prefix removed.
        """
        args_dict = args.copy()    # Copy to avoid inplace-changes outside method
        
        
        if not isinstance(prefix, str):
            raise TypeError("'prefix' must be a str.")
        
        if prefix not in self.POSSIBLE_PREFIXES:
            raise ValueError(f"Specified 'prefix' must be one of {self.POSSIBLE_PREFIXES}")
        
        if not isinstance(args_dict, dict):
            raise TypeError("'args' must be a dictionary.")

        classified_args = {}    # Args with specified prefix container
        for k, v in args_dict.items():
            
            # Check if the key starts with the specified prefix
            if k.startswith(prefix):
                if not re.match(r"\w+_", k):
                    raise ValueError(f"'args' must be prefixed by at least 1 alphanum char and an underscore (i.e. 'prefix_argument'): {k}")
                
                arg_prefix = k.split("_", 1)[0] + "_"    # First underscore encounter for args with native underscores
                if arg_prefix not in self.POSSIBLE_PREFIXES:
                    raise ValueError(f"{arg_prefix} (derived from {k}) must be one of {self.POSSIBLE_PREFIXES}")
                if arg_prefix == prefix:    # Prefix elimination
                    new_key = k[len(arg_prefix):]    # strip prefix by slicing the string
                    classified_args[new_key] = v
        
        return classified_args


class _LayerPlotter:
    """
    A helper class for the visualization of layers from the Madaclim database. This class handles the plotting 
    of either categorical or continuous raster data, along with their respective distribution plots.
    
    Attributes:
        layer_num (int): The number of the layer to be plotted.
        madaclim_layers (MadaclimLayers): An instance of the MadaclimLayers class, which contains layers' information 
            and their raster data.
        plot_args (dict): A dictionary containing optional arguments for the plot.
    """
    def __init__(self, layer_num: int, madaclim_layers: py_madaclim.info.MadaclimLayers, plot_args: dict) -> None:
        """
        Initialize _LayerPlotter with a specific layer number, a MadaclimLayers object, and plot arguments.
        
        Args:
            layer_num (int): The number of the layer to be plotted.
            madaclim_layers (MadaclimLayers): An instance of the MadaclimLayers class.
            plot_args (dict): A dictionary containing optional arguments for the plot.
        """
        self.layer_num = layer_num
        self._madaclim_layers = self._validate_madaclim_layers(madaclim_layers)
        self.plot_args = plot_args

    @property
    def madaclim_layers(self) -> py_madaclim.info.MadaclimLayers:
        """
        Gets or sets for the '_madaclim_layers' attribute.

        Args:
            new_instance (MadaclimLayers): New instance of the MadaclimLayers class.

        Returns:
            MadaclimLayers: An instance of the MadaclimLayers class.
        """
        return self._madaclim_layers
    
    @madaclim_layers.setter
    def madaclim_layers(self, new_instance: py_madaclim.info.MadaclimLayers):
        self._validate_madaclim_layers(madaclim_layers=new_instance)
        self._madaclim_layers = new_instance

    def plot_layer(self) -> Tuple[matplotlib.figure.Figure, List[matplotlib.axes.Axes]]:
        """
        Plots the raster map of the specified geoclimatic layer, along with a 
        distribution histogram of the layer values. The distribution can be plotted 
        as categorical or continuous based on the nature of the layer data.

        For a categorical layer, the function creates a bar plot displaying the
        frequency of each category in the layer. The raster map is plotted with 
        different colors representing each category. A custom legend with categorical 
        labels is also created.

        For a continuous layer, the function plots a histogram representing the 
        frequency distribution of the layer's values. The raster map is plotted with
        a color gradient based on the range of the layer values. A colorbar is 
        provided for reference.

        Both types of plots also include a histogram showing the distribution of 
        raster values at 1km resolution, plotted as a bar graph for categorical 
        data, and a histogram with optional kernel density estimation for continuous data.

        The appearance of the plots can be customized using the 'plot_args' attribute (see _PlotConfig).

        Returns:
            fig (matplotlib.figure.Figure): The top-level container for all plot elements.
            axes (List[matplotlib.axes.Axes]): An array containing the Axes objects 
                of the subplots.

        Raises:
            TypeError: If 'subplots_figsize' is not a tuple.
            ValueError: If 'subplots_figsize' does not contain 2 elements.
            TypeError: If elements of 'subplots_figsize' are not integers.
            ValueError: If the keys from `get_categorical_combinations` do not match 
                the raster values for categorical layer data.

        """
        
        # Create an instance of PlotConfig to separate categories of kwargs plot config
        plot_cfg = _PlotConfig(self.plot_args)
        
        # Set defaults for subplots
        subplots_figsize = plot_cfg.subplots_args.pop("figsize", None)
        subplots_nrows = plot_cfg.subplots_args.pop("nrows", 1)
        subplots_ncols = plot_cfg.subplots_args.pop("ncols", 2)
        
        if subplots_figsize:    # Better exception hints than matplotlib
            if not isinstance(subplots_figsize, tuple):
                raise TypeError("'subplots_figsize' must be a tuple.")
            if not len(subplots_figsize) == 2:
                raise ValueError("'subplots_figsize' must be a tuple of 2 elements.")
            for ele in subplots_figsize:
                if not isinstance(ele, int):
                    raise TypeError(f"{ele} is not an int. 'subplots_figsize' is a tuple of integers.")    
        
        # Set default for imshow
        imshow_cmap = plot_cfg.imshow_args.pop("cmap", "inferno")

        # Set defaults cbar
        cax_position = plot_cfg.cax_args.pop("position", "right")
        cax_size = plot_cfg.cax_args.pop("size", "5%")
        cax_pad = plot_cfg.cax_args.pop("pad", 0.10)

        # Set defaults for histplot
        histplot_color = plot_cfg.histplot_args.pop("color", "grey")
        histplot_bins = plot_cfg.histplot_args.pop("bins", "auto")
        histplot_kde = plot_cfg.histplot_args.pop("kde", True)
        histplot_stat = plot_cfg.histplot_args.pop("stat", "percent")
        histplot_line_kws = plot_cfg.histplot_args.pop("line_kws", {"linestyle" : "--"})
        
        # Plotting the raster map with distribution data for the selected layer
        band_num = self._madaclim_layers.get_bandnums_from_layers(self.layer_num)[0]
        layer_info = self._madaclim_layers.fetch_specific_layers(self.layer_num, "layer_description", "is_categorical", "units")
        layer_key = list(layer_info.keys())[0]
        layer_info = layer_info[layer_key]    # Eliminate redundancy

        layer_units = "Categorical values" if layer_info["is_categorical"] else layer_info["units"]    # Simplify xlabel for categorical units
        layer_description = layer_info["layer_description"]

        # Get geoclim type from layer_num for raster IO selection            
        geoclim_types = ["clim", "env"]
        geoclim_type = [
            geotype for geotype in geoclim_types 
            if self.layer_num in self._madaclim_layers.select_geoclim_type_layers(geotype)["layer_number"].values
        ][0]
        chosen_raster = self._madaclim_layers.clim_raster if geoclim_type == "clim" else self._madaclim_layers.env_raster

        with rasterio.open(chosen_raster) as raster:
            band_data = raster.read(band_num, masked=True)

            # Categorical data multi plots
            if layer_info["is_categorical"]:
            
                subplots_figsize = (20, 10) if subplots_figsize is None else subplots_figsize
                fig, axes = plt.subplots(nrows=subplots_nrows, ncols=subplots_ncols, figsize=subplots_figsize)

                # Custom categorical palette for over 20 categories
                categ_combinations = self._madaclim_layers.get_categorical_combinations(self.layer_num)
                categ_combinations = categ_combinations[self._madaclim_layers.get_layers_labels(self.layer_num)[0]]
                categ_combinations = dict(sorted(categ_combinations.items(), key=lambda item: item[0]))   # Get categ dict {<numerical_val>:<category>}
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

                # Calculate extent of raster
                left, bottom, right, top = raster.bounds

                # Plot the raster map
                axes[0].imshow(band_data.squeeze(), cmap=cmap, norm=norm, extent=[left, right, bottom, top], **plot_cfg.imshow_args)
                
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
                subplots_figsize = (12, 6) if subplots_figsize is None else subplots_figsize
                fig, axes = plt.subplots(nrows=subplots_nrows, ncols=subplots_ncols, figsize=subplots_figsize)

                # Calculate extent of raster
                left, bottom, right, top = raster.bounds
                # Raster map with cbar
                imshow_vmin = plot_cfg.imshow_args.pop("vmin", np.nanmin(band_data.squeeze()))
                imshow_vmax = plot_cfg.imshow_args.pop("vmax", np.nanmax(band_data.squeeze()))
                rasterio.plot.show(band_data.squeeze(), ax=axes[0], cmap=imshow_cmap, vmin=imshow_vmin, vmax=imshow_vmax, extent=[left, right, bottom, top], **plot_cfg.imshow_args)    # Use rasterio.plot.show() instead                
                im = axes[0].get_images()[0]  # get the first image
                
                # Colorbar customization
                divider = make_axes_locatable(axes[0])
                cax = divider.append_axes(position=cax_position, size=cax_size, pad=cax_pad, **plot_cfg.cax_args)
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
                    **plot_cfg.histplot_args
                )
                axes[1].lines[0].set_color("black")
                axes[1].set_title("Distribution of raster values at 1km resolution", fontsize=10)
                axes[1].set_xlabel(layer_units)

            fig.suptitle(
                f"Layer {self.layer_num}: {layer_description}",
                fontsize=16,
                weight="bold",
                ha="center"
            )

            fig.tight_layout()

            return fig, axes

    def _validate_madaclim_layers(self, madaclim_layers: py_madaclim.info.MadaclimLayers) -> py_madaclim.info.MadaclimLayers:
        if not isinstance(madaclim_layers, py_madaclim.info.MadaclimLayers):
            raise TypeError("'madaclim_layers' must be of an instance of py_madaclim.info.MadaclimLayers")
        
        if madaclim_layers.clim_raster is None or madaclim_layers.env_raster is None:
            raise ValueError("The MadaclimLayers instance must have defined 'clim_raster' and 'env_raster attributes.")
        
        return madaclim_layers


class MadaclimRasters:
    """
    Handles operations on Madaclim climate and environmental raster files. 
    Also provides a method to visualize the raster layers (map) and distribution of the raster values.

    Attributes:
        clim_raster (pathlib.Path): Path to the climate raster file.
        clim_crs (pyproj.crs.crs.CRS): The CRS derived of the climate raster file.
        clim_nodata_val (float): The nodata value from the climate raster file
        clim_bounds (tuple): The bounds of the climate raster in order of (left, bottom, right, top)
        env_raster (pathlib.Path): Path to the environmental raster file.
        env_crs (pyproj.crs.crs.CRS): The CRS derived of the environmental raster file.
        env_nodata_val (float): The nodata value from the environmental raster file
        env_bounds (tuple): The bounds of the environmental raster in order of (left, bottom, right, top)
    """

    def __init__(self, clim_raster: pathlib.Path, env_raster: pathlib.Path) -> None:
        """
        Initializes a MadaclimRasters object with climate and environmental raster files.

        Args:
            clim_raster (pathlib.Path): Path to the climate raster file.
            env_raster (pathlib.Path): Path to the environmental raster file.
        
        Example:
            >>> from py_madaclim.raster_manipulation import MadaclimRasters
            >>> mada_rasters = MadaclimRasters(clim_raster="madaclim_current.tif", env_raster="madaclim_enviro.tif")
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
        self._clim_raster = self._validate_raster(clim_raster)
        self._env_raster = self._validate_raster(env_raster)
    
    @property
    def clim_raster(self) -> pathlib.Path:
        """
        Retrieves or sets the climate raster file path.

        Args:
            value (pathlib.Path): The new climate raster file path.

        Returns:
            pathlib.Path: The climate raster file path.
        """
        return self._clim_raster
    
    @clim_raster.setter
    def clim_raster(self, value: pathlib.Path) -> None:
        self._clim_raster = self._validate_raster(value)
        
    @property
    def env_raster(self) -> pathlib.Path:
        """
        Retrieves or sets the environmental raster file path.

        Args:
            value (pathlib.Path): The new environmental raster file path.

        Returns:
            pathlib.Path: The environmental raster file path.
        """
        return self._env_raster
    
    @env_raster.setter
    def env_raster(self, value: pathlib.Path) -> None:
        self._env_raster = self._validate_raster(value)
        
    @property
    def clim_crs(self) -> pyproj.crs.crs.CRS:
        """
        Retrieves the Coordinate Reference System (CRS) from the Madaclim climate raster.

        This property opens the raster file and retrieves the CRS in EPSG format. The EPSG code is used to 
        create and return a pyproj CRS object.

        Returns:
            pyproj.crs.crs.CRS: The CRS derived of the climate raster file.
        """
        # Get epsg from clim_raster
        with rasterio.open(self._clim_raster) as clim_raster:
            clim_epsg = clim_raster.crs.to_epsg()  # Get the EPSG code of the CRS
            clim_crs = pyproj.CRS.from_epsg(clim_epsg)  # Create a pyproj CRS object
        return clim_crs
    @property
    def clim_nodata_val(self) -> float:
        """
        Retrieves the nodata value from the Madaclim climate raster.

        This property opens the raster file and retrieves nodata value from the raster.

        Returns:
            float: The nodata value from the climate raster file
        """
        
        with rasterio.open(self._clim_raster) as clim_raster:
            clim_nodata = clim_raster.nodata
        return clim_nodata
    
    @property
    def clim_bounds(self) -> Tuple[float]:
        """
        Retrieves the bounds of the climate raster.

        This property opens the raster file and retrieves bounds values (left, bottom, right, top) from the raster.

        Returns:
            Tuple[float]: The bounds values in order (left, bottom, right, top)
        """
        
        with rasterio.open(self._clim_raster) as clim_raster:
            clim_bounds = clim_raster.bounds
        return clim_bounds
    
    @property
    def env_crs(self) -> pyproj.crs.crs.CRS:
        """
        Retrieves the Coordinate Reference System (CRS) from the Madaclim environmental raster.

        This property opens the raster file and retrieves the CRS in EPSG format. The EPSG code is used to 
        create and return a pyproj CRS object.

        Returns:
            pyproj.crs.crs.CRS: The CRS derived of the environmental raster file.
        """
        # Get epsg from clim_raster
        with rasterio.open(self._env_raster) as env_raster:
            env_epsg = env_raster.crs.to_epsg()  # Get the EPSG code of the CRS
            env_crs = pyproj.CRS.from_epsg(env_epsg)  # Create a pyproj CRS object
        return env_crs
    
    @property
    def env_nodata_val(self) -> float:
        """
        Retrieves the nodata value from the Madaclim environmental raster.

        This property opens the raster file and retrieves nodata value from the raster.

        Returns:
            float: The nodata value from the environmental raster file
        """
        
        with rasterio.open(self._env_raster) as env_raster:
            env_nodata = env_raster.nodata
        return env_nodata
    
    @property
    def env_bounds(self) -> Tuple[float]:
        """
        Retrieves the bounds of the climate raster.

        This property opens the raster file and retrieves bounds values (left, bottom, right, top) from the raster.

        Returns:
            Tuple[float]: The bounds values in order (left, bottom, right, top)
        """
        
        with rasterio.open(self._env_raster) as env_raster:
            env_bounds = env_raster.bounds
        return env_bounds
        
    def __str__(self) -> str:
        info = (
            f"MadaclimRasters(\n\tclim_raster = {self._clim_raster.name},"
            f"\n\tclim_crs = {self.clim_crs},"
            f"\n\tclim_nodata_val = {self.clim_nodata_val}"
            f"\n\tenv_raster = {self._env_raster.name},"
            f"\n\tenv_crs = {self.env_crs},"
            f"\n\tenv_nodata_val = {self.env_nodata_val}\n)"
    )
        return info
    
    def __repr__(self) -> str:
        return self.__str__()
    
    def plot_layer(self, layer: Union[str, int], **kwargs) -> Tuple[matplotlib.figure.Figure, List[matplotlib.axes.Axes]]:
        """
        Method to plot a specific layer from the Madagascan climate/environmental raster datasets. The layer 
        is displayed as a raster map and its distribution is plotted in a histogram. 

        It accepts layer labels in the following formats: `layer_<num>` (e.g. "layer_1") and `<descriptive_layer_label>` 
        (e.g. "annual_mean_temperature"). Alternatively, the layer number can be supplied directly as an integer.
        
        Depending on whether the layer is categorical or continuous, the visualization will be different. For categorical 
        layers, it will display a map using different colors for each category and a legend mapping categories to colors. 
        For continuous layers, it will display a color gradient map with a color bar.
        

        Args:
            layer (Union[str, int]): Layer to plot. Accepts layer numbers as integers, or layer labels in 
                                    descriptive or `layer_<num>` format.
            **kwargs: Additional arguments to customize the subplots, imshow, colorbar and histplot from matplotlib. 
                    Use "subplots_<arg>", "imshow_<arg>", "cax_<arg>", and "histplot_<arg>" formats to customize corresponding 
                    matplotlib/sns arguments.

        Returns:
            fig (matplotlib.figure.Figure): The top-level container for all plot elements.
            axes (List[matplotlib.axes.Axes]): An array containing the Axes objects 
                of the subplots.
        
        Raises:
            TypeError: If 'layer' is not a str or an int.
            ValueError: If 'layer' is not found within the range of layers.
        
        Note:
            This method returns the fig and axes object for further customization when used by other classes.
            It uses the private _PlotConfig and _LayerPlotter utility classes for the checks and vizualisation.

        Example:
            >>> # Extract environmental layers labels
            >>> from py_madaclim.info import MadaclimLayers
            >>> mada_info = MadaclimLayers(clim_raster="madaclim_current.tif", env_raster="madaclim_enviro.tif")
            >>> env_labels = mada_info.get_layers_labels("env", as_descriptive_labels=True)
            >>> >>> env_labels[0]    # Using altitude as our example
            'env_71_altitude (Altitude in meters)'
            
            >>> # Default visualization of the raster map
            >>> from py_madaclim.raster_manipulation import MadaclimRasters
            >>> mada_rasters = MadaclimRasters(clim_raster=mada_info.clim_raster, env_raster=mada_info.env_raster)    # Using common attr btw the instances
            >>> mada_rasters.plot_layer(env_layers_labels[0])

            >>> # Pass in any number of kwargs to the imshow or cax (raster + colorbarax) or histplot for customization
            >>> mada_rasters.plot_layer(env_labels[0], imshow_cmap="terrain", histplot_binwidth=100, histplot_stat="count")

            >>> # Some layers are categorical data so the figure formatting will change (no cbar)
            >>> # Rock types example
            >>> geo_rock_label = next(label for label in env_labels if "geo" in label)
            >>> mada_rasters.plot_layer(geo_rock_label, subplots_figsize=(12, 8))

            >>> # For numerical features with highly skewed distribution, specify vmin or vmax for the raster map
            >>> mada_rasters.plot_layer(env_labels[3], imshow_vmin=6000)

            >>> # To know which are the categorical data, use the MadaclimLayers utilities
            >>> mada_info.categorical_layers    # as df
            >>> mada_info.get_categorical_combinations()    # As dict, default selects all possibilities
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
        
        # Use LayerPlotter instance helper class to handle customization + plot layer
        plotter = _LayerPlotter(layer_num=layer_num, madaclim_layers=madaclim_info, plot_args=kwargs)
        fig, axes = plotter.plot_layer()
        return fig, axes

    def _validate_raster(self, raster: Union[str, pathlib.Path]) -> pathlib.Path:
        """
        Validates given the raster file.

        Args:
            raster (Union[str, pathlib.Path]): The raster file name or path to validate
        
        Returns:
            pathlib.Path: The path to the validated raster file.
        Raises:
            TypeError: If the input value is not a pathlib.Path or str.
            ValueError: If a pathlib.Path object could not be created from the input value.
            FileExistsError: If the file does not exist.
            RasterioIOError: If the raster file could not be opened.
            NotGeoreferencedWarning: If the raster file is not georeferenced.

        """
            # Validate type
        if not isinstance(raster, (pathlib.Path, str)):
            raise TypeError("'env_raster' must be a pathlib.Path object or str.")
        
        # Validate path and file
        try:
            raster = Path(raster)
        except:
            raise ValueError(f"Could not create a pathlib.Path object from {raster}")
            
        # Check if raster file exists
        if not raster.exists():
            raise FileExistsError(f"Could not find raster file: {raster}")
        
        # Save current warning filters
        original_filters = warnings.filters[:]
        
        # Convert NotGeoreferencedWarning to an error
        warnings.filterwarnings("error", category=rasterio.errors.NotGeoreferencedWarning)
        # Catch any RasterIO errors
        try:
            with rasterio.open(raster) as raster_file:
                return raster
        except rasterio.errors.RasterioIOError as e:
            raise OSError(f"Could not open file: {raster.name}.\nError: {e}")
        except rasterio.errors.NotGeoreferencedWarning as e:
            raise UserWarning(f"Raster file {raster.name} is not georeferenced: {e}")
        finally:
            # Restore original warning filters
            warnings.filters = original_filters
        

class MadaclimPoint:
    #TODO IMPLEMENT VIZ METHODS
    """
    A class representing a specimen as a geographic point with a specific coordinate reference system (CRS)
    and additional attributes. The class provides methods for validating the point's coordinates
    and CRS, as well as sampling values from climate and environmental rasters of the Madaclim database.
    
    Attributes:
        specimen_id (str): An identifier for the point.
        latitude (float): The latitude of the point.
        longitude (float): The longitude of the point.
        source_crs (pyproj.crs.crs.CRS): The coordinate reference system of the point.
        mada_geom_point (shapely.geometry.point.Point): A Shapely Point object representing the point projected 
            in the Madaclim rasters' CRS.
        sampled_layers (Optional[Dict[str, int]]): A dictionary containing the layers labels as keys 
            and their values from the sampled raster at the Point's position as int.
            None if data has not been sampled yet.
        nodata_layers (Optional[List[str]]): A list containing the layers labels.
            None if data has not been sampled yet or no layers sampled containined `nodata` values.
        is_categorical_encoded (bool): The state of the `binary_encode_categorical` method. If True, the method has been called 
            and a new set of binary features has been generated.
            Otherwise, either the layers have not been sampled or the categorical features not been encoded.
        encoded_categ
        encoded_categ_layers (Optional[Dict[str, int]]): A dictionary containing the set of 
            binary encoded categorical features.
        gdf (gpd.GeoDataFrame): A Geopandas GeoDataFrame generated from instance attributes 
            and mada_geom_point geometry. Updates along any changes to the instance's attributes.
        
    """

    def __init__(
            self, 
            specimen_id: str, 
            longitude: float, 
            latitude: float, 
            source_crs: pyproj.crs.crs.CRS=pyproj.CRS.from_epsg(4326), 
            **kwargs
        ) -> None:
        """
        Initialize a MadaclimPoint object with the given `specimen_id`, `latitude`, `longitude`, and `source_crs`.
        The coordinates provided should respect the nature of the given source CRS native units' (i.e. degrees WGS84 or meters for EPSG:3857).
        Upon instantiation, a `mada_geom_point` is created in the projection of the rasters of Madaclim db if `source_crs` differs.
        Optionally, provide additional keyword arguments to store as instance attributes.
        A GeoDataFrame is also created and updated taking the object's attribute and using the `mada_geom_point` for geometry.

        Args:
            specimen_id (str): An identifier for the point.
            latitude (float): The latitude of the point.
            longitude (float): The longitude of the point.
            source_crs (pyproj.crs.crs.CRS, optional): The coordinate reference system of the point. Defaults to WGS84 (EPSG:4326).
            **kwargs: Additional keyword arguments to store as instance attributes.
        Examples:
            Create a MadaclimPoint instance with required parameters `specimen_id`, `longitude`, `latitude`
            
            >>> from py_madaclim.raster_manipulation import MadaclimPoint
            >>> specimen_1 = MadaclimPoint(specimen_id="spe1_aren", latitude=-18.9333, longitude=48.2)    # Default CRS of EPSG:4326
            >>> print(specimen_1)
            MadaclimPoint(
                specimen_id = spe1_aren,
                source_crs = 4326,
                longitude = 48.2,
                latitude = -18.9333,
                mada_geom_point = POINT (837072.9150244407 7903496.320897499),
                sampled_layers = None (Not sampled yet),
                nodata_layers = None (Not sampled yet),
                is_categorical_encoded = False,
                gdf.shape = (1, 8)
            )

            Also accepts any other kwargs and saves them as attributes with specific typing

            >>> specimen_1 = MadaclimPoint(specimen_id="spe1_aren", latitude=-18.9333, longitude=48.2, genus="Coffea", species="arenesiana", has_sequencing=True)
            >>> print(specimen_1)
            MadaclimPoint(
                specimen_id = spe1_aren,
                source_crs = 4326,
                longitude = 48.2,
                latitude = -18.9333,
                mada_geom_point = POINT (837072.9150244407 7903496.320897499),
                sampled_layers = None (Not sampled yet),
                nodata_layers = None (Not sampled yet),
                is_categorical_encoded = False,
                genus = Coffea,
                species = arenesiana,
                has_sequencing = 1.0,
                gdf.shape = (1, 11)
            )
            Core and additionnal attributes are save to a GeoDataFrame with `mada_geom_point` for geometry

            >>> specimen1.gdf
            specimen_id	source_crs	longitude	latitude	mada_geom_point	sampled_layers	nodata_layers	is_categorical_encoded	genus	species	has_sequencing
            0	spe1_aren	4326	48.2	-18.9333	POINT (837072.915 7903496.321)	None	None	False	Coffea	arenesiana	1.0
        """
        self.specimen_id = specimen_id
        self._source_crs = self._validate_crs(source_crs)
        self._longitude, self._latitude = self._validate_lonlat(longitude, latitude)
        self._mada_geom_point = self._construct_point(
            latitude=self.latitude, 
            longitude=self.longitude,
            source_crs=self.source_crs
        )
        
        # Store base attributes names/vals
        self.__base_attr = {k:v for k,v in self.__dict__.items()}

        # Sampled/encoding states and values
        self._sampled_layers = None
        self._nodata_layers = None

        # Rasters updated post-sampling
        self.__clim_raster = None
        self.__env_raster = None
        
        self._is_categorical_encoded = False
        self._encoded_categ_layers = None
        self._encoded_categ_labels = None

        # Store any additional keyword arguments as instance attributes
        base_args = self.get_args_names()[0] + [self.get_args_names()[1]]
        additional_args = [key for key in kwargs if key not in base_args]
        for key in additional_args:
            value = kwargs[key]
            # type conversion especially for csv/json data used on construction
            try:
                value = float(value)
            except ValueError:
                pass
            value = value.strip() if isinstance(value, str) else value
            setattr(self, key, value)

        # Construct the GeoPandas DataFrame with all attrs and `mada_geom_point` for geometry
        self._gdf = self._construct_geodataframe()
        self.__initial_attributes = self.__dict__.copy()    # For dealing with newly set attr post-instantiation
    
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
        self._source_crs = self._validate_crs(value)
        self._update_mada_geom_point()
        
    @property
    def longitude(self) -> float:
        """
        Gets or sets the longitude attribute.

        Args:
            value (float): The longitude value of the point.

        Returns:
            float: The current longitude value of the point.
        """
        return self._longitude
    
    @longitude.setter
    def longitude(self, value: float):
        if hasattr(self, "_longitude") and hasattr(self, "_latitude"):
            self._longitude, _ = self._validate_lonlat(longitude=value, latitude=self._latitude)
            self._update_mada_geom_point()

    @property
    def latitude(self) -> float:
        """
        Gets or sets the latitude attribute.

        Args:
            value (float): The latitude value of the point.

        Returns:
            float: The current latitude value of the point.
        """
        return self._latitude
    
    @latitude.setter
    def latitude(self, value: float):
        if hasattr(self, "_longitude") and hasattr(self, "_latitude"):
            _, self._latitude = self._validate_lonlat(longitude=self._longitude, latitude=value)
            self._update_mada_geom_point()
    

    @property
    def mada_geom_point(self) -> shapely.geometry.point.Point:
        return self._mada_geom_point

    @property
    def base_attr(self) -> dict:
        """
        Get the base attributes when constructing the instance.

        Returns:
            dict: A dictionary containing the base attributes names as keys and their values as values.
        """
        return self.__base_attr
    
    @property
    def sampled_layers(self) -> Optional[Dict[str, int]]:
        """
        Get the instance's data obtained from the `sampled_from_rasters` method.

        Returns:
            Optional[Dict[str, int]]: A dictionary containing the layers labels as keys and their values as int.
                Returns None if data has not been sampled yet.
        """
        return self._sampled_layers
    
    @property
    def nodata_layers(self) -> Optional[List[str]]:
        """
        Get the layers labels containing `nodata` values when calling the`sampled_from_rasters` method.

        Returns:
            Optional[List[str]]: A list containing the layers labels (either layer_<num> or more descriptive).
                Returns None if data has not been sampled yet or no layers sampled containined `nodata` values.
        """
        return self._nodata_layers
    
    @property
    def is_categorical_encoded(self) -> bool:
        """
        Get the state of the binary encoding of the categorical layers.

        Returns:
            bool: The state of the `binary_encode_categorical` method. If True, the method has been called and a new set of binary features has been generated.
                Otherwise, either the layers have not been sampled or the categorical features have not been encoded.
        """
        return self._is_categorical_encoded 
    
    @property
    def encoded_categ_layers(self) -> Optional[Dict[str, int]]:
        """
        Get the binary encoded categorical layers values.

        Returns:
            Optional[Dict[str, int]]: A dictionary containing the set of binary encoded categorical features.
                Keys are contain the layer number (or more description) and the categorical feature.
                Values are the binary encoded value for that given category.
        """
        return self._encoded_categ_layers
    
    @property
    def encoded_categ_labels(self) -> Optional[List[str]]:
        """
        Get the labels from the binary encoded categorical layers of the instance.

        Returns:
            Optional[List[str]]: A list containing the set of the labels from the binary encoding 
                of categorical features. None if the categorical layers have not been encoded yet.
        """
        return self._encoded_categ_labels
    
    @property
    def gdf(self) -> gpd.GeoDataFrame:
        """
        Get the GeoPandas DataFrame using `mada_geom_point` as geometry.

        Returns:
            gpd.GeoDataFrame: A Geopandas GeoDataFrame generated from instance attributes and Point geometry.
        """
        return self._gdf

    def __str__(self) -> str:
        # Get the current state of all attributes (remove recursiveness of __base/initial attrs)
        show_attributes = {}
        for k, v in self.__dict__.items():
            
            if k not in [
                "_MadaclimPoint__base_attr", 
                "_MadaclimPoint__initial_attributes",
                "_MadaclimPoint__clim_raster",
                "_MadaclimPoint__env_raster",
                "_encoded_categ_layers",
                "_encoded_categ_labels"
                ]:
                
                if k == "_source_crs":    # Display EPSG code for readability
                    show_attributes[k] = v.to_epsg()
                
                elif k == "_sampled_layers":    # Display number of sampled layers or status of sampling
                    if self._sampled_layers:
                        show_attributes["_len(sampled_layers)"] = f"{len(v)} layer(s)"
                    else:
                        show_attributes["_sampled_layers"] = "None (Not sampled yet)"

                elif k == "_nodata_layers":    # Display number of nodata layers or status of sampling
                    if self._sampled_layers:
                        show_attributes["_len(nodata_layers)"] = f"{len(v)} layer(s)" if self._nodata_layers else "None (0 layers)"
                    else:
                        show_attributes["_nodata_layers"] = "None (Not sampled yet)"

                elif k == "_gdf":    # gdf shape readability
                    show_attributes["_gdf.shape"] = (self._gdf).shape
                
                else:
                    show_attributes[k] = v

        # Pretty format
        show_attr_list = [f"{k.lstrip('_')} = {v}" for k, v in show_attributes.items()]
        show_attr_str = ",\n\t".join(show_attr_list)
        madapoint_obj = f"MadaclimPoint(\n\t{show_attr_str}\n)"
        
        return madapoint_obj

    def __repr__(self) -> str:
        return self.__str__()
    
    def __setattr__(self, name, value):
        """
        Overrides the __setattr__ method to update GeoDataFrame when attributes are set.

        This method overrides the standard `__setattr__` method. When an attribute 
        is set on the instance, it checks if it's a new attribute (not one of 
        the instance's initial attributes). If it's new, the GeoDataFrame 
        (represented by the `_gdf` attribute) is updated to reflect the new attribute. 

        Args:
            name (str): The name of the attribute being set.
            value (Any): The value being assigned to the attribute.

        Raises:
            AttributeError: If the attribute being set is '_initial_attributes', 
                as this attribute is meant to remain constant after object creation.

        """
        super().__setattr__(name, value)  # Call the parent class's __setattr__ first
        if (hasattr(self, "_MadaclimPoint__initial_attributes") and 
            name != "_gdf"):
            # Update the `gdf` attribute with the newly added attribute
            self._update_gdf()
    
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
    
    def sample_from_rasters(
            self,
            clim_raster: pathlib.Path, 
            env_raster: pathlib.Path,
            layers_to_sample: Union[int, str, List[Union[int, str]]]="all", 
            layer_info: bool=False,
        ) -> None:
        """
        Samples geoclimatic data from raster files for specified layers at the location of the instances's 
        lat/lon coordinates from the `mada_geom_point` attribute.

        Calling this method will also update the `sampled_layers` attributes with the data
        extracted from the layers_to_sample. If sampled data containing 'nodata' values,
        the `nodata_layers` attribute will be updated with the name of the layers accordingly.
        Also, the `gdf` attribute GeoDataFrame will be updated with the `sampled_layers`.

        Args:
            clim_raster_path (pathlib.Path): Path to the climate raster file.
            env_raster_path (pathlib.Path): Path to the environment raster file.
            layers_to_sample (Union[int, str, List[Union[int, str]]], optional): The layer number(s) to sample from the raster files.
                Can be a single int, a single string in the format 'layer_<num>', or a list of ints or such strings. Defaults to 'all'.
            layer_info (bool, optional): Whether to use descriptive labels for the returned dictionary keys. Defaults to False.

        Raises:
            TypeError: If the layers_to_sample is not valid, or if the `mada_geom_point` attribute is not a Point object.
            ValueError: If the layer_number is out of range or if the `mada_geom_point` object is empty.

        Returns:
            None

        Examples:
            >>> # Fetching bioclim layers from the MadaclimLayers class
            >>> from py_madaclim.info import MadaclimLayers
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
            #TODO FIX EXAMPLE
            >>> spe1_l68_l71
            {'clim_68_pet (Annual potential evapotranspiration from the Thornthwaite equation (mm))': 891, 'env_71_altitude (Altitude in meters)': 899}

            >>> # Sample all layers with less descriptive layer names
            >>> specimen_2 = MadaclimPoint(specimen_id="spe2_humb", latitude=-12.716667, longitude=45.066667, source_crs=4326, genus="Coffea", species="humblotiana", has_sequencing=True)
            >>> spe2_all_layers = specimen_2.sample_from_rasters("madaclim_current.tif", "madaclim_enviro.tif")

            ######################################## Extracting data for: spe2_humb ########################################

            Sampling 70 layer(s) from madaclim_current.tif (geoclim_type=clim)...
            Extracting layer 70: ndm (Number of dry months in the year):  100%|████| layer 70/70 [Time remaining: 00:00]

            Sampling 9 layer(s) from madaclim_enviro.tif (geoclim_type=env)...
            Extracting layer 79: forestcover (Percentage of forest cover in 1 km by 1 km grid cells):  100%|████| layer 9/9 [Time remaining: 00:00]
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

            >>> # Calling the sample_from_rasters method also updates the 'sampled_layers' and 'nodata_layers' attributes
            >>> specimen_2.sample_from_rasters(
            ...     clim_raster="madaclim_current.tif",
            ...     env_raster="madaclim_enviro.tif",
            ...     layers_to_sample=[37, 75]
            ... )
            MadaclimPoint(
                    specimen_id = spe2_humb,
                    source_crs = 4326,
                    latitude = -12.716667,
                    longitude = 45.066667,
                    mada_geom_point = POINT (507237.57495924993 8594195.741515966),
                    len(sampled_layers) = 2 layer(s),
                    len(nodata_layers) = 1 layer(s),
                    is_categorical_encoded = False
                    genus = Coffea,
                    species = humblotiana,
                    has_sequencing = True,
                    gdf.shape = (1, 10)
            )
            >>> specimen_2.sampled_layers
            {'layer_37': 238, 'layer_75': -32768}
            >>> specimen_2.nodata_layers
            ['layer_75']

            >>> # Calling the `sample_from_rasters` dynamically updates the `gdf` attributes with the extracted data for easier viz and manipulation
            >>> specimen_2.gdf
            specimen_id  source_crs   latitude  longitude  ...      species  has_sequencing  clim_37_bio1_Annual mean temperature (degrees) env_75_geo_Rock types (categ_vals: 1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13)
            0   spe2_humb        4326 -12.716667  45.066667  ...  humblotiana            True                                             238                                             -32768                     

            [1 rows x 13 columns]
            >>> {k:v for k,v in specimen_2.gdf[["sampled_layers", "nodata_layers"]].iterrows()}[0]    # Updated values for sample status attr
            sampled_layers     2
            nodata_layers    1

        """
        # Reset rasters attributes bound to sampled_state
        self._MadaclimPoint__clim_raster = None
        self._MadaclimPoint__env_raster = None

        # Create a MadaclimRasters to validate both rasters
        mada_rasters = MadaclimRasters(clim_raster=clim_raster, env_raster=env_raster)
        if mada_rasters.clim_crs != Constants.MADACLIM_CRS:
            raise ValueError(f"The provided clim_raster's CRS does not corresponds to Madaclim db's expected crs: {Constants.MADACLIM_CRS}")

        if mada_rasters.env_crs != Constants.MADACLIM_CRS:
            raise ValueError(f"The provided env_raster's CRS does not corresponds to Madaclim db's expected crs: {Constants.MADACLIM_CRS}")
        
        # Create a MadaclimLayers instance to get layers labels and validate layers to sample
        madaclim_info = MadaclimLayers(clim_raster=mada_rasters.clim_raster, env_raster=mada_rasters.env_raster)
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
        has_clim_data = True if clim_raster_sample_info["layers"] else False

        env_raster_sample_info["layers"] = [layer_num for layer_num in layers_numbers if layer_num in geoclim_layer_ranges["env"]]
        env_raster_sample_info["bands"] = madaclim_info.get_bandnums_from_layers(env_raster_sample_info["layers"])
        has_env_data = True if env_raster_sample_info["layers"] else False

        # # Validate if mada_geom_point attribute is of Point geom and not empty
        if not isinstance(self._mada_geom_point, shapely.geometry.point.Point):
            raise TypeError("The 'mada_geom_point' attribute must be a shapely.geometry.point.Point object.")
        
        if self._mada_geom_point.is_empty:
            raise ValueError("The 'mada_geom_point' object cannot be empty.")
        
        # Sample climate and env raster on demand
        sampled_layers = {}
        nodata_layers = []
        
        print("\n" + "#" * 40 + f" \033[1mExtracting data for: {self._specimen_id}\033[0m " +"#" * 40)
        start_time = time.perf_counter()
        
        if clim_raster_sample_info["layers"]:
            total_clim_layers = len(clim_raster_sample_info["layers"])
            
            with rasterio.open(mada_rasters.clim_raster) as clim_raster:
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
                        # Get layer metadata for pbar display and label key for sampled_layers
                        layer_metadata = madaclim_info.fetch_specific_layers(layer_num)
                        pbar.set_description(f"Extracting layer {layer_num}: {layer_metadata['layer_description'].values[0]}")
                        pbar.update()
                        
                        # Sample using the `mada_geom_point`` attributes coordinates for the current layer
                        data = list(clim_raster.sample([(self._mada_geom_point.x, self._mada_geom_point.y)], indexes=band_num))[0]
                        
                        # Save extracted data with specified layer info/name
                        layer_label = (
                            madaclim_info.get_layers_labels(layer_num)[0] 
                            if not layer_info 
                            else madaclim_info.get_layers_labels(layer_num, as_descriptive_labels=True)[0]
                        )
                        sampled_layers[layer_label] = data[0]

                        if data[0] == nodata_clim:    # Save layers where nodata at specimen location
                            nodata_layers.append(layer_label)
        
        if env_raster_sample_info["layers"]:
            total_env_layers = len(env_raster_sample_info["layers"])
            
            with rasterio.open(mada_rasters.env_raster) as env_raster:
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
                        # Get layer metadata for pbar display and label key for sampled_layers
                        layer_metadata = madaclim_info.fetch_specific_layers(layer_num)
                        pbar.set_description(f"Extracting layer {layer_num}: {layer_metadata['layer_description'].values[0]}")
                        pbar.update()
                        
                        # Sample using the `mada_geom_point`` attributes coordinates for the current layer
                        data = list(env_raster.sample([(self._mada_geom_point.x, self._mada_geom_point.y)], indexes=band_num))[0]
                        
                        # Save extracted data with specified layer info/name
                        layer_label = (
                            madaclim_info.get_layers_labels(layer_num)[0] 
                            if not layer_info 
                            else madaclim_info.get_layers_labels(layer_num, as_descriptive_labels=True)[0]
                        )
                        sampled_layers[layer_label] = data[0]

                        if data[0] == nodata_env:    # Save layers where nodata at specimen location
                            nodata_layers.append(layer_label)
                                    

        if len(nodata_layers) > 0:    # No raising exception, just warning print
            print(f"BEWARE! {len(nodata_layers)} layer(s) contain a nodata value at the specimen location")

        end_time = time.perf_counter()
        elapsed_time = end_time - start_time    # Total raster sampling time
        print(f"\nFinished raster sampling operation in {elapsed_time:.2f} seconds.\n")
        
        # Reset categorical encoding attributes
        self._is_categorical_encoded = False   
        self._encoded_categ_layers = None
        self._encoded_categ_labels = None

        # Update private rasters attributes for plot_method
        self._MadaclimPoint__clim_raster = mada_rasters.clim_raster
        self._MadaclimPoint__env_raster = mada_rasters.env_raster
        
        # Update related-instance attributes
        self._nodata_layers = nodata_layers if len(nodata_layers) > 0 else None
        nodata_val = nodata_clim if has_clim_data else nodata_env    # OK because clim_nodata == env_nodata
        
        converted_nan_sampled_layers = {    # Replace nodata int with NaN
            layer_label: np.nan if val == nodata_val else val for layer_label, val in sampled_layers.items()
        }
        self._sampled_layers = converted_nan_sampled_layers
        self._update_gdf()

            
    def binary_encode_categorical(self) -> None:
        """
        Binary encodes the categorical layers contained in the `sampled_layers` attribute.

        This function performs binary encoding of categorical layers
        found in the raster data. It uses the `MadaclimLayers` object 
        to get information about possible categorical layers. If no
        categorical layers are found in the data, a ValueError is raised.
        After the encoding, the function updates the respective instance
        attributes for the categorical encoding status, the encoded
        layers and the gdf replacing the categorical columns by the binary encoded features.

        Raises:
            ValueError: If no categorical layers are found in the raster
            data or if the raster data has not been sampled yet.

        Returns:
            None
        """
        if not self._sampled_layers:
            raise ValueError("Raster data have not been sampled yet. Use `sample_from_rasters` method prior.")
        
        # Check if sampled_layers contains categorical data
        madaclim_info = MadaclimLayers()
        possible_categ_labels = madaclim_info.get_categorical_combinations()
        possible_descriptive_categ_labels = madaclim_info.get_categorical_combinations(as_descriptive_keys=True)
        
        sampled_layers_labels = set(self._sampled_layers.keys())
        categ_layers = set(possible_categ_labels.keys()) | set(possible_descriptive_categ_labels.keys())
        intersect_layers = sampled_layers_labels & categ_layers

        if not intersect_layers:
            raise ValueError(
                f"No categorical data to encode in 'sampled_layers'. It must contain at least one categorical layer.\n"
                f"Either one of:\n{possible_categ_labels.keys()}\n Or one of:\n{possible_descriptive_categ_labels.keys()}\n"
                f"See `get_categorical_combinations` from `MadaclimLayers`"
            )

        # Extract categorical data from sampled_layers and encode them into binary features
        encoded_categ = {}

        for layer in sorted(intersect_layers):
            # Remove the category units from the descriptive label
            layer_label = re.sub(r"\(categ_vals:.*?\)", "", layer)
            layer_label = layer_label.strip()    
            value = self._sampled_layers[layer]
            categories_dict = possible_categ_labels.get(layer, {}) or possible_descriptive_categ_labels.get(layer, {})
            
            for dummy, category in categories_dict.items():
                encoded_categ[f"{layer_label}_{category}"] = 1 if dummy == value else 0
                
            # Create an additional 'nodata' category for layers with missing vals
            if self._nodata_layers:
                if layer in self._nodata_layers:
                    encoded_categ[f"{layer_label}__nodata"] = 1   
                else:
                    encoded_categ[f"{layer_label}__nodata"] = 0
            else:
                encoded_categ[f"{layer_label}__nodata"] = 0
        
        encoded_categ = {k: v for k, v in sorted(encoded_categ.items())}

        # Update categorical layers-related attributes
        self._is_categorical_encoded = True
        self._encoded_categ_layers = encoded_categ    # Overridden `settr` will update `gdf`
        self._encoded_categ_labels = list(encoded_categ.keys())

    def plot_on_layer(self, layer: Union[str, int], **kwargs) -> None:
        """
        **kwargs: Additional arguments to customize the subplots, imshow, colorbar, histplot and Point objects from matplotlib. 
                    Use "subplots_<arg>", "imshow_<arg>", "cax_<arg>", and "histplot_<arg>" formats to customize corresponding 
                    the base Raster and histogram plots as matplotlib/sns arguments. 
                    Use "point_<arg>" to customize the Point objects on the raster.
        """
        def layer_name_range_validation(
                clim_raster: Union[str, pathlib.Path],
                env_raster: Union[str, pathlib.Path],
                layer: Union[str, int]
            ) -> py_madaclim.info.MadaclimLayers:
            """
            Validates the input layer before the visualization. Checks for valid types and values
            based on the MadaclimLayers properties.

            Args:
                clim_raster (Union[str, pathlib.Path]): Path to the climate raster file.
                env_raster (Union[str, pathlib.Path]): Path to the environmental raster file.
                layer (Union[str, int]): Layer to validate. Accepts layer numbers as integers, or layer labels in 
                    descriptive or `layer_<num>` format.descriptive or `layer_<num>` format.

            Raises:
                TypeError: If 'layer' is not a str or an int.
                TypeError: If 'layer' is a number and cannot be converted to an int.
                ValueError: If 'layer' is not found within the range of layers.

            Returns:
                py_madaclim.info.MadaclimLayers: An 'MadaclimLayers' instance if layer is valid.
            """
            
            # Fetch metadata and layer nums with a MadaclimLayers instance
            mada_info = MadaclimLayers(clim_raster=clim_raster, env_raster=env_raster)
            all_layers_df = mada_info.all_layers
            
            # Validate layers to sample
            possible_layers_num_format = mada_info.get_layers_labels()
            possible_layers_desc_format = mada_info.get_layers_labels(as_descriptive_labels=True)


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

            return mada_info

        # Validate object's raster-sampled state
        if (self._MadaclimPoint__clim_raster is None or
            self._MadaclimPoint__env_raster is None):
            raise AttributeError(
                "No 'clim_raster' and 'env_raster' found, use 'sample_from_rasters' prior to setup a reference to the raster files location."
            )
        
        if self._sampled_layers is None:
            raise ValueError("No sampled layers to plot.")

        # Layer validation pre-plotting
        mada_rasters = MadaclimRasters(self._MadaclimPoint__clim_raster, self._MadaclimPoint__env_raster)
        mada_info = layer_name_range_validation(mada_rasters.clim_raster, mada_rasters.env_raster, layer)

        layer_number = mada_info.fetch_specific_layers(layer)["layer_number"].values[0]
        simple_layer_label = mada_info.get_layers_labels(layer_number)[0]
        detailed_layer_label = mada_info.get_layers_labels(layer_number, as_descriptive_labels=True)[0]

        if not (simple_layer_label in self._sampled_layers or
            detailed_layer_label in self._sampled_layers):
            raise ValueError(
                f"Input layer '{simple_layer_label}' or '{detailed_layer_label}'\n"
                f"cannot be found within the 'sampled_layers' at this current state. "
                f"Use the 'sample_from_rasters' method to address that prior to 'plot_on_layer'."
                )

        # Extract gdf Point plotting configs + defaults
        plot_cfg = _PlotConfig(plot_element_args=kwargs)
        rasterpoint_args = {
            "color": plot_cfg.rasterpoint_args.pop("color", "black"),
            "markersize": plot_cfg.rasterpoint_args.pop("markersize", 50),
            "marker": plot_cfg.rasterpoint_args.pop("marker", "o"),
            "edgecolor": plot_cfg.rasterpoint_args.pop("edgecolor", "black"),
            "linestyle": plot_cfg.rasterpoint_args.pop("linestyle", "-"),
            "legend": plot_cfg.rasterpoint_args.pop("legend", True),
            "legend_kwds": plot_cfg.rasterpoint_args.pop("legend_kwds", {"loc": (0.1, 0.9)}),
            "customlabel": plot_cfg.rasterpoint_args.pop("customlabel", self._specimen_id)

        }

        # Plot base raster and distribution histograms
        fig, axes = mada_rasters.plot_layer(layer=layer, **kwargs)
        raster_legend = axes[0].get_legend()
        
        # Overlay with mada_geom_point
        gdf = self._gdf.copy()
        gdf.plot(
            ax=axes[0], 
            color=rasterpoint_args["color"], 
            markersize=rasterpoint_args["markersize"], 
            marker=rasterpoint_args["marker"],
            edgecolor=rasterpoint_args["edgecolor"],
            linestyle=rasterpoint_args["linestyle"],
            legend=rasterpoint_args["legend"], 
            legend_kwds=rasterpoint_args["legend_kwds"],
            **plot_cfg.rasterpoint_args
        )

        # Create a custom legend based on the rasterpoint_args
        if rasterpoint_args["legend"]:
            legend_handle = Line2D(
                [0], [0], 
                marker=rasterpoint_args["marker"],
                markersize=rasterpoint_args["markersize"] / 20,
                markeredgecolor=rasterpoint_args["edgecolor"],
                color=rasterpoint_args["color"], 
                linestyle=rasterpoint_args["linestyle"],
                label=rasterpoint_args["customlabel"]
            )
            axes[0].legend(handles=[legend_handle], loc=rasterpoint_args["legend_kwds"]["loc"])

        # Add the raster legend back to the plot in a different location
        if raster_legend is not None:
            axes[0].add_artist(raster_legend)
            raster_legend.set_bbox_to_anchor((1.05, 1))

        # Extract value to plot vline on dist plot        
        if simple_layer_label in self._sampled_layers:
            sampled_layer_val = self._sampled_layers[simple_layer_label]
        elif detailed_layer_label in self._sampled_layers:
            sampled_layer_val = self._sampled_layers[detailed_layer_label]
        else:
            raise ValueError(    # Failsafe
                f"Could not find {simple_layer_label} or {detailed_layer_label} in {self._sampled_layers}"
            )
            
        # default kwargs for vline
        vline_args = {
            "color": plot_cfg.vline_args.pop("color", rasterpoint_args["color"]),
            "linestyle": plot_cfg.vline_args.pop("linestyle", rasterpoint_args["linestyle"]),
            "linewidth": plot_cfg.vline_args.pop("linewidth", rasterpoint_args["markersize"] / 25),
            "label": plot_cfg.vline_args.pop("label", self._specimen_id)
        }
        
        axes[1].axvline(
            x=sampled_layer_val,
            color=vline_args["color"],
            linestyle=vline_args["linestyle"],
            linewidth=vline_args["linewidth"],
            label=vline_args["label"],
            **plot_cfg.vline_args
        )
        axes[1].legend()
                

    def _validate_crs(self, crs) -> pyproj.crs.crs.CRS:
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
    
    def _validate_lonlat(self, longitude: float, latitude: float) -> Tuple[float]:
        """
        Validate if the given longitude and latitude fall within the bounds of the Madaclim rasters.

        If the source CRS differs from the Madaclim rasters' CRS, the provided longitude and latitude are 
        reprojected to the Madaclim rasters' CRS for validation. If the source CRS and the rasters' CRS are the same, 
        the provided coordinates are directly validated against the rasters' bounds.

        Args:
            longitude (float): The longitude to validate.
            latitude (float): The latitude to validate.

        Returns:
            Tuple[float]: T original longitude and latitude values if valid.

        Raises:
            TypeError: If the provided longitude or latitude can't be converted to float.
            ValueError: If the reprojected or original longitude or latitude don't fall within the bounds of the Madaclim rasters.
        """
        try:
            longitude = float(longitude)
        except:
            raise TypeError(f"Could not convert {longitude} to float. Longitude must be a float.")
        
        try:
            latitude = float(latitude)
        except:
            raise TypeError(f"Could not convert {latitude} to float. Latitude must be a float.")
        
        # Define Madaclim rasters' crs and native units bounds
        madarasters_crs = Constants.MADACLIM_CRS
        madarasters_native_bounds = Constants.MADACLIM_RASTERS_BOUNDS    # (298000.0, 7155000.0, 1101000.0, 8683000.0)
        raster_min_x, raster_max_x = madarasters_native_bounds[0], madarasters_native_bounds[2]
        raster_min_y, raster_max_y = madarasters_native_bounds[1], madarasters_native_bounds[3]

        if self._source_crs != madarasters_crs:
            
            # Reproject the lonlat coords to the Madaclim rasters' crs
            mada_transformer = Transformer.from_crs(self.source_crs, madarasters_crs, always_xy=True)
            reproj_lon, reproj_lat = mada_transformer.transform(longitude, latitude)

            # Get madarasters' bounds in the source_crs' projection
            mada_bounds_reproj_source_crs = mada_transformer.transform_bounds(*madarasters_native_bounds, direction="INVERSE")

            if not raster_min_x <= reproj_lon <= raster_max_x:
                raise ValueError(
                    f"{longitude=} is out of bounds of the Madaclim rasters' for {self._specimen_id}.\n"
                    f"Longitude must fall between {mada_bounds_reproj_source_crs[0]:.6f} "
                    f"and {mada_bounds_reproj_source_crs[2]:.6f} (according to `source_crs`)."
                )
            if not raster_min_y <= reproj_lat <= raster_max_y:
                raise ValueError(
                    f"{latitude=} is out of bounds of the Madaclim rasters' for {self._specimen_id}.\n"
                    f"Latitude must fall between {mada_bounds_reproj_source_crs[1]:.6f} "
                    f"and {mada_bounds_reproj_source_crs[3]:.6f} (according to `source_crs`)."
                )
            return longitude, latitude
        
        else:
            if not raster_min_x <= longitude <= raster_max_x:
                raise ValueError(
                    f"{longitude=} is out of bounds of the Madaclim rasters' for {self._specimen_id}.\n"
                    f"Longitude must fall between {raster_min_x} and {raster_max_x}."
                )
            if not raster_min_y <= latitude <= raster_max_y:
                raise ValueError(
                    f"{latitude=} is out of bounds of the Madaclim rasters' for {self._specimen_id}.\n"
                    f"Latitude must fall between {raster_min_y} and {raster_max_y}."
                )
            return longitude, latitude
    
    #!DEPRECATED METHOD
    def _validate_lat(self, latitude) -> float:
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
            raise TypeError(f"Could not convert {latitude} to float. latitude must be a float.")
        
        # Validate max lon according to crs bounds
        if self._source_crs.is_geographic:
            bounds = self._source_crs.area_of_use.bounds
        else:
            # Extract bounds in native units of projection
            transformer = Transformer.from_crs(self._source_crs.geodetic_crs, self._source_crs, always_xy=True)
            bounds = transformer.transform_bounds(*self._source_crs.area_of_use.bounds)
        
        if bounds[0] < bounds[2]:
            min_lat, max_lat = bounds[0], bounds[2]
        else:
            max_lat, min_lat = bounds[0], bounds[2]
        if not min_lat <= latitude <= max_lat:
            raise ValueError(f"{latitude=} is out of bounds for the crs=EPSG:{self._source_crs.to_epsg()}. latitude must be between {min_lat} and {max_lat}")
        
        return latitude
    
    #!DEPRECATED METHOD
    def _validate_lon(self, longitude) -> float:
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
        if self._source_crs.is_geographic:
            bounds = self._source_crs.area_of_use.bounds
        else:
            # Extract bounds in native units of projection
            transformer = Transformer.from_crs(self._source_crs.geodetic_crs, self._source_crs, always_xy=True)
            bounds = transformer.transform_bounds(*self._source_crs.area_of_use.bounds)
        
        if bounds[0] < bounds[2]:
            min_lon, max_lon = bounds[0], bounds[2]
        else:
            max_lon, min_lon = bounds[0], bounds[2]
        if not min_lon <= longitude <= max_lon:
            raise ValueError(f"{longitude=} is out of bounds for the crs=EPSG:{self._source_crs.to_epsg()}. Longitude must be between {min_lon} and {max_lon}")
        
        return longitude
        
    
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
        Update the `mada_geom_point` attribute by reconstructing the point 
        with the current latitude, longitude, and source_crs.
        """
        if hasattr(self, "_latitude") and hasattr(self, "_longitude"):
            self._mada_geom_point = self._construct_point(
                latitude=self._latitude,
                longitude=self._longitude,
                source_crs=self._source_crs
            )

    def _construct_geodataframe(self) -> gpd.geopandas:
        """
        Constructs a GeoPandas DataFrame using the instance's attributes.

        This method constructs a GeoPandas DataFrame using the attributes of 
        the instance, excluding base attributes and initial attributes. The 
        DataFrame uses `mada_geom_point` as its geometry. 
        
        Some of the object's attributes will be modified or ommitted completely
        to improve readability and increase relevancy of the data saved to
        the GeoDataFrame.

        Note:
            This is a private method, indicated by the underscore prefix. 
            It's intended for internal use within the class, not for use 
            by external code.

        Returns:
            gpd.geopandas: A GeoPandas DataFrame constructed using the 
                instance's attributes.

        Example:
            >>> from py_madaclim.raster_manipulation import MadaclimRasters, MadaclimPoint
            >>> mada_rasters = MadaclimRasters("madaclim_current.tif", "madaclim_enviro.tif")
            >>> spe2 = MadaclimPoint(
            ... specimen_id="spe2_humb", 
            ... latitude=-12.716667, 
            ... longitude=45.066667, 
            ... source_crs=4326
            ... genus="Coffea", 
            ... species="humblotiana", 
            ... has_sequencing=True,
            )
            >>> spe2.gdf
            specimen_id  source_crs   latitude  ...  clim_37_bio1_Annual mean temperature (degrees) clim_69_cwd_Annual climatic water deficit (mm)  env_75_geo_Rock types (categ_vals: 1, 2, 4, 5, 6, 7, 9, 10, 11, 12, 13)
            0   spe2_humb        4326 -12.716667  ...                                             238                                            321                                             -32768                      

            [1 rows x 14 columns]
        """
         # Get the current state of all attributes (remove recursiveness of __base/initial attrs)
        point_attributes = {}
        for k, v in self.__dict__.items():
            
            if k not in [
                "_gdf", 
                "_encoded_categ_layers",
                "_encoded_categ_labels",
                "_MadaclimPoint__base_attr", 
                "_MadaclimPoint__initial_attributes",
                "_MadaclimPoint__clim_raster",
                "_MadaclimPoint__env_raster"
            ]:      
                if k == "_source_crs":    # Display EPSG code for readability
                    point_attributes[k.lstrip("_")] = [v.to_epsg()]
                
                elif k == "_sampled_layers":    # Display number of sampled layers or None if not sampled yet
                    point_attributes[k.lstrip("_")] = [len(v) if self._sampled_layers else v]    

                elif k == "_nodata_layers":    # Display number of nodata layers or status of sampling
                    if self._sampled_layers:
                        point_attributes[k.lstrip("_")] = [len(v) if self._nodata_layers else 0]
                    else:
                        point_attributes[k.lstrip("_")] = [v]

                else:    # Remaining custom defined attributes at construction
                    point_attributes[k.lstrip("_")] = [v]
        
        # Construct gdf with raster sampling state/len
        gdf = gpd.GeoDataFrame(data=point_attributes, geometry="mada_geom_point")

        # Add separate cols for sampled layers in gdf
        if self._sampled_layers:
            for layer, data in self._sampled_layers.items():
                gdf[layer] = data
        
        # Replace categorical layers with binary encoding
        if self._is_categorical_encoded and self._encoded_categ_layers:
            # Find categorical layers within the `sampled_layers`
            madaclim_info = MadaclimLayers()
            possible_categ_labels = madaclim_info.get_categorical_combinations()
            possible_descriptive_categ_labels = madaclim_info.get_categorical_combinations(as_descriptive_keys=True)
            
            sampled_layers_labels = set(self._sampled_layers.keys())
            categ_layers = set(possible_categ_labels.keys()) | set(possible_descriptive_categ_labels.keys())
            intersect_layers = sampled_layers_labels & categ_layers
            
            # Replace sampled_layers with binary encoded layers in categ-only
            gdf = gdf.drop(columns=intersect_layers)
            categ_encoded_df = pd.DataFrame([self._encoded_categ_layers])
            gdf = pd.concat([gdf, categ_encoded_df], axis=1)
           
        return gdf
    
    def _update_gdf(self):
        """
        Updates the instance's `_gdf` attribute using `_construct_geodataframe`.

        This method updates the _gdf attribute of the instance by calling 
        the `_construct_geodataframe` method, which constructs a new 
        GeoPandas DataFrame using the instance's attributes.

        Note:
            This is a private method, indicated by the underscore prefix. 
            It's intended for internal use within the class, not for use 
            by external code.
        """
        self._gdf = self._construct_geodataframe()
    
    
class MadaclimCollection:
    
    #TODO DOCSTRINGS CLS

    #TODO CHECK IF SAMPLED_LAYERS + CATEGORICAL ENCODING WORKS
    
    def __init__(self, madaclim_points: Optional[Union[MadaclimPoint, List[MadaclimPoint]]]=None) -> None:
        """
        Instantiate a collection of MadaclimPoint objects. By default, the MadaclimCollection is empty.
        It will populate the collection with a single instance or a list of MadaclimPoint instances by calling the add_points method with the given madaclim_points.

        Args:
            madaclim_points (Optional[Union[MadaclimPoint, List[MadaclimPoint]]], optional): A single MadaclimPoint object or a list of MadaclimPoint objects to be added to the MadaclimCollection. Initialize an empty MadaclimCollection by default (None).
        Examples:
            >>> from py_madaclim.geoclim.raster_manipulation import MadaclimPoint, MadaclimCollection
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
                sampled_layers = None,
                nodata_layers = None
            )

        """
        self._all_points = []
        if madaclim_points:
            self.add_points(madaclim_points)
        
        # Sampled/encoding states and values
        self._sampled_layers = None
        self._nodata_layers = None

        self._is_categorical_encoded = False
        self._encoded_categ_layers = None
        self._encoded_categ_labels = None
        
        # GeoDataframe construction and updating
        self._gdf = self._construct_geodataframe()



    @property
    def all_points(self) -> list:
        """
        Get the `all_points` attribute.

        Corresponds to a list of each object in the collection.

        Returns:
            list: A list of all the MadaclimPoint objects in the MadaclimCollection.

        Examples:
            >>> # MadaclimPoints are stored in the .all_points attributes in a list
            >>> sample_A = MadaclimPoint()
            >>> collection = MadaclimCollection(sample_A)
            #TODO FIX EXAMPLE
            >>> collection.all_points[0]
            MadaclimPoint(
                specimen_id = sample_A,
                source_crs = 4326,
                latitude = -18.9333,
                longitude = 48.2,
                mada_geom_point = POINT (837072.9150244407 7903496.320897499),
                sampled_layers = None,
                nodata_layers = None
            )
        """
        return self._all_points
    
    @property
    def sampled_layers(self) -> Optional[Dict[str, Dict[str, int]]]:
        """
        Get the sampled_layers attribute of the collection generated from the `sampled_from_rasters` method.

        This attribute is a nested dictionary. The outer dictionary uses the MadaclimPoint.specimen_id as keys. 
        The corresponding value for each key is another dictionary, which uses layer_names as keys and sampled values from rasters as values.

        Returns:
            Optional[Dict[str, Dict[str, int]]]: A dictionary with MadaclimPoint.specimen_id as keys and a dictionary of layer_names (str) and sampled values (int) as values.
                None if Collection has not been sampled yet.
        """
        return self._sampled_layers

    @property
    def nodata_layers(self) -> Optional[Dict[str, Union[str, List[str]]]]:
        """
        Get the nodata_layers attribute of the collection generated from the `sampled_from_rasters` method.

        This attribute is a dictionary that contains the MadaclimPoint.specimen_id as keys and the values as the 'nodata_layers' as str or list of str.
        
        Returns:
            Optional[Dict[str, Union[str, List[str]]]]: A dictionary with MadaclimPoint.specimen_id as keys and values of str or list of str of the layers_name with nodata values.
                None if Collection has not been sampled yet or all layers sampled contained valid data.
        """
        return self._nodata_layers
    
    @property
    def is_categorical_encoded(self) -> bool:
        """
        Get the state of the binary encoding of the categorical layers of the collection.

        Returns:
            bool: The state of the `binary_encode_categorical` method. If True, the method has been called 
                and a new set of binary features has been generated for the whole collection.
                Otherwise, either the layers have not been sampled or the categorical features havenot been encoded.
        """
        return self._is_categorical_encoded 
    
    @property
    def encoded_categ_layers(self) -> Optional[Dict[str, Dict[str, int]]]:
        """
        Get the binary encoded categorical layers values.

        Returns:
            Optional[Dict[str, int]]: A nested dictionary containing the set of binary encoded categorical features.
                The outer dictionary uses the MadaclimPoint.specimen_id as keys. 
                The corresponding value for each key is another dictionary where the keys are layer number
                (or a more descriptive label) with the categorical feature.
                Values correspond to the binary encoded value for that given category
        """
        return self._encoded_categ_layers
    
    @property
    def encoded_categ_labels(self) -> Optional[List[str]]:
        """
        Get the labels from the binary encoded categorical layers for the whole collection.

        Returns:
            Optional[List[str]]: A list containing the set of the labels from the binary encoding 
                of categorical features from the whole collection.
                None if the categorical layers have not been encoded yet.
        """
        return self._encoded_categ_labels
    
    @property
    def gdf(self) -> gpd.GeoDataFrame:
        """
        Get the GeoPandas DataFrame of the collection (concat of all points' gdfs).

        Returns:
            gpd.GeoDataFrame: A Geopandas GeoDataFrame generated from instance attributes and Point geometry.
        """
        self._update_gdf()
        return self._gdf
    
    def __str__(self) -> str:
        if len(self._all_points) == 0:
            return "No MadaclimPoint inside the collection."
        else:
            all_points_short = []
            for point in self._all_points:
                point_info = (
                    f"MadaclimPoint(specimen_id={point.specimen_id}, lat={point.latitude:.6f}, lon={point.longitude:.6f}, "
                    f"sampled_layers={True if point.sampled_layers else False}, "
                    f"categ_encoded={point.is_categorical_encoded})"
                )
                all_points_short.append(point_info)

            return "MadaclimCollection = [\n" + "\t" + ",\n\t".join(all_points_short) + "\n]"
        
    def __repr__(self) -> str:
        return self.__str__()
    
    # def __setattr__(self, name, value):
    #     """
    #     Overrides the __setattr__ method to update GeoDataFrame when attributes are set.

    #     This method overrides the standard `__setattr__` method. When an attribute 
    #     is set on the instance, it checks if it's a new attribute (not one of 
    #     the instance's initial attributes). If it's new, the GeoDataFrame 
    #     (represented by the `_gdf` attribute) is updated to reflect the new attribute. 

    #     Args:
    #         name (str): The name of the attribute being set.
    #         value (Any): The value being assigned to the attribute.

    #     Raises:
    #         AttributeError: If the attribute being set is '_initial_attributes', 
    #             as this attribute is meant to remain constant after object creation.

    #     """
    #     super().__setattr__(name, value)  # Call the parent class's __setattr__ first
    #     if hasattr(self, "_MadaclimCollection__initial_attributes") and name != "_gdf":
    #         # Update the `gdf` attribute with the newly added attribute
    #         self._update_gdf()


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
                sampled_layers = None,
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
                sampled_layers = None,
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
                sampled_layers = None,
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

            # Replace incompatible attribute names
            sanitized_col_names = [MadaclimCollection.sanitize_attr_name(col_name) for col_name in col_names]
            csv_data.fieldnames = sanitized_col_names


            # Contruct points with all req+opt attributes from sanitized headers
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
                sampled_layers = None,
                nodata_layers = None
            )

        """
        df_copy = df.copy()
        if not isinstance(df_copy, pd.DataFrame):
            raise TypeError("'df' is not a pd.DataFrame.")
        
        # Get the required + default args to use to construct the MadaclimPoint objects
        madapoint_required_args, madapoint_default_crs_arg = MadaclimPoint.get_args_names()
        madapoint_default_crs_val = MadaclimPoint.get_default_source_crs()

        # Check for required args present in df
        missing_args = [req_arg for req_arg in madapoint_required_args if req_arg not in df_copy.columns]
        if len(missing_args) > 0:
            raise ValueError(f"df is missing the following required args to construct the MadaclimPoint objects:\n{missing_args}")
        
        # Warn for default EPSG if not provided
        if madapoint_default_crs_arg not in df_copy.columns:
            print(f"Warning! No {madapoint_default_crs_arg} column in the df. Using the default value of EPSG:{madapoint_default_crs_val}...")

        # Sanitize col names for attrs + construct points to add to collection
        df_copy = df_copy.rename(columns={col: MadaclimCollection.sanitize_attr_name(col) for col in df_copy.columns})
        points = []

        for _, row in df_copy.iterrows():
            # If source_crs not in row use default val
            if madapoint_default_crs_arg not in row or not row[madapoint_default_crs_arg]:
                row[madapoint_default_crs_arg] = madapoint_default_crs_val

            print(f"Creating MadaclimPoint(specimen_id={row['specimen_id']}...)")
            point = MadaclimPoint(**row)
            points.append(point)
        
        new_collection = cls(points)
        print(f"Created new MadaclimCollection with {len(points)} samples.")
        return new_collection
    
    @staticmethod
    def sanitize_attr_name(attribute: str):
            """Strips incompatible char from attribute names

            Args:
                attribute (str): The attribute header from the csv file

            Returns:
                str: The compatible attribute name.
            """
            return attribute.replace(" ", "_").replace(".", "").replace("(", "").replace(")", "").replace("(", "")

    def add_points(self, madaclim_points: Union[MadaclimPoint, List[MadaclimPoint]]) -> None:
        """
        Adds one or more MadaclimPoint objects to the MadaclimCollection.

        Args:
            madaclim_points (Union[MadaclimPoint, List[MadaclimPoint]]): A single MadaclimPoint object or a list of MadaclimPoint objects to be added to the MadaclimCollection.

        Raises:
            TypeError: If the input is not a MadaclimPoint object or a list of MadaclimPoint objects.
            ValueError: If the input MadaclimPoint(s) is/are already in the MadaclimCollection or if their specimen_id(s) are not unique.
        Examples:
            >>> from py_madaclim.geoclim.raster_manipulation import MadaclimPoint, MadaclimCollection

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
            File ".../coffeaphylogeo/src/py_madaclim/geoclim/raster_manipulation.py", line 1013, in add_points
                MadaclimPoint(specimen_id=spe2, mada_geom_point=POINT (610233.867750987 7772846.143786541))
                ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            ValueError: MadaclimPoint(
                specimen_id = spe1,
                source_crs = 4326,
                latitude = -23.574583,
                longitude = 46.419806,
                mada_geom_point = POINT (644890.8921103649 7392153.658976035),
                sampled_layers = None,
                nodata_layers = None
            ) is already in the current MadaclimCollection instance.

        """
        # Add multiple MadaclimPoint objects
        if isinstance(madaclim_points, list):
            for point in madaclim_points:
                
                if not isinstance(point, MadaclimPoint):
                    raise TypeError(f"{point} is not a MadaclimPoint object. Accepted types are a single MadaclimPoint and a list of MadaclimPoint objects.")
                
                if point in self._all_points:
                    raise ValueError(f"{point} is already in the current MadaclimCollection instance.")
                
                if point.specimen_id in [mada_pt.specimen_id for mada_pt in self._all_points]:    # specimen_id unique id validation
                    raise ValueError(f"specimen_id={point.specimen_id} already exists inside current MadaclimCollection. Every MadaclimPoint must have a unique id.")
                
                self._all_points.append(point) 
        else:
            # Add single MadaclimPoint object
            if not isinstance(madaclim_points, MadaclimPoint):
                raise TypeError("The madaclim_point to add is not a MadaclimPoint object.")
            
            if madaclim_points in self._all_points:
                raise ValueError(f"{madaclim_points} is already in the current MadaclimCollection instance.")
            
            if madaclim_points.specimen_id in [mada_pt.specimen_id for mada_pt in self._all_points]:    # specimen_id unique id validation
                    raise ValueError(f"specimen_id={madaclim_points.specimen_id} already exists inside current MadaclimCollection. Every MadaclimPoint must have a unique id.")
                
            self._all_points.append(madaclim_points)
        # Update collection gdf after point addition
        self._update_gdf()

    def remove_points(
            self, *, 
            madaclim_points: Optional[Union[MadaclimPoint, List[MadaclimPoint]]]=None, 
            indices: Optional[Union[int, List[int]]]=None, 
            clear:bool=False
        ) -> None:
        """
        Removes MadaclimPoint objects from the MadaclimCollection based on specified criteria.

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
                sampled_layers = None,
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
        if not len(self._all_points) > 0:
            raise ValueError("No points to delete since the MadaclimCollection is empty.")
        
        # Drop all points
        if clear:
            if madaclim_points is not None or indices is not None:
                raise ValueError("When using 'clear', do not provide 'madaclim_points' or 'indices'.")
            else:
                self._all_points.clear()
                self._update_gdf()    # Return string when empty points
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
                        if point_to_rm not in self._all_points:
                            raise ValueError(f"{point_to_rm} does not exists in the current MadaclimCollection instance.")
                        madaclim_points_objects.append(point_to_rm)
                    else:
                        if point_to_rm not in [point.specimen_id for point in self._all_points]:
                            raise ValueError(f"{point_to_rm} is not a valid MadaclimPoint.specimen_id for all points of the collection.")
                        madaclim_points_ids.append(point_to_rm)

                self._all_points = [point for point in self._all_points if point not in madaclim_points_objects and point.specimen_id not in madaclim_points_ids]

            # Single object removal
            else:    
                if not isinstance(madaclim_points, (MadaclimPoint, str)):
                    raise TypeError(f"{madaclim_points} is not a MadaclimPoint object or str. Accepted types are a single MadaclimPoint, a str, a list of MadaclimPoint objects or a list of str.")
                
                if isinstance(madaclim_points, MadaclimPoint):
                    if madaclim_points not in self._all_points:
                        raise ValueError(f"{madaclim_points} does not exists in the current MadaclimCollection instance.")
                    madaclim_points_objects.append(madaclim_points)
                else:
                    if madaclim_points not in [point.specimen_id for point in self._all_points]:
                        raise ValueError(f"{madaclim_points} is not a valid MadaclimPoint.specimen_id for all points of the collection.")
                    madaclim_points_ids.append(madaclim_points)

                self._all_points = [point for point in self._all_points if point not in madaclim_points_objects and point.specimen_id not in madaclim_points_ids]

        
        
        # Remove multiple points from an indices list
        if indices is not None:
            if isinstance(indices, list):
                try:
                    indices = [int(index) if index >= 0 else len(self._all_points) + int(index) for index in indices]
                except:
                    raise TypeError("Indices should be integers.")
                for index in indices:
                    if index < 0 or index >= len(self._all_points):
                        raise ValueError(f"Index {index} is out of bounds.")
                self._all_points = [point for i, point in enumerate(self._all_points) if i not in indices]

            else:    # Remove single point
                try:
                    index = int(indices) if indices >= 0 else len(self._all_points) + int(indices)
                except:
                    raise TypeError("Single index must be an integer.")
                
                if index < 0 or index >= len(self._all_points):
                    raise IndexError("Index out of range.")
                
                else:
                    self._all_points.pop(index)
        
        self._update_gdf()    # Update collection geodataframe with removed points


    def sample_from_rasters(
            self,
            clim_raster: pathlib.Path, 
            env_raster: pathlib.Path,
            layers_to_sample: Union[int, str, List[Union[int, str]]]="all", 
            layer_info: bool=False,
        ) -> None:
        
        """
        Samples geoclimatic data from raster files for specified layers at the location of each point belonging to the MadaclimCollection's instance.

        Calling this method will also update the `sampled_layers` attributes with the data
        extracted from the layers_to_sample for every point in the collection. If sampled data containing 'nodata' values,
        the `nodata_layers` attribute will be updated with the name of the layers accordingly.
        Also, the `gdf` attribute GeoDataFrame will be updated with the `sampled_layers`.
        
        Args:
            clim_raster_path (pathlib.Path): Path to the climate raster file.
            env_raster_path (pathlib.Path): Path to the environment raster file.
            layers_to_sample (Union[int, str, List[Union[int, str]]], optional): The layer number(s) to sample from the raster files.
                Can be a single int, a single string in the format 'layer_<num>', or a list of ints or such strings. Defaults to 'all'.
            layer_info (bool, optional): Whether to use descriptive labels for the returned dictionary keys. Defaults to False.
            

        Returns:
            None
            
        Raises:
            ValueError: If the MadaclimCollection doesn't contain any MadaclimPoints.

        Notes:
            This method also updates the 'sampled_layers' and 'nodata_layers' attributes of the MadaclimCollection instance.

        Examples:
        #TODO FIX EXAMPLES
            >>> # Start with a collection
            >>> from py_madaclim.geoclim.raster_manipulation import MadaclimPoint, MadaclimCollection
            >>> specimen_1 = MadaclimPoint(specimen_id="spe1_aren", latitude=-18.9333, longitude=48.2, genus="Coffea", species="arenesiana", has_sequencing=True)
            >>> specimen_2 = MadaclimPoint(specimen_id="spe2_humb", latitude=-12.716667, longitude=45.066667, source_crs=4326, genus="Coffea", species="humblotiana", has_sequencing=True)
            >>> collection = MadaclimCollection()
            >>> collection.add_points([specimen_1, specimen_2])

            >>> # Fetch a specific set of layers to sample(using the MadaclimLayers class utilities)
            >>> from py_madaclim.info import MadaclimLayers
            >>> madaclim_info = MadaclimLayers()
            >>> bioclim_labels = [label for label in madaclim_info.get_layers_labels(as_descriptive_labels=True) if "bio" in label]
            >>> bio1 = bioclim_labels[0]
            >>> bio1
            'clim_37_bio1 (Annual mean temperature)'

            >>> # Validating the rasters
            mada_rasters = MadaclimRasters(clim_raster="madaclim_current.tif", env_raster="madaclim_enviro.tif")

            >>> # Sample the current_climate raster
            >>> collection    # sampled status set to False
            MadaclimCollection = [
                    MadaclimPoint(specimen_id=spe1_aren, mada_geom_point=POINT (837072.9150244407 7903496.320897499),sampled=False),
                    MadaclimPoint(specimen_id=spe2_humb, mada_geom_point=POINT (507237.57495924993 8594195.741515966),sampled=False)
            ]
            >>> collection_bioclim_data = collection.sample_from_rasters(
                    mada_rasters.clim_raster, 
                    mada_rasters.env_raster,
                    layers_to_sample=bioclim_labels
                )

            ######################################## Extracting data for: spe1_aren ########################################

            Sampling 19 layer(s) from madaclim_current.tif (geoclim_type=clim)...
            Extracting layer 55: Precipitation of coldest quarter:  100%|████| layer 19/19 [Time remaining: 00:00]

            Finished raster sampling operation in 0.01 seconds.


            ######################################## Extracting data for: spe2_humb ########################################

            Sampling 19 layer(s) from madaclim_current.tif (geoclim_type=clim)...
            Extracting layer 55: bio19 Precipitation of coldest quarter:  100%|████| layer 19/19 [Time remaining: 00:00]

            Finished raster sampling operation in 0.02 seconds.

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
            >>> # Results also stored in the `sampled_layers` attribute
            >>> collection.sampled_layers["spe2_humb"]["layer_55"]
            66

            >>> # Sample all layers and examine nodata layers with more informative layers names
            >>> collection_all_layers, collection_nodata_layers = collection.sample_from_rasters(layer_info=True, return_nodata_layers=True)

            ######################################## Extracting data for: spe1_aren ########################################

            Sampling 70 layer(s) from madaclim_current.tif (geoclim_type=clim)...
            Extracting layer 70: ndm (Number of dry months in the year):  100%|████| layer 70/70 [Time remaining: 00:00]

            Sampling 9 layer(s) from madaclim_enviro.tif (geoclim_type=env)...
            Extracting layer 79: forestcover (None):  100%|████| layer 9/9 [Time remaining: 00:00]

            Finished raster sampling operation in 0.43 seconds.


            ######################################## Extracting data for: spe2_humb ########################################

            Sampling 70 layer(s) from madaclim_current.tif (geoclim_type=clim)...
            Extracting layer 70: ndm (Number of dry months in the year):  100%|████| layer 70/70 [Time remaining: 00:00]

            Sampling 9 layer(s) from madaclim_enviro.tif (geoclim_type=env)...
            Extracting layer 79: forestcover (None):  100%|████| layer 9/9 [Time remaining: 00:00]
            BEWARE! 5 layer(s) contain a nodata value at the specimen location

            Finished raster sampling operation in 0.41 seconds.

            >>> len(collection_nodata_layers)
            1
            >>> len(collection_nodata_layers["spe2_humb"])
            5
            >>> collection_nodata_layers["spe2_humb"][-1]
            'env_79_forestcover (None)'
            >>> # Also saved in the 'nodata_layers' attribute
            >>> collection.nodata_layers["spe2_humb"]
            ['env_75_geology (1=Alluvial_&_Lake_deposits, 2=Unconsolidated_Sands, 4=Mangrove_Swamp, 5=Tertiary_Limestones_+_Marls_&_Chalks, 6=Sandstones, 7=Mesozoic_Limestones_+_Marls_(inc._"Tsingy"), 9=Lavas_(including_Basalts_&_Gabbros), 10=Basement_Rocks_(Ign_&_Met), 11=Ultrabasics, 12=Quartzites, 13=Marble_(Cipolin))', 'env_76_soil (None)', 'env_77_vegetation (None)', 'env_78_watersheds (None)', 'env_79_forestcover (None)']


            >>> # layers_to_sample also accepts a single layer, or multiple layers as the output from the `get_layers_labels` method in MadaclimLayers
            >>> collection.sample_from_rasters(37)
            {'spe1_aren': {'layer_37': 196}, 'spe2_humb': {'layer_37': 238}}
        """
        
        if not len(self._all_points) > 0:
            raise ValueError("No MadaclimPoint to sample from in the Collection.")

        sampled_layers, nodata_layers = {}, {}

        # Sample rasters for whole collection
        for point in self._all_points:
            point.sample_from_rasters(
                clim_raster=clim_raster,
                env_raster=env_raster,
                layers_to_sample=layers_to_sample,
                layer_info=layer_info,
            )

            # Save sampled raster data (nested dicts)
            sampled_layers[point.specimen_id] = point.sampled_layers

            # Save layers name only when val is nodata
            nodata_layers[point.specimen_id] = point.nodata_layers if point.nodata_layers else None
                
        # Reset categorical encoding
        self._is_categorical_encoded = False    
        self._encoded_categ_layers = None
        self._encoded_categ_labels = None

        # Update instance attributes
        self._sampled_layers = sampled_layers
        self._nodata_layers = nodata_layers if len(nodata_layers) > 0 else None
        self._update_gdf()

    def binary_encode_categorical(self) -> None:
        """
        Binary encodes the categorical layers contained in the `sampled_layers` attribute.

        This function performs binary encoding of categorical layers
        found in the raster data for each of the Point in the collection.
        After the encoding, the function updates the Collection instance's
        attributes for the categorical encoding status, the encoded
        layers and the gdf replacing the categorical columns by the binary encoded features.

        Raises:
            ValueError: If the MadaclimCollection doesn't contain any MadaclimPoints
                or if the raster data has not been sampled yet.

        Returns:
            None

        Notes:
            See the `binary_encode_categorical` method for logic and possible raised exceptions.
        """
        if not len(self._all_points) > 0:
            raise ValueError("No MadaclimPoint to sample from in the Collection.")
        
        if not self._sampled_layers:
            raise ValueError(
                "Raster data for each point of the collection have not been sampled yet"
                "Use `sample_from_rasters` method prior to `binary_encode_categorical`."
            )
        
        encoded_categ = {}
        encoded_labels = []

        for point in self._all_points:
            point.binary_encode_categorical()
            encoded_categ[point.specimen_id] = point.encoded_categ_layers   # Extract label: value pairs
            labels = point.encoded_categ_layers.keys()
            encoded_labels.extend(label for label in labels if label not in encoded_labels)

        # Update categorical layers-related attributes
        self._is_categorical_encoded = True
        self._encoded_categ_layers = encoded_categ
        self._encoded_categ_labels = encoded_labels
        self._update_gdf()
        
    def _construct_geodataframe(self) -> Union[str, gpd.GeoDataFrame]:
        """
        Constructs a GeoDataFrame from the MadaclimPoint objects stored in the MadaclimCollection.

        The method fills in missing columns across the collection to standardize the GeoDataFrame.
        It also reorders the columns based on the presence of sampled layers and their encoding states.
        The method aims to keep a logical display with sampled layers and/or encoded layers at the right end of the GeoDataFrame.

        Returns:
            Union[str, gpd.GeoDataFrame]: A GeoDataFrame that contains all MadaclimPoint data in a standardized format.
            In case there are no MadaclimPoint objects in the collection, a string message is returned instead.

        Raises:
            ValueError: If there are no MadaclimPoint objects in the collection.

        """
        if not len(self._all_points) > 0:
            collection_gdf = "No MadaclimPoint to sample from in the Collection."
        else:
            # Get all unique and keep the order
            all_points_gdf_cols = [list(point.gdf.columns) for point in self._all_points]
            # unique_cols_whole_collection = set().union(*all_points_gdf_cols)
            unique_cols_whole_collection = list(dict.fromkeys(sum(all_points_gdf_cols, [])))

            # Fill-in missing cols accross collection to standardize collection gdf
            points_gdfs = []
            for point in self._all_points:
                point_df = point.gdf.copy()
                missing_cols = set(unique_cols_whole_collection) - set(point_df.columns)
                for col in missing_cols: 
                    point_df[col] = np.nan
                points_gdfs.append(point_df)

            collection_gdf = pd.concat(points_gdfs, ignore_index=True)
            
            # Check in both collection an individual objects for presence sampled_layers and NO categ_encoding state
            if hasattr(self, "_sampled_layers") and hasattr(self, "_nodata_layers"):    # Avoid block-trigger when constructing with points
                has_sampled_layers = self._sampled_layers is not None or any([point.sampled_layers for point in self._all_points])
                has_not_encoded_layers = not self._is_categorical_encoded or not any([point.is_categorical_encoded for point in self._all_points])

                # Reorder the cols to keep logical display + sampled_layers at right end
                if has_sampled_layers and has_not_encoded_layers:
                    all_sampled_layers_cols = [list(point.sampled_layers.keys()) for point in self._all_points if point.sampled_layers]
                    unique_sampled_layers_cols = list(dict.fromkeys(sum(all_sampled_layers_cols, [])))
                    unsampled_cols = [col for col in collection_gdf.columns if col not in unique_sampled_layers_cols]
                    new_cols_order = list(unsampled_cols) + unique_sampled_layers_cols
                    collection_gdf = collection_gdf.loc[:, new_cols_order]

            # Check in both collection an individual objects for sampled_layers and categ_encoding state
            if hasattr(self, "_is_categorical_encoded") and hasattr(self, "_encoded_categ_layers"):    # Avoid block-trigger when constructing with points
                has_sampled_layers = self._sampled_layers is not None or any([point.sampled_layers for point in self._all_points])
                has_encoded_layers = self._is_categorical_encoded or any([point.is_categorical_encoded for point in self._all_points])

                # Reorder the cols to keep logical display + sampled_layers at right end
                if has_sampled_layers and has_encoded_layers:
                    # Remove the categorical keys from the sampled_layers
                    madaclim_info = MadaclimLayers()
                    possible_categ_labels = list(madaclim_info.get_categorical_combinations().keys())
                    possible_descriptive_categ_labels = list(madaclim_info.get_categorical_combinations(as_descriptive_keys=True).keys())
                    
                    all_sampled_layers_cols = [list(point.sampled_layers.keys()) for point in self._all_points if point.sampled_layers]
                    noncateg_sampled_layers_cols = [ele for sublist in all_sampled_layers_cols for ele in sublist if ele not in possible_categ_labels + possible_descriptive_categ_labels]
                    all_categ_layers_cols = [list(point.encoded_categ_layers.keys()) for point in self._all_points if point.is_categorical_encoded]
                    
                    unique_sampled_layers_cols = list(dict.fromkeys(noncateg_sampled_layers_cols, []))
                    unique_encoded_layers_cols = list(dict.fromkeys(sum(all_categ_layers_cols, [])))
                    
                    non_sampled_categ_cols = [col for col in collection_gdf.columns if col not in unique_sampled_layers_cols + unique_encoded_layers_cols]
                    reordered_cols = list(non_sampled_categ_cols) + unique_sampled_layers_cols + unique_encoded_layers_cols
                    collection_gdf = collection_gdf.loc[:, reordered_cols]

        return collection_gdf
    
    def _update_gdf(self) -> None:
        """
        Updates the internal GeoDataFrame of the MadaclimCollection.

        The method sets the internal `_gdf` attribute to the GeoDataFrame constructed by the `_construct_geodataframe` method.
        """
        self._gdf = self._construct_geodataframe()
        