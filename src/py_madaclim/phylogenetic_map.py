"""
Phylogenetic to Geographic Map Visualization Module

This module provides tools for visualizing phylogenetic trees alongside geographic
data and raster layers, allowing for simultaneous inspection of evolutionary
relationships and spatial distributions.

Classes:
    RasterMetadata: Extract and manage metadata from raster files
    PhylogeneticMap: Main class for creating phylogenetic-geographic maps
"""

import json
import pathlib
from pathlib import Path
from typing import Optional, Union, List, Dict

import pandas as pd
import numpy as np
import rasterio
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
from matplotlib.patches import ConnectionPatch
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from Bio import Phylo

from py_madaclim.info import MadaclimLayers


class RasterMetadata:
    """
    A class to extract and manage metadata from raster files.
    
    Extracts band names and descriptions from raster files, with optional support
    for external JSON metadata files. Useful for labeling plots with meaningful
    band descriptions instead of generic band numbers.
    
    Attributes:
        raster_path (Path): Path to the raster file
        metadata_file (Path): Optional path to JSON metadata file
        bands_info (pd.DataFrame): DataFrame with band information
    """
    
    def __init__(self, raster_path: Union[str, Path], metadata_file: Optional[Union[str, Path]] = None):
        """
        Initialize RasterMetadata by reading the raster file and extracting band information.
        
        Args:
            raster_path (Union[str, Path]): Path to the raster file
            metadata_file (Union[str, Path], optional): Path to JSON metadata file with band descriptions
            
        Raises:
            FileNotFoundError: If raster file does not exist
            IOError: If raster file cannot be opened
        """
        self.raster_path = Path(raster_path)
        self.metadata_file = Path(metadata_file) if metadata_file else None
        self.metadata = self._load_metadata() if self.metadata_file else {}
        self._validate_raster()
        self.bands_info = self._extract_bands_metadata()
    
    def _validate_raster(self):
        """Validate that the raster file exists and is readable."""
        if not self.raster_path.exists():
            raise FileNotFoundError(f"Raster file not found: {self.raster_path}")
        
        try:
            with rasterio.open(self.raster_path) as src:
                pass
        except Exception as e:
            raise IOError(f"Cannot open raster file: {e}")
    
    def _load_metadata(self) -> dict:
        """
        Load band metadata from a JSON file.
        
        Returns:
            dict: Metadata dictionary loaded from JSON file
        """
        if not self.metadata_file or not self.metadata_file.exists():
            return {}
        
        try:
            with open(self.metadata_file, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"Warning: Could not load metadata file: {e}")
            return {}
    
    def _extract_bands_metadata(self) -> pd.DataFrame:
        """
        Extract metadata for all bands in the raster file.
        Uses loaded metadata JSON file if available.
        
        Returns:
            pd.DataFrame: DataFrame with columns:
                - band_number: Band index (1-indexed)
                - band_name: Short band name
                - band_description: Full description
        """
        # Default band mapping
        band_mapping = {
            1: {"name": "alt", "description": "Altitude (meters)"},
            2: {"name": "slo", "description": "Slope (degrees)"},
            3: {"name": "asp", "description": "Aspect; clockwise from North (degrees)"},
            4: {"name": "solrad", "description": "Solar radiation (Wh.m-2.day-1)"},
            5: {"name": "geo", "description": "Rock types (Geology)"},
            6: {"name": "soi", "description": "Soil types"},
            7: {"name": "veg", "description": "Vegetation types"},
            8: {"name": "wat", "description": "Watersheds"},
            9: {"name": "forcov", "description": "Forest cover (%)"},
        }
        
        # Update with loaded metadata if available
        if self.metadata and "table_4" in self.metadata:
            table_4 = self.metadata["table_4"]
            layer_names = table_4.get("layer_name", [])
            layer_descs = table_4.get("layer_description", [])
            
            for i, (name, desc) in enumerate(zip(layer_names, layer_descs)):
                if (i + 1) in band_mapping:
                    band_mapping[i + 1] = {"name": name, "description": desc}
        
        bands_list = []
        
        with rasterio.open(self.raster_path) as src:
            print(f"\nRaster file: {self.raster_path.name}")
            print(f"Total bands: {src.count}")
            print(f"CRS: {src.crs}")
            print(f"Shape: {src.height} x {src.width}")
            print("-" * 80)
            
            for band_idx in range(1, src.count + 1):
                if band_idx in band_mapping:
                    band_name = band_mapping[band_idx]["name"]
                    band_description = band_mapping[band_idx]["description"]
                else:
                    band_name = f"Band {band_idx}"
                    band_description = f"Band {band_idx}"
                
                bands_list.append({
                    'band_number': band_idx,
                    'band_name': band_name,
                    'band_description': band_description,
                })
                
                print(f"Band {band_idx}: {band_name:8} | {band_description}")
        
        print("-" * 80)
        df = pd.DataFrame(bands_list)
        return df
    
    def get_band_name(self, band_number: int) -> str:
        """
        Get the short name of a specific band.
        
        Args:
            band_number (int): Band number (1-indexed)
            
        Returns:
            str: Band short name (e.g., 'alt', 'slo')
            
        Raises:
            ValueError: If band number is out of range
        """
        if band_number < 1 or band_number > len(self.bands_info):
            raise ValueError(f"Band {band_number} is out of range (1-{len(self.bands_info)})")
        
        return self.bands_info.iloc[band_number - 1]['band_name']
    
    def get_band_description(self, band_number: int) -> str:
        """
        Get the full description of a specific band.
        
        Args:
            band_number (int): Band number (1-indexed)
            
        Returns:
            str: Band full description (e.g., 'Altitude (meters)')
            
        Raises:
            ValueError: If band number is out of range
        """
        if band_number < 1 or band_number > len(self.bands_info):
            raise ValueError(f"Band {band_number} is out of range (1-{len(self.bands_info)})")
        
        return self.bands_info.iloc[band_number - 1]['band_description']
    
    def get_all_bands_info(self) -> pd.DataFrame:
        """
        Get complete metadata for all bands.
        
        Returns:
            pd.DataFrame: DataFrame with band information
        """
        return self.bands_info.copy()
    
    def display_bands(self):
        """Print a formatted table of all bands and their information."""
        print("\n" + "=" * 80)
        print("RASTER BANDS METADATA")
        print("=" * 80)
        for _, row in self.bands_info.iterrows():
            print(f"Band {row['band_number']}: {row['band_name']:8} | {row['band_description']}")
        print("=" * 80 + "\n")


class PhylogeneticMap:
    """
    Create combined phylogenetic-geographic visualizations with raster overlays.
    
    This class combines phylogenetic tree visualization with geographic mapping,
    allowing simultaneous inspection of evolutionary relationships and spatial
    distributions. Raster layers (elevation, climate data, etc.) can be overlaid
    on the geographic map.
    
    Attributes:
        tree (Bio.Phylo.BaseTree.Tree): Phylogenetic tree
        gps_data (pd.DataFrame): GPS coordinates with specimen IDs
        offsets (dict): X-axis offsets for tree node positioning
        raster_path (Path): Path to raster file
        metadata_file (Path): Path to raster metadata file
        env_raster (Path): Path to environmental raster (for MadaclimLayers)
        clim_raster (Path): Path to climate raster (for MadaclimLayers)
    """
    
    def __init__(
        self,
        tree_path: Union[str, Path],
        gps_data: pd.DataFrame,
        offsets: Dict[str, float],
        raster_path: Union[str, Path],
        metadata_file: Optional[Union[str, Path]] = None,
        env_raster: Optional[Union[str, Path]] = None,
        clim_raster: Optional[Union[str, Path]] = None,
    ):
        """
        Initialize PhylogeneticMap.
        
        Args:
            tree_path (Union[str, Path]): Path to Newick format tree file
            gps_data (pd.DataFrame): DataFrame with columns: specimen_id, latitude, longitude, 
                                    and any numeric column for coloring
            offsets (Dict[str, float]): Dictionary mapping node names to x-axis offsets
            raster_path (Union[str, Path]): Path to raster file for geographic map
            metadata_file (Union[str, Path], optional): Path to raster metadata JSON
            env_raster (Union[str, Path], optional): Path to environmental raster (MadaclimLayers)
            clim_raster (Union[str, Path], optional): Path to climate raster (MadaclimLayers)
        """
        self.tree = Phylo.read(tree_path, "newick")
        self.gps_data = gps_data.copy()
        self.offsets = offsets
        
        self.raster_path = Path(raster_path)
        self.metadata_file = Path(metadata_file) if metadata_file else None
        self.env_raster = Path(env_raster) if env_raster else None
        self.clim_raster = Path(clim_raster) if clim_raster else None
        
        # Initialize RasterMetadata
        self.raster_metadata = RasterMetadata(self.raster_path, self.metadata_file)
        
        # Initialize MadaclimLayers if rasters provided
        if self.env_raster and self.clim_raster:
            self.mada_layers = MadaclimLayers(clim_raster=self.clim_raster, env_raster=self.env_raster)
        else:
            self.mada_layers = None
    
    def create_visualization(
        self,
        raster_band: int = 1,
        raster_cmap: str = "terrain",
        extent: Optional[List[float]] = None,
        figsize: tuple = (35, 12),
        title: Optional[str] = None,
        save_path: Optional[Union[str, Path]] = None,
    ):
        """
        Create the combined phylogenetic-geographic visualization.
        
        Args:
            raster_band (int): Which raster band to display (default: 1)
            raster_cmap (str): Colormap for raster (default: "terrain")
            extent (List[float], optional): Map extent [west, east, south, north]
            figsize (tuple): Figure size in inches (default: (35, 12))
            title (str, optional): Figure title
            save_path (Union[str, Path], optional): Path to save the figure
            
        Returns:
            tuple: (fig, axes) - matplotlib figure and axes objects
        """
        # If extent not provided, use GPS data bounds
        if extent is None:
            extent = [
                self.gps_data['longitude'].min() - 0.5,
                self.gps_data['longitude'].max() + 0.5,
                self.gps_data['latitude'].min() - 0.5,
                self.gps_data['latitude'].max() + 0.5,
            ]
        
        # Create figure
        fig = plt.figure(figsize=figsize)
        
        # TODO: Implement the full visualization
        # For now, return basic structure that user can build upon
        ax_tree = fig.add_subplot(121)
        ax_map = fig.add_subplot(122, projection=ccrs.PlateCarree())
        
        print("PhylogeneticMap visualization created.")
        print(f"Raster band {raster_band}: {self.raster_metadata.get_band_description(raster_band)}")
        
        if save_path:
            fig.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"Figure saved to: {save_path}")
        
        return fig, (ax_tree, ax_map)


__all__ = ['RasterMetadata', 'PhylogeneticMap']
