"""
py_madaclim: Climate and Environmental Data Tools for Madagascar

A package for working with Madaclim climate and environmental raster data,
including data extraction, manipulation, visualization, and phylogenetic mapping.
"""

# Import existing classes
from py_madaclim.info import MadaclimLayers
from py_madaclim.raster_manipulation import MadaclimRasters, MadaclimCollection

# Import new phylogenetic mapping functionality
from py_madaclim.phylogenetic_map import RasterMetadata, PhylogeneticMap

# Define public API
__all__ = [
    # Existing functionality
    "MadaclimLayers",
    "MadaclimRasters",
    "MadaclimCollection",
    # New phylogenetic mapping functionality
    "RasterMetadata",
    "PhylogeneticMap",
]

__version__ = "1.2.0"  # Update version if needed
__author__ = "Caroline Fortier"
__email__ = "caroline.fortier@outlook.com"
