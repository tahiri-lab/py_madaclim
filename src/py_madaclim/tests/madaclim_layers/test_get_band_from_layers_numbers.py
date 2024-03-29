from pathlib import Path
from py_madaclim.info import MadaclimLayers
from time import perf_counter

data_dir = Path.cwd().parents[3] / "data"
mada_info = MadaclimLayers(clim_raster= data_dir / "coffea_example/madaclim_current.tif")

# Time for validation with double I/O raster operations
valid_start = perf_counter()
print("Evaluating time performance for both O(n) methods")
print("starting with ")
for i in range(1, 71):
    print(mada_info._get_band_from_layer_number(i, "clim"))

valid_end = perf_counter()
valid_time = valid_end - valid_start
print(f"With _validate_raster (2X I/O): {valid_time} s")

# Time for halving the I/O operations
valid_start = perf_counter()
for i in range(1, 71):
    # mada_info.get_band_from_layer_number_novalid(i, "clim")
    #! DEPRECATED METHOD
    pass

valid_end = perf_counter()
valid_time = valid_end - valid_start
print(f"NO raster validation (1X I/O): {valid_time} s")

Results = """
With _validate_raster (2X I/O): 27.55771644799097 s
NO raster validation (1X I/O): 26.27205237200542 s
"""

# Performance issue is with calling the .read() function. Because this is static for our MadaclimDB, let's just remove the .read operation and hardcode the number
