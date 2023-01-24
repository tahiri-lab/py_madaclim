# IMPORT TESTING
from geospatial.extract_raster import SpecimenGeoPoint
from definitions import ROOT_DIR

specimen_a = SpecimenGeoPoint("test",1,2)
specimen_b = SpecimenGeoPoint("test_b",4,6)

print(specimen_a)
print(specimen_b)

print(ROOT_DIR)