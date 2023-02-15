import unittest
from geospatial.madaclim import MadaClimLayers

class TestMadaclimLayers(unittest.TestCase):
    
    def test_get_clim_layers(self):
        madalayers = MadaClimLayers()
        other = "test"
        self.assertEqual(madalayers.get_clim_layers(), other)

if __name__ == "__main__":
    mada_test = MadaClimLayers()
    print(mada_test.get_clim_layers())