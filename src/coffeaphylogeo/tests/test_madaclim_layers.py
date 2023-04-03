import unittest
from geoclim.madaclim_layers import MadaclimLayers

class TestMadaclimLayers(unittest.TestCase):
    
    def test_empty_instantiation(self):
        madalayers = MadaclimLayers()
        empty_dict = {}
        self.assertEqual(madalayers.madaclim_layers, empty_dict)

if __name__ == "__main__":
    unittest.main()