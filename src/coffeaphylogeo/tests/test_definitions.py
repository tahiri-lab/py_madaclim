import unittest
from coffeaphylogeo.definitions import Definitions

class TestDefinitions(unittest.TestCase):
    
    def test_config_attribute(self):
        defs = Definitions()
        defs.config_filename = "test_config.yaml"
        self.assertEqual(defs.config_filename, "test_config.yaml")

if __name__ == "__main__":
    unittest.main()