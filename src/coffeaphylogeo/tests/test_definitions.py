import unittest
from coffeaphylogeo.definitions import Definitions

class TestDefinitions(unittest.TestCase):
    
    def test_configfilename_attribute_name(self):
        defs = Definitions()
        print(defs.config_filename)
        defs.config_filename = "bioclim_chelsa.yaml"
        print(f"After change : {defs.config_filename}")

        self.assertEqual(defs.config_filename, "bioclim_chelsa.yaml")

    def test_config_reference_with_attribute_change(self):
        defs = Definitions()
        init_keys = defs.config.keys()
        print(init_keys)
        
        # Change the config_filename and test if refs change
        defs.config_filename = "bioclim_chelsa.yaml"
        refactor_keys = defs.config.keys()
        print(refactor_keys)

        self.assertNotEqual(init_keys, refactor_keys)


if __name__ == "__main__":
    unittest.main()