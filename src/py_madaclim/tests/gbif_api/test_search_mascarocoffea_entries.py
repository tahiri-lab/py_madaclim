import unittest
from unittest.mock import Mock
from unittest.mock import patch
from py_madaclim.utils.gbif_api import search_mascarocoffea_entries

class TestSearchMascarocofffGbifApi(unittest.TestCase):
    def test_year_range_type(self):
        with self.assertRaises(TypeError):
            search_mascarocoffea_entries(year_range="2020,2021")  # year_range should be a tuple, not a string

    def test_year_range_values(self):
        with self.assertRaises(ValueError):
            search_mascarocoffea_entries(year_range=(2021, 2020))  # first year should be smaller than the second
    
    @patch("requests.get")
    def test_year_range_none(self, mock_get):
        mock_response = mock_get.return_value
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "offset": 0,
            "limit": 300,
            "count": 600,
            "endOfRecords": True,
            "results": ["result"]
        }
        search_mascarocoffea_entries(year_range=None)
        args, kwargs = mock_get.call_args
        self.assertNotIn('year', kwargs['params'])

    @patch("requests.get")
    def test_single_response(self, mock_get):
        mock_response = mock_get.return_value
        mock_response.json.return_value = {
            "offset": 0,
            "limit": 300,
            "count": 600,
            "endOfRecords": True,
            "results": ["result"]
        }

        # Call the function
        result = search_mascarocoffea_entries()

        # Check that the list contains one item
        self.assertEqual(len(result), 1)

        # Check that the item in the list is the expected result
        self.assertEqual(result, ["result"])

if __name__ == "__main__":
    unittest.main()
