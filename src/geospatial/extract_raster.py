import pyproj
import csv

# List of all EPSG reference codes
EPSG_codes = [int(code) for code in pyproj.get_codes('EPSG', 'CRS')]

class SpecimenGeoPoint:
    """
    #TODO class docstring
    """
    all = []
    """
    list: initialize list to later get all class instances
    """

    def __init__(self, id, x, y, epsg=4326):
        """
        #TODO docstring for constructor
        """
        self._id = id

        # Float check for x,y
        try:
            x = float(x)
            y = float(y)
        except ValueError:
            pass
        if not (isinstance(x, float) or isinstance(x, int)):
            raise TypeError("x coordinate must be float or int")
        elif not (isinstance(y, float) or isinstance(y, int)):
            raise TypeError("y coordinate must be float or int")
        else:
            self._x = x
            self._y = y

        # EPSG code in db as str/int
        try : 
            epsg = int(epsg)
        except ValueError:
            print("If epsg code entered as str, only [0-9] allowed")
        # epsg code validation
        if epsg not in EPSG_codes:
            raise ValueError(
                "Input EPSG code not valid, see"
                "https://pyproj4.github.io/pyproj/stable/api/database.html#pyproj.database.get_codes"
            )
        else:
            self._epsg = epsg

    @classmethod
    def load_csv(cls, csvfile):
        """
        Instantiate SpecimenGeoPoint objects into a list by parsing a csv file containing the attributes.
        
        Parameters
        ----------
        csvfile : .csv
            .csv file containing attributes for SpecimenGeoPoint. Header be included and follow this order : 
            id, epsg, x, y,
        
        Returns
        -------
        specimens : `list` of `SpecimenGeoPoint`
            List of all SpecimenGeoPoint instances created from the csv file containing data. 
        
        Examples
        --------
        >>> from scripts.data_extraction import SpecimenGeoPoint
        >>> from pathlib import Path
        >>> csv_file = Path("./data/cities.csv")
        >>> data = SpecimenGeoPoint.load_csv(csv_file)
        >>> print(data)
        >>> #TODO COMPLETE EXAMPLE
        """

        with open(csvfile, 'r') as file:
            reader = csv.DictReader(file)
            # Loop through csv and instantiate objects with proper types
            specimens = [SpecimenGeoPoint(
                id = row['id'],
                epsg = int(row['epsg']),
                x = float(row['x']),
                y = float(row['y']))
                for row in reader]
        return specimens

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):
        self._id = value

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, value):
        if not isinstance(float(value), float):
            raise TypeError("x coordinate must be a float or int")
        else :
            self._x = value
    
    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, value):
        if not isinstance(float(value), float):
            raise TypeError("y coordinate must be a float or int")
        else :
            self._y = value
    
    @property
    def epsg(self):
        return self._epsg

    @epsg.setter
    def epsg(self, value):
        if not isinstance(value, int):
            raise TypeError("epsg must be an int")

        elif value not in EPSG_codes:
            raise ValueError(
                "Input EPSG code not valid, see"
                "https://pyproj4.github.io/pyproj/stable/api/database.html#pyproj.database.get_codes"
            )
        else:
            self._epsg = value
    
    def __repr__(self) -> str:
        """
        #TODO docstring
        """
        return f"{self.__class__.__name__}(id={self._id}, x={self._x}, y={self._y}, epsg={self._epsg})"

    def get_all(self):
        pass
   
class MadaclimDataPoint(SpecimenGeoPoint):
    """
    #TODO docstring
    """
    def __init__(self, id, x, y, epsg=32738):
        """
        #TODO docstring
        """
        super().__init__(id, x, y, epsg=32738)

    
