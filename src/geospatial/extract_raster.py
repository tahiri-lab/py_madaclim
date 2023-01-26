import pyproj
import csv
from pyproj import Transformer, CRS

# List of all EPSG reference codes
EPSG_codes = [int(code) for code in pyproj.get_codes('EPSG', 'CRS')]

#TODO def function to calclulate min/max bounds of a given CRS
    
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
        self.id = id
        self.x = x
        self.y = y
        self.epsg = epsg
    
        # Append each instances of SpecimenGeoPoint class to all list when an instance is born
        SpecimenGeoPoint.all.append(self)

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
        # Do int/valid str conversion
        try:
            value = float(value)
        except ValueError:
            print(f"{value} coordinate must be float or int")
        # Failed try
        if not (isinstance(value, float) or isinstance(value, int)):
            raise TypeError(f"Attribute '{value}' coordinate must be float or int")
        else:
            self._x = value
            #TODO IMPLEMENT COODINATE SAFE CHECKING DEPENDING ON EPSG
            #https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.set_crs.html
    
    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, value):
        # Do int/valid str conversion
        try:
            value = float(value)
        except ValueError:
            print(f"{value} coordinate must be float or int")
        # Failed try
        if not (isinstance(value, float) or isinstance(value, int)):
            raise TypeError(f"Attribute '{value}' coordinate must be float or int")
        else:
            self._y = value
            #TODO IMPLEMENT COODINATE SAFE CHECKING DEPENDING ON EPSG
            #https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoDataFrame.set_crs.html
    
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

    def __repr__(self) -> str:
        """
        #TODO docstring
        """
        return f"{self.__class__.__name__}(id={self._id}, x={self._x}, y={self._y}, epsg={self._epsg})"

    def get_all(self):
        pass

    def get_info(self):
        """
        #TODO docstring
        """
        attributes = [
            att for att in dir(self)
            if att != "all"    # Remove all att
            if not att.startswith("__" and "_")     # Filter special methods and setters 
            and not callable(getattr(self, att))    # Filter methods
        ]
        attributes = sorted(attributes)
        return attributes

    def transform_crs(self, epsg_out):
        """
        Transforms the coordinates to the desired EPSG coordinate reference system. 
        Default argument is EPSG:4326 which is the reference for the WorldClim and Chelsa datasets
        
        Parameters
        ----------
        epsg_out : int
            EPSG Geodetic Parameter Dataset code of the coordinate reference system (CRS) (default is 4326)
        #TODO FINISH DOCSTRING
        """
        # Check if code valid
        if epsg_out not in EPSG_codes:
            raise ValueError("Input EPSG code not valid, see https://pyproj4.github.io/pyproj/stable/api/database.html#pyproj.database.get_codes")
        
        # Call transform method from pyproj
        transformer = Transformer.from_crs(self.epsg, CRS.from_epsg(epsg_out), always_xy=True)
        x_out, y_out = transformer.transform(self.x, self.y)
        self.x = x_out
        self.y = y_out
   
class MadaclimDataPoint(SpecimenGeoPoint):
    """
    #TODO docstring
    """
    def __init__(self, id, x, y, epsg=32738):
        """
        #TODO docstring
        """
        super().__init__(id, x, y, epsg=32738)

    def get_climate_attributes(self):
        """
        #TODO docstring
        """
        attributes = [
            att for att in dir(self)
            if att not in ["all", "id", "x", "y", "epsg"]    # Remove all att
            if not att.startswith("__" and "_")     # Filter special methods and setters 
            and not callable(getattr(self, att))    # Filter methods
        ]
        attributes = sorted(attributes)
        return attributes 
        

    
