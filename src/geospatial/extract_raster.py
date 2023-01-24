import pyproj

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

    def __init__(self, id : str, x, y, epsg=4326):
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

    def __str__(self) -> str:
        return f"SpecimenGeoPoint(id={self._id}, x={self._x}, y={self._y}, epsg={self._epsg})"

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

   
    