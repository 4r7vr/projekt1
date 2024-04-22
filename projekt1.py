# - XYZ (geocentryczne) -> BLH (elipsoidalne fi, lambda, h) 
# - BLH -> XYZ 
# - XYZ -> NEUp 
# - BL(GRS80, WGS84, ew. Krasowski) -> 2000 
# - BL(GRS80, WGS84, ew. Krasowski) -> 1992 

from math import sin, cos, sqrt, atan, atan2, degrees, radians

class Transformacje:
   def __init__(self, model):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            e2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """  
        if model == "wgs84":
            self.a = 6378137.0 
            self.b = 6356752.31424518 
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "krasowski":
            self.a = 6378245.000
            self.b = 6356863.019
        else:
            raise NotImplementedError(f"{model} model not implemented")
        
        self.e = sqrt(self.a**2 - self.b**2)/self.a 
        self.e2 = self.e**2
        