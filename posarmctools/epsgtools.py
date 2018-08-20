import mpl_toolkits.basemap.pyproj as pyproj

# Define some common projections using EPSG codes
wgs84 = pyproj.Proj("+init=EPSG:4326") # LatLon with WGS84 datum used by GPS units and Google Earth
epsg3857 = pyproj.Proj("+init=EPSG:3857") # WGS84 / Pseudo Mercator
epsg3395 = pyproj.Proj("+init=EPSG:3395") # WGS84 / World Mercator
epsg3948 = pyproj.Proj("+init=EPSG:3948") # RGF93 / CC48 Projected coordinate system

def wgs84ToEpsg( Long_Lat, epsg3xxx, shift=0, orig=(0,0) ):
    
    if shift:
        epsg = pyproj.transform(wgs84, epsg3xxx, Long_Lat[0], Long_Lat[1])
        epsg = (epsg[0] - orig[0], epsg[1] - orig[1])
    else:
        epsg = pyproj.transform(wgs84, epsg3xxx, Long_Lat[0], Long_Lat[1])
        
    return epsg

def epsgToWgs84( X_Y, epsg3xxx, shift=0, orig=(0,0)):

	if shift:
		X_Y = (X_Y[0] + orig[0], X_Y[1] + orig[1])
		wgs = pyproj.transform(epsg3xxx, wgs84, X_Y[0], X_Y[1])
	else:
		wgs = pyproj.transform(epsg3xxx, wgs84, X_Y[0], X_Y[1])

	return wgs
