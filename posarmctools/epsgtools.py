import configparser, importlib, logging, os, re, sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import scipy.interpolate as interp

import pyproj
from osgeo import gdal

logger = logging.getLogger(__name__)

# Define some common projections using EPSG codes
#wgs84     = pyproj.Proj("+init=EPSG:4326")  # LatLon with WGS84 datum used by GPS units and Google Earth
#epsg3857  = pyproj.Proj("+init=EPSG:3857")  # WGS84 / Pseudo Mercator does not seem to keep distances
#epsg3395  = pyproj.Proj("+init=EPSG:3395")  # WGS84 / World Mercator
#epsg3948  = pyproj.Proj("+init=EPSG:3948")  # RGF93 / CC48 Projected coordinate system
#epsg32630 = pyproj.Proj("+init=EPSG:32630") # # WGS 84 / UTM zone 30N
wgs84     = pyproj.Proj("EPSG:4326")  # LatLon with WGS84 datum used by GPS units and Google Earth
epsg3857  = pyproj.Proj("EPSG:3857")  # WGS84 / Pseudo Mercator does not seem to keep distances
epsg3395  = pyproj.Proj("EPSG:3395")  # WGS84 / World Mercator
epsg3948  = pyproj.Proj("EPSG:3948")  # RGF93 / CC48 Projected coordinate system
epsg32630 = pyproj.Proj("EPSG:32630") # # WGS 84 / UTM zone 30N

def plotEpsgWithRemarkablePoints(x, y, log, prefixes, ptsEpsg, epsgStr, day):
	plt.figure()
	for hour, sl in zip(log.loader.hours, log.slices):
		plt.plot( x[sl], y[sl], label=hour )

	ax = plt.gca()
	addRemarkablePoints(ax, prefixes, ptsEpsg)
	plt.xlabel("x")
	plt.ylabel("y")
	plt.grid()
	plt.legend()
	title = f"{epsgStr} {day}"
	plt.title(title)
	plt.gca().set_aspect('equal')

def getxy(logs, epsg3xxx, idx_long, idx_lat):
	x_y = [wgs84LongLatToEpsg( (log[:,idx_long], log[:,idx_lat]), epsg3xxx ) 
	for log in logs]
	return zip(*x_y)

def addRemarkablePoint(ax, pt):
	ax.plot(pt[0], pt[1], 'o', markerfacecolor="None", markeredgecolor = 'r' )

def addRemarkablePoints(ax, prefixes, ptsEpsg):
	cmap = get_cmap("Set1")  # type: matplotlib.colors.ListedColormap
	colors = cmap.colors  # type: list
	for num, prefix in enumerate(prefixes):
		color = colors[num %len(colors)]
		for idx in prefixes[prefix]:
			if idx == prefixes[prefix][0]:
				label = prefix
				try:
					pt = ptsEpsg[f"{prefix}{idx}"]
				except:
					pt = ptsEpsg[prefix]
			else:
				label=""
				pt = ptsEpsg[prefix + f"{idx}"]
			ax.plot( pt[0], pt[1], 'o', color=color, markeredgecolor = 'black', label=label)

def getReferencePointPrefixes(varDict):
	prefixList = set([var.rstrip('0123456789') for var in varDict])
	prefixDict = {key: [] for key in prefixList}
	for var in varDict:
		prefix = var.rstrip('0123456789')
		num = var.lstrip(prefix)
		if num != "":
			prefixDict[prefix].append(int(num))
		else:
			prefixDict[prefix].append(-1)
	for var in prefixDict:
		prefixDict[prefix].sort()
	return prefixDict

def getReferencePoints(filename, epsg3xxx, file="py"):
	if file == "py":
		module = importlib.import_module(filename.split('.py')[0])
		myVars = [var for var in dir(module) if var[:2] != '__' and var[-2:] != '__']
		pts = {var: module.__dict__[var] for var in myVars if not callable(module.__dict__[var])}
		
	elif file == "ini":
		with open(filename) as f:
			config = configparser.ConfigParser()
			config.read_file(f)
			pts = {}
			for key in config["points"]:
				pts[key] = tuple(map(eval, config["points"][key].split(', ')))
	else:
		print(f"error, file type not known: {file}")
	
	ptsLongLat = {pt: pts[pt][::-1] for pt in pts}
	ptsDict = getReferencePointPrefixes(pts)
	ptsEpsg = {pt: wgs84LongLatToEpsg(ptsLongLat[pt], epsg3xxx) for pt in ptsLongLat}

	return pts, ptsEpsg, ptsDict

def getLatLongFromStr(a, b):
    Da = float(a.split("°")[0])
    r = a.split("°")[1]
    Ma = float(r.split("'")[0])
    r = r.split("'")[1]
    Sa = float(r.split("N")[0])
    
    Db = float(b.split("°")[0])
    r = b.split("°")[1]
    Mb = float(r.split("'")[0])
    r = r.split("'")[1]
    Sb = float(r.split("W")[0])
    
    return ( Da + Ma/60 + Sa/3600, -(Db + Mb/60 + Sb/3600) )

def wgs84LongLatToEpsg(Long_Lat, epsg3xxx, shift=0, orig=(0,0)):
	# The axis order may be swapped if the source and destination CRS’s are defined as having the 
	# first coordinate component point in a northerly direction (See PROJ FAQ on axis order). You can 
	# check the axis order with the pyproj.crs.CRS class. If you prefer to keep your axis order as 
	# always x,y, you can use the always_xy option when creating the pyproj.transformer.Transformer.
	# pyproj.crs.CRS("epsg:3948")
	# pyproj.crs.CRS("epsg:4326")
    if shift:
        epsg = pyproj.transform(wgs84, epsg3xxx, Long_Lat[0], Long_Lat[1], always_xy=True)
        epsg = (epsg[0] - orig[0], epsg[1] - orig[1])
    else:
        epsg = pyproj.transform(wgs84, epsg3xxx, Long_Lat[0], Long_Lat[1], always_xy=True)
    return epsg

def epsgToWgs84LongLat(X_Y, epsg3xxx, shift=0, orig=(0,0)):
	if shift:
		X_Y = (X_Y[0] + orig[0], X_Y[1] + orig[1])
		wgs = pyproj.transform(epsg3xxx, wgs84, X_Y[0], X_Y[1], always_xy=True)
	else:
		wgs = pyproj.transform(epsg3xxx, wgs84, X_Y[0], X_Y[1], always_xy=True)
	return wgs

def interpolate_xyz(record, log, logEpsg):
	x = np.interp(record.timestamps_allRamps, log.timestamps, logEpsg[0])
	y = np.interp(record.timestamps_allRamps, log.timestamps, logEpsg[1])
	z = np.interp(record.timestamps_allRamps, log.timestamps, log.alt)
	xyz_a = np.stack((record.rampNumber_allRamps, record.timestamps_allRamps, x, y, z), -1)

	timestamps_allRamps_b = record.timestamps_allRamps + record.parameters.rampPeriod / 2

	x_b = np.interp(timestamps_allRamps_b, log.timestamps, logEpsg[0])
	y_b = np.interp(timestamps_allRamps_b, log.timestamps, logEpsg[1])
	z_b = np.interp(timestamps_allRamps_b, log.timestamps, log.alt)
	xyz_b = np.stack((record.rampNumber_allRamps, timestamps_allRamps_b, x_b, y_b, z_b), -1)

	return xyz_a, xyz_b

def interpolate_lla(record, log):
	Lat = np.interp(record.timestamps_allRamps, log.timestamps, log.lat )
	Long = np.interp(record.timestamps_allRamps, log.timestamps, log.long )
	Alt = np.interp(record.timestamps_allRamps, log.timestamps, log.alt )
	lla_a = np.stack(
		(record.rampNumber_allRamps, record.timestamps_allRamps, Lat, Long, Alt), -1)

	timestamps_allRamps_b = record.timestamps_allRamps + record.parameters.rampPeriod / 2

	Lat_b  = np.interp(timestamps_allRamps_b, log.timestamps, log.lat)
	Long_b = np.interp(timestamps_allRamps_b, log.timestamps, log.long)
	Alt_b  = np.interp(timestamps_allRamps_b, log.timestamps, log.alt)
	lla_b = np.stack(
		(record.rampNumber_allRamps, timestamps_allRamps_b, Lat_b, Long_b, Alt_b), -1)

	return lla_a, lla_b
