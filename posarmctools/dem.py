from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from epsgtools import *

def loadDem( src_filename, withPlot=0 ):

	dataset = gdal.Open(src_filename, gdal.GA_ReadOnly)

	print("Driver: {}/{}".format(dataset.GetDriver().ShortName,
		dataset.GetDriver().LongName))
	print("Size is {} x {} x {}".format(dataset.RasterXSize,
		dataset.RasterYSize,
		dataset.RasterCount))
	print("Projection is {}".format(dataset.GetProjection()))

	# Fetch the coefficients for transforming between 
	# pixel/line (P,L) raster space => projection coordinates (Xp,Yp) space
	# Xp = GT[0] + P*GT[1] + L*GT[2]
	# Yp = GT[3] + P*GT[4] + L*GT[5]

	GT = dataset.GetGeoTransform()
	if GT:
		print("Origin = ({}, {})".format( GT[0], GT[3] ) )
		print("Pixel Size = ({}, {})".format( GT[1], GT[5] ) )

	band = dataset.GetRasterBand(1)
	print("Band Type={}".format(gdal.GetDataTypeName(band.DataType)))

	XSize = band.XSize
	YSize = band.YSize
	dted_elevations = band.ReadAsArray(0, 0, XSize, YSize )

	dted_long = GT[0] + np.arange(XSize) * GT[1]
	dted_lat = GT[3] + np.arange(YSize) * GT[5]

	vmin = 0
	vmax = 300
	left = GT[0]
	right = GT[0] + XSize * GT[1]
	bottom = GT[3] + YSize * GT[5]
	top = GT[3]

	if withPlot:
		plt.figure()
		plt.imshow(dted_elevations, extent=(left, right, bottom, top), vmin=vmin, vmax=vmax, cmap='terrain')
		plt.colorbar()
		plt.grid()
		ax = plt.gca()

	rectBivariateSpline = interp.RectBivariateSpline( dted_lat[::-1], dted_long, dted_elevations[::-1,:] )

	return (dted_long, dted_lat, dted_elevations, rectBivariateSpline)