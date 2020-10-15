import json, os

from osgeo import gdal, osr
import numpy as np
import matplotlib.pyplot as plt

# **CPLErr GDALDataset::SetGeoTransform 	( 	double *  	padfTransform	)**
# Set the affine transformation coefficients. 

# **CPLErr GDALDataset::GetGeoTransform 	( 	double *  	padfTransform	)**
# Fetch the affine transformation coefficients.

# Fetches the coefficients for transforming between pixel/line (P,L) raster space, and projection coordinates (Xp,Yp) space.

# Xp = padfTransform[0] + P \* padfTransform[1] + L * padfTransform[2];
# Yp = padfTransform[3] + P \* padfTransform[4] + L * padfTransform[5];

# In a north up image, padfTransform[1] is the pixel width, and padfTransform[5] is the pixel height.
# The upper left corner of the upper left pixel is at position (padfTransform[0],padfTransform[3]).

# The default transform is (0,1,0,0,0,1) and should be returned even when a CE_Failure error is returned,
# such as for formats that don't support transformation to projection coordinates.
# This method does the same thing as the C GDALGetGeoTransform() function.

class Tiff(object):
    def __init__(self, filename, focusing, imgAbs, cmap='gray'):
        self.filename = filename
        self.focusing = focusing

        tiffName = os.path.join(focusing.out_dir, self.filename)
        
        plt.imsave(tiffName, imgAbs, cmap=cmap, format='tiff')

        self.upperleftx = focusing.scene.X[0][0]
        self.upperlefty = focusing.scene.Y[0][0]

        self.src_ds = gdal.Open(tiffName, gdal.GA_Update)
        self.setGeoTransform(focusing)
        self.setProjection()

    def setGeoTransform(self, focusing):
        # Specify raster location through geotransform array
        # (upperleftx, scalex, skewx, upperlefty, skewy, scaley)
        # Scale = size of one pixel in units of raster projection

        #Xp = padfTransform[0] + P * padfTransform[1] + L * padfTransform[2];
        #Yp = padfTransform[3] + P * padfTransform[4] + L * padfTransform[5];
        #Xp = upperleftx + P * scalex + L * skewx;
        #Yp = upperlefty + P * skewy  + L * scaley;

        if focusing.alongTrack:
            scalex = self.focusing.conf.d_x * (-np.cos(self.focusing.track.theta))
            skewx = self.focusing.conf.d_x  * (-np.sin(self.focusing.track.theta))
            skewy = self.focusing.conf.d_y  * (-np.sin(self.focusing.track.theta))
            scaley = self.focusing.conf.d_y * (+np.cos(self.focusing.track.theta))
            gt = [self.upperleftx, scalex, skewx, self.upperlefty, skewy, scaley]
        else:
            scalex = self.focusing.conf.d_x
            skewx = 0
            skewy = 0
            scaley = self.focusing.conf.d_y
            
            gt = [self.upperleftx, scalex, skewx, self.upperlefty, skewy, scaley]

        # Set location
        self.src_ds.SetGeoTransform(gt)
        self.src_ds.FlushCache()                     # write to disk

    def setProjection(self):
        srs = osr.SpatialReference()            # establish encoding
        srs.ImportFromEPSG(3948)                # RGF93 / CC48 Projected coordinate system
        self.src_ds.SetProjection(srs.ExportToPrettyWkt()) # export coords to file
        self.src_ds.FlushCache()                     # write to disk
