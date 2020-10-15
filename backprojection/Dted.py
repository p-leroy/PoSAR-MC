import os, pickle

from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp

from posarmctools import epsgtools as epsg

def saveDted(dted, out_dir):
    filename = os.path.join(out_dir, 'dted')
    with open(filename, 'wb') as file:
        pickle.dump(dted, file)

def loadDted(out_dir):
    filename = os.path.join(out_dir, 'dted')
    with open(filename, 'wb') as file:
        dted = pickle.load(file)
    return dted

class Dted(object):
    def __init__(self, filename):
        self.filename = filename
        print("get dataset")
        self.ds = self.__getDataset()

        self.band = self.ds.GetRasterBand(1)
        self.XSize = self.band.XSize
        self.YSize = self.band.YSize
        self.elevation = self.band.ReadAsArray(0, 0, self.XSize, self.YSize )

        # Fetch the coefficients for transforming between 
        # pixel/line (P,L) raster space => projection coordinates (Xp,Yp) space
        # Xp = GT[0] + P*GT[1] + L*GT[2]
        # Yp = GT[3] + P*GT[4] + L*GT[5]

        self.GT = self.ds.GetGeoTransform()
        self.long = self.GT[0] + np.arange(self.XSize) * self.GT[1]
        self.lat = self.GT[3] + np.arange(self.YSize) * self.GT[5]
        self.meshgrid_long, self.meshgrid_lat = np.meshgrid(self.long, self.lat)

        print("build meshgrids")
        self.X, self.Y = epsg.wgs84LongLatToEpsg((self.meshgrid_long, self.meshgrid_lat), epsg.epsg3948)
        print("build rectBivariateSpline")
        self.rectBivariateSpline = interp.RectBivariateSpline(self.lat[::-1], self.long, self.elevation[::-1,:] )

    def __getDataset(self):
        ds = gdal.Open(self.filename, gdal.GA_ReadOnly)
        return ds

    def info(self):
        print("Driver: {}/{}".format(self.ds.GetDriver().ShortName,
        self.ds.GetDriver().LongName))
        print("Size is {} x {} x {}".format(self.ds.RasterXSize,
        self.ds.RasterYSize, self.ds.RasterCount))
        print("Projection is {}".format(self.ds.GetProjection()))

        print("Origin = ({}, {})".format( self.GT[0], self.GT[3] ) )
        print("Pixel Size = ({}, {})".format( self.GT[1], self.GT[5] ) )
        print("Band Type={}".format(gdal.GetDataTypeName(self.band.DataType)))

    def plot(self):
        vmin = 0
        vmax = 300
        left = self.GT[0]
        right = self.GT[0] + self.XSize * self.GT[1]
        bottom = self.GT[3] + self.YSize * self.GT[5]
        top = self.GT[3]

        plt.figure()
        plt.imshow(self.elevation, extent=(left, right, bottom, top), vmin=vmin, vmax=vmax, cmap='terrain')
        plt.colorbar()
        plt.grid()

    def plotXY(self, sceneX=None, sceneY=None):
        fig, ax = plt.subplots(1,1)
        vmin = 0
        vmax = 130
        im = ax.pcolormesh(self.X, self.Y, self.elevation, cmap='terrain', vmin=vmin, vmax=vmax)
        if sceneX is not None and sceneY is not None:
            ax.plot( sceneX[::10,::10], sceneY[::10,::10], ".k", label = "scene" )
        ax.set_aspect("equal")
        ax.grid()
        plt.colorbar(im)
