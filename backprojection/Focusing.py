import logging, os

import numpy as np

from .Conf import Conf
from .Track import Track
from .Scene import Scene
from .Signal import Signal
from .Dted import Dted
from .loadbackprojection import MyParameters_LETG

from posarmctools import posar
import posarmctools.epsgtools as epsg

logger = logging.getLogger(__name__)

def getPhi_a(res, lambda_c):
    phi_a = 2 * np.arcsin(lambda_c / (4 * res)) * 180 / np.pi
    return phi_a

class Focusing(object):
    def __init__(self, filename, tag="Sa", alongTrack=1, pointsFile="py", corr=False):
        self.filename = filename
        self.pointsFile = pointsFile
        self.tag = tag
        self.alongTrack = alongTrack
        self.corr = corr

        print(f"======= Conf")
        self.conf = Conf(filename, pointsFile=pointsFile)

        print("======= Track")
        self.track = Track(self.conf, filename=None)

        print("======= Load analytic signal...")
        self.signal = Signal(self.conf, tag=tag, filename=None)
        if self.signal.load() is False: # load signal.sa and signal.coupling
            print("/!\\ signal.load() failed, try to build the analytic signal /!\\")
            record = posar.Record(self.conf.rec_dir, "record", "bin", version="X_v1")
            A_reshaped = record.readBins("ADLINKCh0")
            self.signal.build(A_reshaped)
            self.signal.save()
        else:
            print("existing analytic signal found and loaded")
        self.sa = self.signal.sa
        self.avg = self.signal.avg

        print("======= Load antenna positions...")
        self.loadAntennaPositions()

        print("======= Load Digital Terrain Elevation Data...")
        if self.conf.hScene == -1:
            self.dted = Dted(self.conf.dted)
        else:
            print(f"DTED not loaded, conf.hScene {self.conf.hScene:.2f}")

        print("======= Scene")
        self.buildScene()

        print("======= Set focusingParameters")
        self.parameters = MyParameters_LETG()
        self.setFocusingParameters()

        self.out_dir = os.path.join(self.conf.root_dir, "FOCUSING", self.conf.data_date)
        try:
            os.makedirs(self.out_dir, exist_ok=True) 
            logger.info(f"{self.out_dir} created successfully")
        except OSError as error:
            logger.error(f"{self.out_dir} can not be created, error: {error}")

    def setFocusingParameters(self):
        Nover = self.conf.overSamplingRatio * self.signal.Nf
        rangeResolution = self.conf.c / (2 * self.conf.B0)
        self.r_over = np.arange(Nover) * rangeResolution / self.conf.overSamplingRatio
        dr_over = self.r_over[1] - self.r_over[0]
        
        self.parameters.Nx = self.scene.nX
        self.parameters.Ny = self.scene.nY
        self.parameters.Nover = self.r_over.size
        self.parameters.dx = dr_over
        self.parameters.Naz = self.signal.Naz
        self.parameters.Nf = self.signal.Nf
        self.parameters.hScene = self.conf.hScene

        self.parameters.phi_a_deg = self.conf.phi_a_deg
        self.parameters.uxx = self.track.ux[0]
        self.parameters.uxy = self.track.ux[1]
        self.parameters.meanX = self.scene.X_mean
        self.parameters.meanY = self.scene.Y_mean
        self.parameters.kc = self.conf.params.kc

        self.azimuthResolution = self.conf.lambda_c / (4 * np.sin( self.conf.phi_a_deg * np.pi / 180 / 2 ) )

        print(f"Nf {self.signal.Nf}, Naz {self.signal.Naz}, phi_a_deg {self.conf.phi_a_deg}, azRes {self.azimuthResolution:.3f}")
        print("range from {:.2f}m to {:.2f}m, resolution = {}m, oversampled = {}m, ".format(
            self.r_over[0], self.r_over[-1], rangeResolution, rangeResolution / self.conf.overSamplingRatio))

    def loadAntennaPositions(self):
        # load positions for all ramps
        if self.corr:
            tag = "_corr"
        else:
            tag = ""
        if self.conf.nav:
            filename = os.path.join(self.conf.out_dir, f"n_time_xyz_nav_a{tag}.npy")
        else:
            filename = os.path.join(self.conf.out_dir, f"n_time_xyz_gps_a{tag}.npy")
        
        print(f"load {filename}")
        self.xyz = np.load(filename)
        self.xa = self.xyz[:,2]
        self.ya = self.xyz[:,3]
        self.za = self.xyz[:,4]
        self.xa_mean = np.mean(self.xa)
        self.ya_mean = np.mean(self.ya)
        self.za_mean = np.mean(self.za)
        print(f"xa.shape {self.xa.shape}")
        print(f"xa_mean = {self.xa_mean:.2f}, ya_mean = {self.ya_mean:.2f}, za_mean = {self.za_mean:.2f}")

    def addPtsEpsg(self, ax, refs, color='r', marker='x'):
        for ref in refs:
            ax.plot(self.conf.ptsEpsg[ref][0], self.conf.ptsEpsg[ref][1],
            marker=marker, color=color)

    def buildScene(self):
        if self.conf.groundRange == 1:
            if self.alongTrack:
                self.scene = Scene(self.conf, self.conf.proj, track=self.track, target_ll=self.conf.center)
            else:
                self.scene = Scene(self.conf, self.conf.proj, target_ll=self.conf.center)
            if self.conf.hScene == -1:
                self.sceneZ = self.dted.rectBivariateSpline.ev(self.scene.lat, self.scene.long)
            else:
                self.sceneZ = np.ones(self.scene.X.shape) * self.conf.hScene
                print(f"hScene {self.conf.hScene:.2f}")
        else:
            self.scene = Scene(self.conf, self.conf.proj,\
                track=self.track, target_ll=self.conf.center, xyz=self.xyz)
            self.sceneZ = np.ones(self.scene.X.shape) * self.scene.hFocusing

    def rebuildScene(self):
        print(f"======= Conf")
        self.conf = Conf(self.filename, pointsFile=self.pointsFile)

        print("======= Scene")
        self.buildScene()

        print("======= Set focusingParameters")
        self.setFocusingParameters()