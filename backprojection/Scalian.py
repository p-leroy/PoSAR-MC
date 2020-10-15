import configparser, os

import numpy as np
import matplotlib.pyplot as plt

from .Conf import Conf
from .Signal import Signal
from .Track import Track

from posarmctools import posar
import posarmctools.epsgtools as epsg

class Scalian(object):
    def __init__(self, filename, pointsFile="ini", mode="default"):
        self.filename = filename

        print("======= Conf")
        self.conf = Conf(filename, pointsFile=pointsFile)
        
        self.out_dir = os.path.join(self.conf.root_dir, "SCALIAN", self.conf.data_date)
        try:
            os.makedirs(self.out_dir, exist_ok=True)
            print(f"{self.out_dir} created successfully (or already existing)")
        except OSError as error:
            print(f"{self.out_dir} can not be created, error: {error}")

        if mode == "default":
            self.dataShape, self.dataDtype = self.load_build_save("Up")
            self.trackShape, self.trackDtype = self.load_save_nav_euler("Up")
            self.write_ini()
            
            self.dataShape, self.dataDtype = self.load_build_save("Down")
            self.trackShape, self.trackDtype = self.load_save_nav_euler("Down")
            self.write_ini()
        elif mode is None:
            print(f"mode is {None}, nothing read, neither posar data nor navigation data")

    def load_data(self, upOrDown):
        # set window to "None"
        self.conf.upOrDown = upOrDown
        self.conf.window = "None"
        print(f"======= Load analytic signal, upOrDown {upOrDown}, window {self.conf.window}...")

        signal = Signal(self.conf, filename=None)
        if signal.load() is False: # load signal.sa and signal.coupling
            print("/!\\ signal.load() failed, try to build the analytic signal /!\\")
            record = posar.Record(self.conf.rec_dir, "record", "bin", version="X_v1")
            A_reshaped = record.readBins("ADLINKCh0")
            signal.build(A_reshaped)
            signal.save()
        else:
            print("existing analytic signal found and loaded")

        return signal

    def load_build_save(self, upOrDown):
        signal = self.load_data(upOrDown)
        filename = f'{self.out_dir}/Sa_{self.conf.data_date}_{self.conf.upOrDown}.data'
        print(f"save {filename}")
        signal.sa.astype('complex64').tofile(filename)
        return signal.sa.shape, 'complex64'

    def loadAntennaPositionsAndAngles(self):
        if self.conf.upOrDown == "Up":
            if self.conf.rampUpFirst == 1:
                ab = "a"
            elif self.conf.rampUpFirst == 0:
                ab = "b"
            else:
                ab = "a"
                print(f"unknown value for rampUpFirst {self.conf.rampUpFirst}, ab set to {ab}")
        elif self.conf.upOrDown == "Down":
            if self.conf.rampUpFirst == 1:
                ab = "b"
            elif self.conf.rampUpFirst == 0:
                ab = "a"
            else:
                ab = "a"
                print(f"unknown value for rampUpFirst {self.conf.rampUpFirst}, ab set to {ab}")
        else:
            print(f"unknown value for upOrDown {self.conf.upOrDown}")
        # load positions for all ramps
        if self.conf.nav:
            filename = os.path.join(self.conf.out_dir, f"n_time_lla_nav_{ab}.npy")
        else:
            filename = os.path.join(self.conf.out_dir, f"n_time_lla_gps_{ab}.npy")
        print(f"load {filename}")
        lla = np.load(filename)
        lat = lla[:,2]
        lon = lla[:,3]
        alt = lla[:,4]
        lat_mean = np.mean(lat)
        lon_mean = np.mean(lon)
        alt_mean = np.mean(alt)
        print(f"lat.shape {lat.shape}")
        print(f"lat_mean = {lat_mean:.2f}, lon_mean = {lon_mean:.2f}, alt_mean = {alt_mean:.2f}")

        #load roll, pitch, yaw for all ramps
        filename = os.path.join(self.conf.out_dir, f"n_time_rpy_{ab}.npy")
        print(f"load {filename}")
        rpy = np.load(filename)

        return lla[:, 2:], rpy[:, 2:]

    def loadAntennaXYZ(self):
        if self.conf.upOrDown == "Up":
            if self.conf.rampUpFirst == 1:
                ab = "a"
            elif self.conf.rampUpFirst == 0:
                ab = "b"
            else:
                ab = "a"
                print(f"unknown value for rampUpFirst {self.conf.rampUpFirst}, ab set to {ab}")
        elif self.conf.upOrDown == "Down":
            if self.conf.rampUpFirst == 1:
                ab = "b"
            elif self.conf.rampUpFirst == 0:
                ab = "a"
            else:
                ab = "a"
                print(f"unknown value for rampUpFirst {self.conf.rampUpFirst}, ab set to {ab}")
        else:
            print(f"unknown value for upOrDown {self.conf.upOrDown}")
        # load positions for all ramps
        if self.conf.nav:
            filename = os.path.join(self.conf.out_dir, f"n_time_xyz_nav_{ab}.npy")
        else:
            filename = os.path.join(self.conf.out_dir, f"n_time_xyz_gps_{ab}.npy")
        print(f"load {filename}")
        xyz = np.load(filename)
        xa = xyz[:,2]
        ya = xyz[:,3]
        za = xyz[:,4]
        xa_mean = np.mean(xa)
        ya_mean = np.mean(ya)
        za_mean = np.mean(za)
        print(f"xa.shape {xa.shape}")
        print(f"xa_mean = {xa_mean:.2f}, ya_mean = {ya_mean:.2f}, za_mean = {za_mean:.2f}")

        return xyz[:, 2:]

    def load_save_nav_euler(self, upOrDown, dtype='float64'):
        self.conf.upOrDown = upOrDown
        print(f"======= Load antenna positions {upOrDown}...")
        print(f"nav: {self.conf.nav}, upOrDown {upOrDown}, rampUpFirst {self.conf.rampUpFirst}")
        xyz, rpy = self.loadAntennaPositionsAndAngles()
        xyz_rpy = np.c_[xyz, rpy]
        filename = f'{self.out_dir}/Sa_{self.conf.data_date}_{self.conf.upOrDown}.track'
        print(f"save {filename}")
        xyz_rpy.astype(dtype).tofile(filename)
        return xyz.shape, dtype
        
    def write_ini(self, suffix=None):
        config = configparser.ConfigParser()
        config['data'] = {}
        config['data']['shape'] = str(self.dataShape)
        config['data']['type'] = str(self.dataDtype)
        config['track'] = {}
        config['track']['shape'] = str(self.trackShape)
        config['track']['type'] = str(self.trackDtype)
        config['track']['col0'] = "latitude [°]"
        config['track']['col1'] = "longitude [°]"
        config['track']['col2'] = "altitude [°]"
        config['track']['col3'] = "roll [rad]"
        config['track']['col4'] = "pitch [rad]"
        config['track']['col5'] = "yaw [rad]"
        config['parameters'] = {}
        config['parameters']['sampling_frequency'] = str(self.conf.params.samplingFrequency)
        config['parameters']['range_resolution'] = str(self.conf.c / (2 * self.conf.B0))
        config['parameters']['start_frequency'] = str(self.conf.params.startFrequency)
        config['parameters']['stop_frequency'] = str(self.conf.params.stopFrequency)
        config['parameters']['ramp_period'] = str(self.conf.params.rampPeriod)
        
        if suffix is None:
            filename = os.path.join(self.out_dir, f"Sa_{self.conf.data_date}_{self.conf.upOrDown}.ini")
        else:
            filename = os.path.join(self.out_dir, f"Sa_{self.conf.data_date}_{self.conf.upOrDown}_{suffix}.ini")
        with open(filename, 'w') as configfile:
            print(f"write {filename}")
            config.write(configfile)

    def extract_data(self, target_ll, nbPointsToKeep=30000):
        # load antenna positions and angles
        lla, rpy = self.loadAntennaPositionsAndAngles()
        xyz = self.loadAntennaXYZ()
        xa = xyz[:, 0]
        ya = xyz[:, 1]

        print("======= Track")
        # load track
        track = Track(self.conf, filename=None)
        targetLat = target_ll[0]
        targetLong = target_ll[1]
        targetX, targetY = epsg.wgs84LongLatToEpsg((targetLong, targetLat), self.conf.proj)

        # project target on the track
        A = (track.refX, track.refY)
        C = (targetX, targetY)
        print(f"A {A}")
        print(f"C {C}")
        AC = (C[0] - A[0], C[1] - A[1])
        AC_dot_Ux = (AC[0] * track.ux[0] + AC[1] * track.ux[1])
        refX = A[0] + AC_dot_Ux * track.ux[0]
        refY = A[1] + AC_dot_Ux * track.ux[1]

        distanceRefTrackXY = ((refX - xa)**2 + (refY - ya)**2)**0.5
        min_distanceRefTrackXY = np.amin(distanceRefTrackXY)
        idx = np.where(distanceRefTrackXY == min_distanceRefTrackXY)[0][0]
        print(f"min_distanceRefTrackXY {min_distanceRefTrackXY}, idx {idx}")
        
        plt.figure()
        plt.plot(A[0], A[1], 'o', label="A")
        plt.plot(C[0], C[1], 'o', label="C")
        plt.plot(refX, refY, 'o', label="ref")
        plt.plot(xa, ya)
        plt.plot(xa[idx], ya[idx], 'D')
        plt.gca().set_aspect('equal')
        plt.legend()
        plt.grid()

        slice_start = idx - int(nbPointsToKeep / 2)
        slice_stop = idx + int(nbPointsToKeep / 2)
        slice_ = slice(slice_start, slice_stop)
        print(f"{slice_}")

        # save navigation and euler data slice
        lla_rpy = np.c_[lla[slice_, :], rpy[slice_, :]]
        length = lla_rpy.shape[0]
        filename = f'{self.out_dir}/Sa_{self.conf.data_date}_{self.conf.upOrDown}_{length}.track'
        print(f"save {filename}")
        lla_rpy.astype("float64").tofile(filename)

        # load posar data
        signal = self.load_data(self.conf.upOrDown)

        # save posar data slice
        filename = f'{self.out_dir}/Sa_{self.conf.data_date}_{self.conf.upOrDown}_{length}.data'
        print(f"save {filename}")
        signal.sa[slice_, :].astype('complex64').tofile(filename)

        # save ini
        self.dataShape = signal.sa[slice_, :].shape
        self.dataDtype = signal.sa.dtype
        self.trackShape = lla[slice_, :].shape
        self.trackDtype = lla.dtype
        self.write_ini(suffix=length)
