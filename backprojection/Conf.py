import configparser, json, logging, os, sys

import numpy as np
import matplotlib.pyplot as plt

from posarmctools import posar
from posarmctools import epsgtools as epsg

logger = logging.getLogger(__name__)

class Conf(object):
    def __init__(self, filename, pointsFile="py"):
        self.filename = filename
        self.name = filename.split("/")[-1].split(".ini")[0]
        self.pointsFile = pointsFile

        # [posar]
        self.dir_posar, self.day, self.hour = self.__get_posar()
        self.data_date = f"{self.day}_{self.hour}"
        self.root_dir = os.path.join(self.dir_posar, self.day)
        self.rec_dir = os.path.join(self.root_dir, self.data_date)
        self.out_dir = os.path.join(self.root_dir, "OUT", self.data_date)
        self.record = posar.Record(self.rec_dir, "record", "bin", version="X_v1")
        # insert dir_posar in the path to have an access to the references
        sys.path.insert(0, self.dir_posar)

        # [scene]
        self.__get_scene()

        # [sbg]
        self.nav = self.__get_sbg()

        # [dted]
        self.__get_dted()

        # [focusing]
        self.__get_focusing()

        self.params_filename = os.path.join(self.rec_dir, self.data_date + "_parameters.xml")
        self.params = posar.PosarXParameters(self.params_filename)
        self.Tp = self.params.rampPeriod / 1e6
        self.B0 = ( self.params.stopFrequency - self.params.startFrequency ) * 1e6
        self.fc = ( self.params.stopFrequency + self.params.startFrequency ) * 1e6 / 2
        self.fs = self.params.samplingFrequency
        self.c = 3e8
        self.lambda_c = self.c / self.fc

    def __get_posar(self):
        with open(self.filename) as f:
            config = configparser.ConfigParser()
            config.read_file(f)
            # [posar]
            dir_posar = config.get("posar", "dir_posar")
            day = config.get("posar", "day")
            hour = config.get("posar", "hour")
            logger.info(f"dir_posar {dir_posar}")
            logger.info(f"day {day}, hour {hour}")
        return dir_posar, day, hour

    def __get_sbg(self):
        with open(self.filename) as f:
            config = configparser.ConfigParser()
            config.read_file(f)
            # [sbg]
            nav = int(config.get("sbg", "nav"))
            logger.info(f"nav {nav}")
        return nav

    def __get_scene(self):
        with open(self.filename) as f:
            config = configparser.ConfigParser()
            config.read_file(f)
            # [scene]
            self.hScene = float(config.get("scene", "hScene"))
            self.d_x = float(config.get("scene", "d_x"))
            self.d_y = float(config.get("scene", "d_y"))
            self.xMin = float(config.get("scene", "xMin"))
            self.xMax = float(config.get("scene", "xMax"))
            self.yMin = float(config.get("scene", "yMin"))
            self.yMax = float(config.get("scene", "yMax"))
            self.center = eval(config.get("scene", "center"))
            self.shiftY = float(config.get("scene", "shiftY"))
            lines = list(filter(None, config.get("scene", "points").splitlines()))
            points = []
            for line in lines:
                points.extend(line.split(" "))
            logger.info(f"points {points}")
            refs = config.get("scene", "reference_points")
            self.projStr = config.get("scene", "projection")
            if self.projStr == "epsg3948":
                self.proj = epsg.epsg3948
                logger.info(f"projection {self.projStr}")
            else:
                self.proj = None
                logger.error(f"unknown projection {self.projStr}")
        if self.pointsFile == "ini":
            refs = os.path.join(self.dir_posar, refs)
        pts, ptsEpsg, ptsDict = epsg.getReferencePoints(refs, self.proj, file=self.pointsFile)
        self.ptsEpsg = {key: ptsEpsg[key] for key in points}

    def __get_dted(self):
        with open(self.filename) as f:
            config = configparser.ConfigParser()
            config.read_file(f)
            # [scene]
            self.dted = config.get("dted", "filename")

    def __get_focusing(self):
        with open(self.filename) as f:
            config = configparser.ConfigParser()
            config.read_file(f)
            # [focusing]
            self.groundRange = int(config.get("focusing", "groundRange"))
            self.overSamplingRatio = int(config.get("focusing", "overSamplingRatio"))
            self.phi_a_deg = float(config.get("focusing", "phi_a_deg"))
            self.upOrDown = config.get("focusing", "upOrDown")
            self.rampUpFirst = int(config.get("focusing", "rampUpFirst"))
            self.window = config.get("focusing", "window")
            self.firstFile = int(config.get("focusing", "firstFile"))
            self.nBins = int(config.get("focusing", "nBins"))
            if self.nBins == -1:
                self.nBins = self.record.nBins
                logger.info(f"firstFile {self.firstFile}, nBins {self.nBins} [auto determination]")
            else:
                logger.info(f"firstFile {self.firstFile}, nBins {self.nBins} [ini file]")

            if self.groundRange == 1:
                self.groundOrSlant = "ground range"
            elif self.groundRange == 0:
                self.groundOrSlant = "slant range"
            else:
                self.groundOrSlant = "unknown"
