import logging, os

import numpy as np

from posarmctools import epsgtools as epsg

from .LogEkfEuler import LogEkfEuler
from .LogEkfNav import LogEkfNav
from .LogEventB import LogEventB
from .LogGpsPos import LogGpsPos
from .LogLoader import LogLoader
from .LogUtcData import LogUtcData
from .LogGpsVel import LogGpsVel

logger = logging.getLogger(__name__)

class Logs(object):
    """docstring for Logs"""
    def __init__(self, conf, epsg3xxx, delimiter=None):
        self.conf = conf
        self.epsg3xxx = epsg3xxx

        self.loader = LogLoader(conf.dir_sbg, conf.session, conf.day, conf.sbg_hours)

        print("load gps")
        logger.info("load gps")
        self.gps = LogGpsPos("sbgLogGpsPos.csv", self.loader, delimiter=",")
        self.x, self.y = epsg.wgs84LongLatToEpsg((self.gps.long, self.gps.lat), epsg3xxx)
        self.gpsEpsg = epsg.wgs84LongLatToEpsg((self.gps.long, self.gps.lat), epsg3xxx)
        self.vel = LogGpsVel("sbgLogGpsVel.dat", self.loader)

        np.save(conf.out_dir + "/gps_epsg_transform", self.gpsEpsg)

        print("load utc")
        logger.info("load utc")
        self.utc = LogUtcData("sbgLogUtcData.dat", self.loader)

        print("load eventB")
        logger.info("load eventB")
        self.eventB = LogEventB("sbgLogEventB.dat", self.loader, utc=self.utc)

        print("load ekf")
        logger.info("load ekf")
        self.euler = LogEkfEuler("sbgLogEkfEuler.dat", self.loader)
        self.nav = LogEkfNav("sbgLogEkfNav.dat", self.loader)
        self.navEpsg = epsg.wgs84LongLatToEpsg((self.nav.long, self.nav.lat), epsg3xxx)

        np.save(os.path.join(self.conf.out_dir, "nav_epsg_transform"), self.navEpsg)
