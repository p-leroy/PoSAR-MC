from .Log import Log
import numpy as np

import logging
logger = logging.getLogger(__name__)

class LogGpsVel(Log):
    """docstring for LogGpsVel"""
    def __init__(self, name, logLoader, delimiter=None):
        super().__init__(name, logLoader, delimiter)

        self.north = self.data[:, self.names.index('velocityNorth')]
        self.east = self.data[:, self.names.index('velocityEast')]
        self.down = self.data[:, self.names.index('velocityDown')]
        self.course = self.data[:, self.names.index('course')]
        self.statusfield = self.data[:, self.names.index('status')].astype("uint32")
        self.vel = ( self.north**2 + self.east**2 + self.down**2) ** 0.5
        
        # STATUS
        # 0 SBG_ECOM_VEL_SOL_COMPUTED     A valid solution has been computed.
        # 1 SBG_ECOM_VEL_INSUFFICIENT_OBS Not enough valid SV to compute a solution.
        # 2 SBG_ECOM_VEL_INTERNAL_ERROR   An internal error has occurred.
        # 3 SBG_ECOM_VEL_LIMIT            Velocity limit exceeded.
        # TYPE
        # 0 SBG_ECOM_VEL_NO_SOLUTION  No valid velocity solution available.
        # 1 SBG_ECOM_VEL_UNKNOWN_TYPE An unknown solution type has been computed.
        # 2 SBG_ECOM_VEL_DOPPLER      A Doppler velocity has been computed.
        # 3 SBG_ECOM_VEL_DIFFERENTIAL A velocity has been computed between two positions.
        self.status = np.zeros(self.statusfield.shape)
        self.type = np.zeros(self.statusfield.shape)
        for idx, val in enumerate(self.statusfield.shape):
            self.status[idx] = val & 0b000011
            self.type[idx] = np.right_shift( (val & 0b111111000000), 6)
