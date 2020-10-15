from .Log import Log
import numpy as np

import logging
logger = logging.getLogger(__name__)

class LogEkfEuler(Log):
    """docstring for LogEkfEuler"""
    def __init__(self, name, logLoader, delimiter=None):
        super().__init__(name, logLoader, delimiter)

        self.roll = self.data[:, self.names.index('roll')]
        self.pitch = self.data[:, self.names.index('pitch')]
        self.yaw = self.data[:, self.names.index('yaw')]

    def interpolate(self, record):
	    r = np.interp(record.timestamps_allRamps, self.timestamps, self.roll)
	    p = np.interp(record.timestamps_allRamps, self.timestamps, self.pitch)
	    y = np.interp(record.timestamps_allRamps, self.timestamps, self.yaw)
	    rpy_a = np.stack(
            (record.rampNumber_allRamps, record.timestamps_allRamps, r, p, y), -1)

	    timestamps_allRamps_b = record.timestamps_allRamps + record.parameters.rampPeriod / 2

	    r_b  = np.interp(timestamps_allRamps_b, self.timestamps, self.roll)
	    p_b = np.interp(timestamps_allRamps_b, self.timestamps, self.pitch)
	    y_b  = np.interp(timestamps_allRamps_b, self.timestamps, self.yaw)
	    rpy_b = np.stack(
		    (record.rampNumber_allRamps, timestamps_allRamps_b, r_b, p_b, y_b), -1)

	    return rpy_a, rpy_b
