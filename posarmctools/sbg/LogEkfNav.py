from .Log import Log
import numpy as np

import logging
logger = logging.getLogger(__name__)

class LogEkfNav(Log):
    """docstring for LogEkfNav"""
    def __init__(self, name, logLoader, delimiter=None):
        super().__init__(name, logLoader, delimiter)

        self.lat = self.data[:, self.names.index('Lat')]
        self.long = self.data[:, self.names.index('Long')]
        self.alt = self.data[:, self.names.index('Alt')]
        self.undulation = self.data[:, self.names.index('undulation')]
        self.statusfield = self.data[:, self.names.index('status')].astype("uint32")

        self.solution_mode = np.zeros(self.statusfield.shape)
        for idx, val in enumerate(self.statusfield):
            self.solution_mode[idx] = val & 0b1111
