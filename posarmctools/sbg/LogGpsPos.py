from .Log import Log
import numpy as np

import logging
logger = logging.getLogger(__name__)

class LogGpsPos(Log):
	"""docstring for LogGpsPos"""
	def __init__(self, name, logLoader, delimiter=None):
		super().__init__(name, logLoader, delimiter)
		
		self.statusfield = self.data[:, self.names.index('status')].astype("uint32")
		self.lat = self.data[:, self.names.index('latitude')]
		self.long = self.data[:, self.names.index('longitude')]
		self.alt = self.data[:, self.names.index('altitude')]
		self.undulation = self.data[:, self.names.index('undulation')]
		self.status = np.zeros(self.statusfield.shape)
		self.type = np.zeros(self.statusfield.shape)
		for idx, val in enumerate(self.statusfield):
			self.status[idx] = val & 0b000011
			self.type[idx] = np.right_shift( (val & 0b111111000000), 6)
