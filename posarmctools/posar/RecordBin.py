import numpy as np

import logging
logger = logging.getLogger(__name__)

from .Synchro import Synchro

class RecordBin(Synchro):
	"""docstring for RecordBin"""
	def __init__(self, name, utc=None):
		super().__init__(name, utc)

	def read(self, samplesPerFile, timeSerie_A, timeSerie_B=None, version="ADLINKCh0"):
		print(f"reading {self.name}")
		with open(self.name, 'rb') as fd:
			if version == "ADLINKCh0":
				scalingFactor = 2 / 65535
				offset = -32768
				dum = np.fromfile(fd, dtype = np.uint16)
				timeSerie_A[:] = (dum[:] + offset) * scalingFactor
			elif version == "ADLINK":
				dum = np.fromfile(fd, dtype = np.uint16)
				timeSerie_A[:] = dum[ 0 : 2 * samplesPerFile : 2 ]
				timeSerie_B[:] = dum[ 1 : 2 * samplesPerFile : 2 ]
			elif version == "old":
				dum = np.fromfile(fd, dtype = np.int16)
				timeSerie_A[:] = dum[ 0 : 2 * samplesPerFile : 2 ]
				timeSerie_B[:] = dum[ 1 : 2 * samplesPerFile : 2 ]
			else:
				logger.error("readFile: unknown file version [possible options: ADLINKCh0, ADLINK, old]")
