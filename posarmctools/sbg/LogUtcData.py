from .Log import Log
import numpy as np

import logging
logger = logging.getLogger(__name__)

class LogUtcData(Log):
	"""docstring for LogUtcData"""
	def __init__(self, name, logLoader, delimiter=None):
		super().__init__(name, logLoader, delimiter)
		self.addSecOfDay()
		self.printLoopbacksTimes(print)

	def toStr(self, idx):
		h = self.names.index('hour')
		m = self.names.index('minute')
		s = self.names.index('second')
		ns = self.names.index('nanoSecond')
		utc = self.data[idx, :]
		return f"{utc[h]:.0f}:{utc[m]:.0f}:{utc[s] + utc[ns] * 1e-9:.3f}"

	def addSecOfDay(self):
		idxHour = self.names.index("hour")
		idxMin = self.names.index("minute")
		idxSec = self.names.index("second")
		idxNSec = self.names.index("nanoSecond")
		utcSeconds = self.data[:, idxHour] * 3600 + self.data[:, idxMin] * 60 +\
			self.data[:, idxSec] + self.data[:, idxNSec] * 1e-9
		self.secOfDay = utcSeconds

	def getYearMonthDay(self, line=1):
		# sometimes, the first line, index 0, is not valid
		return [self.data[line, idx]\
			for idx in [self.names.index(col) for col in ("year", "month", "day")]]

	def printLoopbacksTimes(self, printer):
		msg = "start ... " + self.toStr(self.loopbacks[0])
		printer(msg)
		logger.info(msg)
		if len(self.loopbacks) > 2:
			for loop, loopback in enumerate(self.loopbacks[1:-1]):
				msg = f"loop {loop} ... " + self.toStr(loopback)
				printer(msg)
				logger.info(msg)
		else:
			msg = "no loopback detected"
			printer(msg)
			logger.info(msg)
		msg = "stop ... " + self.toStr(self.loopbacks[-1] - 1)
		printer(msg)
		logger.info(msg)
