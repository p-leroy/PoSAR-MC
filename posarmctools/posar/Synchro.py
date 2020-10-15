import time, os, math
import numpy as np

import logging
logger = logging.getLogger(__name__)

class Synchro(object):
	"""docstring for Synchro"""
	def __init__(self, name, utc=None, ymdhmns=None):
		self.name = name
		self.utc = utc
		self.ymdhmns = ymdhmns

		self.access       = time.gmtime(os.path.getatime(name))
		self.creation     = time.gmtime(os.path.getctime(name))
		self.modification = time.gmtime(os.path.getmtime(name))
		self.mFloat = os.path.getmtime(name) - math.floor(os.path.getmtime(name))
		self.seconds = -1

		m = self.modification
		logger.debug(f"file {name}")
		logger.debug(f"UTC modification {m.tm_year} {m.tm_mon} {m.tm_mday} {m.tm_hour}:{m.tm_min}:{m.tm_sec} + {self.mFloat:.3f}")

		if utc == None:
			self.utcIndex = -1
			self.loopback = -1
			self.modificationTimestamp = -1
		else:
			utcYear, utcMonth, utcDay = utc.getYearMonthDay()
			if utcYear == m.tm_year and utcMonth == m.tm_mon and utcDay == m.tm_mday:
				logger.debug(f"[UTC {utcYear:.0f} {utcMonth:.0f} {utcDay:.0f}] == [modification {m.tm_year} {m.tm_mon} {m.tm_mday}]")
				self.utcIndex = self.getIdxInUtc()
				self.loopback = self.utc.getLoopback(self.utcIndex)
				self.modificationTimestamp = self.utc.data[self.utcIndex, 0]
			else:
				logger.error(f"[UTC {utcYear:.0f} {utcMonth:.0f} {utcDay:.0f}] != [modification {m.tm_year} {m.tm_mon} {m.tm_mday}]")
				if self.ymdhmns != None:
					logger.error(f"ymdhmns used to estimate the modification time: {self.ymdhmns}")
					self.utcIndex = self.getIdxInUtc(self.ymdhmns)
					self.loopback = self.utc.getLoopback(self.utcIndex)
					self.modificationTimestamp = self.utc.data[self.utcIndex, 0]
				else:
					logger.error(f"ymdhmns is None, impossible to estimate the modification time")
					pass

	def printTimes(self, printer=print):
		a = self.access
		c = self.creation
		m = self.modification
		printer(f"{self.name}")
		printer(f"UTC access _____ {a.tm_year} {a.tm_mon} {a.tm_mday} {a.tm_hour}:{a.tm_min}:{a.tm_sec}")
		printer(f"UTC creation ___ {c.tm_year} {c.tm_mon} {c.tm_mday} {c.tm_hour}:{c.tm_min}:{c.tm_sec}")
		printer(f"UTC modification {m.tm_year} {m.tm_mon} {m.tm_mday} {m.tm_hour}:{m.tm_min}:{m.tm_sec} + {self.mFloat:.3f}")

	def printModificationTime(self, printer=print):
		m = self.modification
		printer(f"{self.name} times")
		printer(f"UTC modification {m.tm_year} {m.tm_mon} {m.tm_mday} {m.tm_hour}:{m.tm_min}:{m.tm_sec} + {self.mFloat:.3f}")

	def getIdxInUtc(self, ymdhmns=None):
		if ymdhmns is None:
			self.seconds = self.modification.tm_hour * 3600 + self.modification.tm_min * 60 +\
				self.modification.tm_sec + self.mFloat
		else:
			self.seconds = ymdhmns[3] * 3600 + ymdhmns[4] * 60 + ymdhmns[5]
		before = np.where(self.utc.secOfDay <= self.seconds)[0][-1]
		logger.debug(f"index {before} found in utc data (size {self.utc.timestamps.size}), corresponding to {self.utc.toStr(before)}")
		return before
