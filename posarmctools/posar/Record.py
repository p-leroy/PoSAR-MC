import re, os
import numpy as np
from .RecordBin import RecordBin
from .Timestamps import Timestamps
from .PosarMCParameters import PosarMCParameters_v2, PosarXParameters

import logging
logger = logging.getLogger(__name__)

class Record(object):
	"""docstring for Record"""
	def __init__(self, rec_dir, prefix="record", extension="bin", utc=None, period=1e6, version="X_v1"):
		print(f"rec_dir {rec_dir}")
		self.rec_dir = rec_dir
		self.day_hour = rec_dir.split('/')[-1]
		logger.info(f"\n\n - RECORD - {self.day_hour}\n")
		y, m, d, h, mn, s = self.day_hour.split('_')
		self.y = int(y)
		self.m = int(m)
		self.d = int(d)
		self.h = int(h)
		self.mn = int(mn)
		self.s = int(s)
		self.hour = self.day_hour[-8:]

		self.out_dir = os.path.join('/', *self.rec_dir.split('/')[:-1], 'OUT', self.day_hour)
		try:
			os.makedirs(self.out_dir, exist_ok = True) 
			logger.info(f"{self.out_dir} created successfully")
		except OSError as error:
			logger.error(f"{self.out_dir} can not be created, error: {error}") 
		self.prefix = prefix
		self.extension = extension
		self.utc = utc

		self.binNums, self.idx = self.getInfo(rec_dir, prefix, extension)

		self.recBins = self.getRecBins(self.utc)
		self.nBins = len(self.recBins)
		self.modificationTimestamps = [recBin.modificationTimestamp for recBin in self.recBins]
		self.modificationSeconds = [recBin.seconds for recBin in self.recBins]
		self.firstRec = self.recBins[0]
		self.firstRec.printTimes(printer=logger.info)
		self.delta, self.UTCShift = self.getDelta()
		if abs(self.delta) > 10:
			logger.error(f"large delta between name and modification time = {self.delta:.3f} s")
		else:
			logger.info(f"delta between name and modification time = {self.delta:.3f} s")

		self.timestamps = Timestamps(rec_dir, utc)
		self.bufferNumber, self.logEvents, self.idxBuffer, self.idxTime =\
			self.timestamps.checkTimestamps(period=period)
		self.shiftedLogEvents = self.__shiftLogEvents()

		self.parameters = self.__getParameters(version)

		self.extendedLogEvents, self.rampNumber_allRamps, self.timestamps_allRamps = self.__getInterpolatedTimestamps()

	def getInfo(self, rec_dir, prefix, extension):
		ret = -1
		num = []
		p = re.compile("^" + prefix + "\d+\." + extension + "$")
		logger.info(f"rec_dir   *** {rec_dir}")
		if os.path.isdir(rec_dir):
			for root, dirs, files in os.walk(rec_dir):
				for file in files:
					if p.match(file) != None:
						num.append(int(file.split(prefix)[1].split("." + extension)[0]))
			if num != []:
				num.sort()
				logger.info(f"{len(num)} records found in directory, " 
					+ f"start = {prefix}{num[0]}.{extension}, stop = {prefix}{num[-1]}.{extension}")
				diff = num[1] - num[0]
				idx = np.where(np.diff(num) != diff)
				if idx[0].shape[0] == 0:
					logger.info("records are contiguous")
				else:
					logger.error("records are not contiguous")
				ret = 0
			else:
				logger.error(f"no file {prefix}.{extension} found in {rec_dir}")
				ret = -1
		else:
			logger.error("rec_dir does not exist")
			ret = -1

		if ret == -1:
			return -1, -1
		else:
			return num, idx

	def print(self):
		print(f"{self.rec_dir}")
		print(f"looking for files {self.prefix}xxx.{self.extension}")
		print(f"{len(self.binNums)} files found")
		print(f"numbers are {self.binNums}")

	def getDelta(self):
		modificationSeconds = self.firstRec.modification.tm_hour * 3600 + self.firstRec.modification.tm_min * 60 +\
			self.firstRec.modification.tm_sec + self.firstRec.mFloat
		nameSeconds = self.h * 3600 + self.mn * 60 + self.s
		delta = nameSeconds - modificationSeconds
		UTCShift = round( delta / 3600 )
		logger.info(f"delta = {delta:.3f}, estimated UTCShift = {UTCShift}")
		nameSeconds = nameSeconds - 3600 * UTCShift
		delta = nameSeconds - modificationSeconds
		return delta, UTCShift

	def getRecBins(self, utc=None):
		recBins = [RecordBin(f"{self.rec_dir}/{self.prefix}{binNum}.{self.extension}", utc) for binNum in self.binNums]
		return recBins

	def __shiftLogEvents(self):
		aux = np.array(self.logEvents)
		loopback = self.recBins[0].loopback
		aux += loopback * 2**32
		innerLoopbacks = np.where(np.diff(aux) < 0)[0] + 1
		if innerLoopbacks.size != 0:
			print(f"inner loopback detected at {innerLoopbacks}")
			for index in innerLoopbacks:
				aux[index:] += 2**32
		return aux

	def __getParameters(self, posarVersion):
		params = None
		filename = os.path.join(self.rec_dir, self.day_hour + "_parameters.xml")
		if posarVersion == "X_v1":
			params = PosarXParameters(filename)
		elif posarVersion == "C_v2":
			params = PosarMCParameters_v2(filename)
		return params

	def __getInterpolatedTimestamps(self):
		nEventsRequested = self.nBins + 1 # requested for the interpolation
		rampNumber_firstRampInEachFile = np.arange(nEventsRequested) * self.parameters.rampsPerFile
		rampNumber_allRamps = np.arange(self.nBins * self.parameters.rampsPerFile)
		
		nEventsToAdd = nEventsRequested - self.shiftedLogEvents.size
		if nEventsToAdd > 0:
			logger.info(f"{self.day_hour} number of added events {nEventsToAdd}")
			lastDelta = self.shiftedLogEvents[-2] + self.shiftedLogEvents[-1]
			eventsToAdd = self.shiftedLogEvents[-1] + np.array(range(nEventsToAdd)) * lastDelta
			extendedLogEvents = np.r_[self.shiftedLogEvents, eventsToAdd]
		else:
			extendedLogEvents = self.shiftedLogEvents[:nEventsRequested]
			
		timestamps_allRamps = np.interp(rampNumber_allRamps, rampNumber_firstRampInEachFile, extendedLogEvents)
		
		return extendedLogEvents, rampNumber_allRamps, timestamps_allRamps

	def readBins(self, version, firstBin=None, lastBin=None):
		if firstBin == None and lastBin == None:
			binsToRead = self.recBins
			firstBin = self.binNums[0]
			lastBin = self.binNums[-1]
		else:
			idxFirst = self.binNums.index(firstBin)
			idxLast = self.binNums.index(lastBin)
			binsToRead = [self.recBins[index] for index in range(idxFirst, idxLast+1)]
		nbFiles = len(binsToRead)

		samplesPerFile = self.parameters.samplesPerRamp * self.parameters.rampsPerFile
		adc_A = np.zeros((nbFiles, samplesPerFile))
		print(f"reading binaries from {firstBin} to {lastBin}")
		print(f"number of files to read {nbFiles}, samplesPerFile {samplesPerFile}")
		print(f"samplesPerRamp {self.parameters.samplesPerRamp}, rampsPerFile {self.parameters.rampsPerFile}")
		for index, rec in enumerate(binsToRead):
			rec.read(samplesPerFile, adc_A[index, :], version=version)
		return adc_A.reshape(nbFiles * self.parameters.rampsPerFile, self.parameters.samplesPerRamp)
