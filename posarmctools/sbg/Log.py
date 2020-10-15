import numpy as np

class Log(object):
	"""docstring for Log"""
	def __init__(self, name, logLoader, delimiter=None):
		self.name = name
		self.loader = logLoader
		self.data, self.names, self.slices = logLoader.load(name, delimiter=delimiter)
		self.timestamps = np.array(self.data[:, 0])
		self.loopbacks = self.__getLoopbacks()
		self.loopbackDict = self.__getLoopbackDict()

	def __getLoopbacks(self):
		# after a call to this function, timestamps are continuous without loopback
		loopbacks = [0]
		loopbacks.extend(k for k in (np.where(np.diff(self.timestamps) < 0)[0] + 1))
		loopbacks.append(len(self.timestamps))
		for num, loopback in enumerate(loopbacks[1:]):
			range_ = range(loopbacks[num], loopbacks[num + 1])
			self.timestamps[range_] = self.timestamps[range_] + num * 2**32
		return np.array(loopbacks)

	def __getLoopbackDict(self):
		myDict = {hour: None for hour in self.loader.hours}
		for loopback in self.loopbacks[1:-1]:
			for num, slice_ in enumerate(self.slices):
				if loopback >= slice_.start and loopback < slice_.stop:
					myDict[self.loader.hours[num]] = loopback
		return myDict

	def getLoopback(self, idx):
		return np.where(idx >= self.loopbacks)[0][-1]
