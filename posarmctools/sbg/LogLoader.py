import numpy as np

class LogLoader(object):
	"""docstring for LogLoader"""
	def __init__(self, data_dir, session, day, hours):
		self.data_dir = data_dir
		self.session = session
		self.day = day
		self.hours = hours

	def load(self, log, delimiter=None):
		prefix = f"{self.data_dir}/{self.session}/{self.day}/"
		logs = [f"{prefix}{h}/{log}" for h in self.hours]
		data = [np.loadtxt(log, skiprows=1, delimiter=delimiter)
				for log in logs]
		for idx, d in enumerate(data):
			if d.ndim ==1:
				data[idx] = d.reshape(1, -1)
		lengths = [d.shape[0] for d in data]
		slices = [slice(sum(lengths[0:idx]), length + sum(lengths[0:idx]))
					for idx, length in enumerate(lengths)]
		data = np.vstack(data)
		# read columns titles
		with open(logs[0], "r") as file:
			names = [string.strip() for string in file.readline().split(delimiter)]
		return data, names, slices
