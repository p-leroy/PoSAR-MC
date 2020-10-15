import numpy as np

class LogSbg(object):
	"""docstring for LogSbg"""
	def __init__(self, data_dir, session, day, hours):
		self.data_dir = data_dir
		self.session = session
		self.day = day
		self.hours = hours
		
	def load(self, log, delimiter=None):
		prefix = f"{self.data_dir}/{self.session}/{self.day}/"
		logs = [f"{prefix}{h}/{log}" for h in self.hours ]
		data = [np.loadtxt(log, skiprows=1, delimiter=delimiter)
		for log in logs]
		# read columns titles
		with open(logs[0], "r") as file:
			names = [string.strip() for string in file.readline().split(delimiter)]
		return data, names
