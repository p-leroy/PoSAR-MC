from .Log import Log
import numpy as np

import logging
logger = logging.getLogger(__name__)

class LogEventB(Log):
    """docstring for LogEventB"""
    def __init__(self, name, logLoader, delimiter=None, utc=None):
        super().__init__(name, logLoader, delimiter)
        self.utc = utc
        self.shiftedTimestamps = np.zeros(self.data.shape[0])
        if utc != None:
            self.__shiftTimestamps()

    def __shiftTimestamps(self):
        # check that hours of the utc data are the same as hours of the current LogEventB
        if set(self.utc.loader.hours) != set(self.loader.hours):
            logging.error("utc.loader.hours {utc.loader.hours} != self.loader.hours {self.loader.hours}")
        else:
            for num, hour in enumerate(self.loader.hours):
                sliceLogEvent = self.slices[num]
                sliceUtc = self.utc.slices[num]
                eventsOfThisHour = self.shiftedTimestamps[sliceLogEvent]
                dataOfThisHour = self.data[sliceLogEvent, 0]
                utcTimestampsOfThisHour = self.utc.data[sliceUtc, 0]
                loopback = self.utc.getLoopback(self.utc.slices[num].start)
    
                if self.utc.loopbackDict[hour] != None:
                    # a loopback has occured during this hour
                    idx = np.where( dataOfThisHour >= utcTimestampsOfThisHour[0] )
                    eventsOfThisHour[idx] = dataOfThisHour[idx] + loopback * 2**32
                    idx = np.where( dataOfThisHour < utcTimestampsOfThisHour[-1] )
                    eventsOfThisHour[idx] = dataOfThisHour[idx] + (loopback+1) * 2**32
                else:
                    # no loopback occured during this hour
                    eventsOfThisHour[:] = dataOfThisHour[:] + loopback * 2**32
