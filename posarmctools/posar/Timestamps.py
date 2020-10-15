import numpy as np

import logging
logger = logging.getLogger(__name__)

from .Synchro import Synchro

class Timestamps(Synchro):
    """docstring for Timestamps"""
    def __init__(self, rec_dir, utc=None, suffix="_timeStamps.data"):
        self.day_hour = rec_dir.split('/')[-1]
        self.name = rec_dir + "/" + self.day_hour + suffix
        self.utc = utc
        y, m, d, h, mn, s = self.day_hour.split('_')[0:6]
        self.ymdhmns = (int(y), int(m), int(d), int(h), int(mn), int(s))
        super().__init__(self.name, utc, self.ymdhmns)
    
    def checkTimestamps(self, period=1e6, epsDiff=10, default=4294967295):
        bufferNumber, timestamp = np.loadtxt(self.name, skiprows=1, unpack=True)
        logger.info(f"{self.day_hour} number of timestamps: {len(timestamp)}")
        diffBuffer = bufferNumber[1] - bufferNumber[0]
        idxBuffer = np.where(np.diff(bufferNumber) != diffBuffer)[0]
        if idxBuffer.shape[0] == 0:
            logger.info("buffers are contiguous")
        else:
            if (bufferNumber[idxBuffer+1] == default).all():
                logger.warning(f"{self.day_hour} *** buffers are not contiguous")
                logger.warning(f"{self.day_hour} *** BUT bufferNumber{idxBuffer} == 4294967295")
            else:
                logger.error(f"{self.day_hour} *** buffers are not contiguous")
                logger.error(f"{self.day_hour} *** AND bufferNumber{idxBuffer} != 4294967295")
                
        diffMax = period + epsDiff
        diffMin = period - epsDiff
        diffTimestamp = np.diff(timestamp)
        idxTime = np.where((diffTimestamp > diffMax) | (diffTimestamp < diffMin))[0]
        if idxTime.shape[0] == 0:
            logger.info(f"OK  timestamps are contiguous")
            logger.info(f"min {diffTimestamp.min()} max {diffTimestamp.max()}"\
            + f" mean {diffTimestamp.mean():.3f} std {diffTimestamp.std():.3f}")
        else:
            logger.error(f"{self.day_hour} *** timestamps are not contiguous")
            logger.error(f"{self.day_hour} *** diffTimestamp{idxTime} = {diffTimestamp[idxTime]}")
            
        return bufferNumber, timestamp, idxBuffer, idxTime
