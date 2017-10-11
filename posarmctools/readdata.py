import numpy as np

def readFile( filename, samplesPerFile, timeSerie_A, timeSerie_B ):
    fd = open(filename,'rb')
    
    dum = np.fromfile(fd, dtype = np.int16)

    timeSerie_A[:] = dum[ 0 : 2 * samplesPerFile : 2 ]
    timeSerie_B[:] = dum[ 1 : 2 * samplesPerFile : 2 ]
    
    fd.close()
