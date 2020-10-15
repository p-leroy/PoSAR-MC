import logging, os

import numpy as np
import scipy.signal as signal

logger = logging.getLogger(__name__)

def real_to_analytic(s):
    Nf  = s.size
    # if Nf is a multiple of 2, fftshift(fftshift(X)) = X
    if ( Nf % 2 ) != 0 :
        logger.error("Nf should be a multiple of 2")
    return np.conjugate(signal.hilbert(s))[::2]

class Signal(object):
    def __init__(self, conf, tag="Sa", filename=None): 
        self.conf = conf
        self.filename = filename
        self.tag = tag
        # analytic signal
        self.sa = None
        self.Naz = None
        self.Nf = None
        self.avg = None

        lastFile = self.conf.firstFile + self.conf.nBins - 1
        firstFile = self.conf.firstFile
        upOrDown = self.conf.upOrDown
        hanning = self.conf.window
        out_dir = self.conf.out_dir
        self.sa_name = f'{out_dir}/{self.tag}_files_{firstFile}_{lastFile}_ramp{upOrDown}_{hanning}.npy'
        self.avg_name = f'{out_dir}/avg_{self.tag}_files_{firstFile}_{lastFile}_ramp{upOrDown}_{hanning}.npy'

    def load(self):
        ret = False

        if self.filename is None:
            if os.path.isfile(self.sa_name):
                # LOAD SA
                logger.info(f"load {self.sa_name}")
                print(f"load {self.sa_name}")
                self.sa = np.load(self.sa_name)
                self.Naz = self.sa.shape[0]
                self.Nf = self.sa.shape[1]
                print(f"sa.shape {self.sa.shape}")
                # LOAD AVG
                logger.info(f"load {self.avg_name}")
                print(f"load {self.avg_name}")
                self.avg = np.load(self.avg_name).reshape(1, -1)
                print(f"avg.shape {self.avg.shape}")
                ret = True
            else:
                print(f"file unfound: {self.sa_name}")
                ret = False

        else:
            print("filename is not None")
            self.sa = None
            self.Naz = None
            self.Nf = None
            ret = False

        return ret

    def load_avg(self):
        ret = False

        if self.filename is None:
            if os.path.exists(self.avg_name) and ret is True:
                logger.info(f"load {self.avg_name}")
                print(f"load {self.avg_name}")
                self.avg = np.load(self.avg_name)
            else:
                print(f"file unfound: {self.avg_name}")
                ret = False

        else:
            print("filename is not None")
            self.avg = None
            ret = False

        return ret

    def build(self, A):
        print(f"shape of the samples matrix = {A.shape}, build analytic signal:")

        Ns_upRamp = int(A.shape[1] / 2)
        Ns = int(A.shape[1] / 4)
        nbRamps = A.shape[0]
    
        sa = np.zeros((nbRamps, Ns), dtype=complex)

        # SHIFT
        shift = 0
        print(f" *** rampUpFirst: {self.conf.rampUpFirst} (1 => the first ramp is a ramp up)")
        if self.conf.upOrDown == "Up" and self.conf.rampUpFirst == 0:
            shift = Ns_upRamp
        if self.conf.upOrDown == "Down" and self.conf.rampUpFirst == 1:
            shift = Ns_upRamp

        # WINDOW
        print(f" *** window: {self.conf.window}")
        if self.conf.window == "hanning":
            win = np.hanning(Ns_upRamp)
        elif self.conf.window == "hamming":
            win = np.hamming(Ns_upRamp)
        elif self.conf.window == "None":
            win = np.ones(Ns_upRamp)
        else:
            logger.error("unknow specified window: {self.confconf.window}")

        # REAL TO ANALYTIC
        print(f" *** ramp: {self.conf.upOrDown}")
        if self.conf.upOrDown == "Up":
            for ramp in range(nbRamps):
                sa[ ramp, :] = real_to_analytic(A[ramp, shift : Ns_upRamp + shift] * win)
        elif self.conf.upOrDown == "Down":
            for ramp in range(nbRamps):
                sa[ramp, :] = real_to_analytic(np.flipud(A[ramp, shift : Ns_upRamp + shift]) * win)
        else:
            logger.error("unknow upOrDown parameter: {self.confconf.upOrDown}")
        print("analytic signal built, perform ifftshift")
            
        self.sa = np.fft.ifftshift(sa, 1)
        self.Naz = self.sa.shape[0]
        self.Nf = self.sa.shape[1]
        self.avg = np.average(sa, 0)

    def save(self):
        firstFile = self.conf.firstFile
        lastFile = self.conf.firstFile + self.conf.nBins - 1
        window = self.conf.window
        upOrDown = self.conf.upOrDown

        print(f"save {self.avg_name}")
        np.save(self.avg_name, self.avg)

        Sa_name = f'{self.conf.out_dir}/Sa_files_{firstFile}_{lastFile}_ramp{upOrDown}_{window}.npy'
        print(f"save {Sa_name}")
        np.save( Sa_name, self.sa)

    def removeAvg(self):
        self.sa = self.sa - self.avg