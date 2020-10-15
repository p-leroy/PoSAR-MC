import logging

import numpy as np

logger = logging.getLogger(__name__)

def RD_realtocomp( RD ):

    Nf  = RD.size

    # if Nf is a multiple of 2, fftshift(fftshift(X)) = X
    if ( ( Nf % 2 ) != 0 ):
        logger.error( "Nf should be a multiple of 2" )

    h = np.zeros( Nf )
    h[ 0 ] = 1
    h[ int( Nf / 2) ] = 1
    h[ 1 : int(Nf / 2) ] = 2

    RDs = np.fft.ifft( np.fft.fftshift( RD ) )

    RDh = RDs * h

    RDc = np.fft.fftshift( np.fft.fft( RDh ) )

    return RDc[::2]

def build_rd_from_npy( params, npyFile ):
    A = np.load( npyFile )
    print( "data loaded from {}".format(npyFile) )
    print( "shape of the samples matrix = {}".format(A.shape) )
    nbRamps = A.shape[0]
    Ns_upRamp = int(params.samplesPerRamp / 2)
    Ns = int(params.samplesPerRamp / 4)

    RDc = np.zeros( (nbRamps, Ns) , dtype = complex )

    for ramp in range(nbRamps):
        RDc[ ramp, :] = RD_realtocomp( A[ramp, 0 : Ns_upRamp ] ) * 2 / 4096

    return RDc

def build_rd_from_data_rampUp( params, A, rampUpFirst=1, withHanning=0, withHamming=0 ):
    print( "build rd from ramp up" )
    print( "shape of the samples matrix = {}".format(A.shape) )
    nbRamps = A.shape[0]
    Ns_upRamp = int(params.samplesPerRamp / 2)
    Ns = int(params.samplesPerRamp / 4)

    RDc = np.zeros( (nbRamps, Ns) , dtype = complex )

    if rampUpFirst:
        print("ramp up first in the data files")
        shift = 0
    else:
        print("ramp down first in the data files")
        shift = Ns_upRamp

    if withHanning:
        print("with Hanning window")
        window = np.hanning( Ns_upRamp )
    elif withHamming:
        print("with Hamming window")
        window = np.hamming( Ns_upRamp )
    else:
        print("no window")
        window = np.ones( Ns_upRamp )
    
    for ramp in range(nbRamps):
        RDc[ ramp, :] = RD_realtocomp( A[ramp, shift : Ns_upRamp + shift ] * window )
    
    return RDc

def build_rd_from_data_rampDown( params, A, rampDownFirst=1, withHanning=0, withHamming=0 ):
    print( "build rd from ramp down" )
    print( "shape of the samples matrix = {}".format(A.shape) )
    nbRamps = A.shape[0]
    Ns_upRamp = int(params.samplesPerRamp / 2)
    Ns = int(params.samplesPerRamp / 4)

    RDc = np.zeros( (nbRamps, Ns) , dtype = complex )

    if rampDownFirst:
        print("ramp down first in the data files")
        shift = 0
    else:
        print("ramp up first in the data files")
        shift = Ns_upRamp
    
    if withHanning:
        print("with Hanning window")
        window = np.hanning( Ns_upRamp )
    elif withHamming:
        print("with Hamming window")
        window = np.hamming( Ns_upRamp )
    else:
        print("no window")
        window = np.ones( Ns_upRamp )

    for ramp in range(nbRamps):
        RDc[ ramp, :] = RD_realtocomp( np.flipud( A[ramp, shift : Ns_upRamp + shift ] ) * window )

    return RDc

def build_rd_from_data_sim( A, samplesPerRamp, shifted=0 ):
    print( "shape of the samples matrix = {}".format(A.shape) )
    nbRamps = A.shape[0]
    Ns_upRamp = int(samplesPerRamp / 2)
    Ns = int(samplesPerRamp / 4)

    RDc = np.zeros( (nbRamps, Ns) , dtype = complex )

    if shifted == 0:
        print("Data are NOT shifted")
        for ramp in range(nbRamps):
            RDc[ ramp, :] = RD_realtocomp( A[ramp, 0 : Ns_upRamp ] ) * 2 / 4096
    else:
        print("Data are shifted")
        for ramp in range(nbRamps):
            RDc[ ramp, :] = RD_realtocomp( A[ramp, Ns_upRamp : 2 * Ns_upRamp ] ) * 2 / 4096
    
    return RDc

def build_rd_from_data_upNdown( params, A, shifted=0, shift=0 ):
    print( "shape of the samples matrix = {}".format(A.shape) )
    nbRamps = A.shape[0]
    Ns_upRamp = int(params.samplesPerRamp / 2)
    Ns = int(params.samplesPerRamp / 4)

    if shift == 0:
        RDc = np.zeros( (nbRamps * 2, Ns) , dtype = complex ) # up ramps and down ramps
        if shifted == 0:
            print("up ramp first in the data files")
            for ramp in range(nbRamps):
                RDc[ 2 * ramp, :] = RD_realtocomp( A[ramp, 0 : Ns_upRamp ] ) * 2 / 4096
                RDc[ 2 * ramp+1, :] = RD_realtocomp( np.flipud( A[ramp, Ns_upRamp : 2 * Ns_upRamp ] ) ) * 2 / 4096
        else:
            print("down ramp first in the data files")
            for ramp in range(nbRamps):
                RDc[ 2 * ramp, :] = RD_realtocomp( np.flipud(A[ramp, 0 : Ns_upRamp ] ) ) * 2 / 4096
                RDc[ 2 * ramp+1, :] = RD_realtocomp( A[ramp, Ns_upRamp : 2 * Ns_upRamp ] ) * 2 / 4096
    else:
        print("down ramps are shifted by {}, first and last ramps removed for reconstruction".format(shift))
        RDc = np.zeros( (nbRamps * 2, Ns) , dtype = complex ) # up ramps and down ramps
        A_vec = A.reshape( nbRamps * params.samplesPerRamp )
        if shifted == 0:
            print("up ramp first in the data files")
            for ramp in range(1,nbRamps-1):
                sampleStart = ramp * params.samplesPerRamp
                s_up_0 = sampleStart
                s_up_1 = sampleStart + Ns_upRamp
                s_down_0 = sampleStart + Ns_upRamp + shift
                s_down_1 = sampleStart + 2 * Ns_upRamp + shift
                RDc[ 2 * (ramp), :] = \
                RD_realtocomp( A_vec[ s_up_0 : s_up_1 ] ) * 2 / 4096
                RDc[ 2 * (ramp)+1, :] = \
                RD_realtocomp( np.flipud( A_vec[ s_down_0 : s_down_1 ] ) ) * 2 / 4096
        else:
            print("NOT CODED YET")
    
    print( "shape of the RDc matrix = {}".format(RDc.shape) )

    return RDc

def read_data_and_build_rd(
    params,
    firstFileToRead,
    lastFileToRead,
    RDc,
    data_dir
):

    fileNumber = 0
    Ns = int(params.samplesPerRamp / 2)
    dum = np.zeros( params.samplesPerRamp * 2, dtype = np.int16 )
    timeSerie_A = np.zeros( Ns, dtype = np.int16 )

    for loop in range (firstFileToRead, lastFileToRead, params.blocksPerFile ):
        print( str(loop) + ' / ' + str(lastFileToRead - params.blocksPerFile) )
        # open the file containing data
        stream = data_dir + '/record' + str(loop) + '.bin'
        fd = open(stream,'rb')
        # get the data contained in the file ramp by ramp

        for ramp in range(params.rampsPerFile):
            rampNumber = fileNumber * params.rampsPerFile + ramp
            # get one data set adc_A, adc_B
            dum = np.fromfile(fd, dtype = np.int16, count = params.samplesPerRamp * 2)
            timeSerie_A = dum[ 0 : 2 * Ns : 2 ]
            RDc[ rampNumber, :] = RD_realtocomp( timeSerie_A ) * 2 / 4096

        fileNumber = fileNumber + 1

        fd.close()

    return timeSerie_A

def read_data_and_build_rd_shrink(
    params,
    firstFileToRead,
    lastFileToRead,
    data_dir,
    firstRange,
    lastRange,
    shifted=0
):

    fileNumber = 0
    Ns = int(params.samplesPerRamp / 2)
    dum = np.zeros( params.samplesPerRamp * 2, dtype = np.int16 )
    timeSerie_A = np.zeros( Ns, dtype = np.int16 )
    nbFiles = int ( (lastFileToRead - firstFileToRead) / params.blocksPerFile + 1)
    nbRamps = params.rampsPerFile * nbFiles
    RDc = RDc = np.zeros( (nbRamps, lastRange - firstRange ) , dtype = complex )

    print("nbFiles = {}, nbRamps = {}".format( nbFiles, nbRamps ) )

    if shifted == 0:
        print("Data are NOT shifted")
    else:
        print("Data are shifted")

    for loop in range (firstFileToRead, lastFileToRead+1, params.blocksPerFile ):
        print( "{} / {}".format( loop, lastFileToRead ) )
        # open the file containing data
        stream = data_dir + '/record' + str(loop) + '.bin'
        fd = open(stream,'rb')
        # get the data contained in the file ramp by ramp

        for ramp in range(params.rampsPerFile):
            rampNumber = fileNumber * params.rampsPerFile + ramp
            # get one data set adc_A, adc_B
            dum = np.fromfile(fd, dtype = np.int16, count = params.samplesPerRamp * 2)
            if shifted == 0:
                timeSerie_A = dum[ 0 : 2 * Ns : 2 ]
            else:
                timeSerie_A = dum[ 2 * Ns : 4 * Ns : 2 ]
            RDc[ rampNumber, :] = RD_realtocomp( timeSerie_A )[firstRange:lastRange] * 2 / 4096

        fileNumber = fileNumber + 1

        fd.close()

    return RDc
