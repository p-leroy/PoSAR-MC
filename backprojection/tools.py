import numpy as np

c = 3e8

def pulse( x ):
    y = np.zeros( x.shape )
    y[ np.where( (-1/2<x) & (x<1/2) ) ] = 1
    return y

def pulse2( x ):
    if (-1/2<x) and (x<1/2):
        y = 1
    else:
        y = 0
    return y

def sb1( t, r, alpha, fc ):
    tau = 2 * r / c
    y = np.exp( -1j * 2  * np.pi * fc * tau ) \
            * np.exp( 1j * np.pi * alpha * tau**2 ) \
            * np.exp(  -1j * 2 * np.pi * alpha * t * tau )
    return y

def sb2( t, r, alpha, fc, T ):
    y = np.zeros( t.shape, dtype=complex )
    tau = 2 * r / c
    idx = np.where( t < tau )
    y[idx] = np.exp( -1j * 2  * np.pi * fc * (T-tau) ) \
            * np.exp( 1j * np.pi * alpha * (T-tau)**2 ) \
            * np.exp(  -1j * 2 * np.pi * alpha * t[idx] * (T-tau) )
    idx = np.where( t > tau )
    y[idx] = np.exp( -1j * 2  * np.pi * fc * tau ) \
            * np.exp( 1j * np.pi * alpha * tau**2 ) \
            * np.exp(  -1j * 2 * np.pi * alpha * (t[idx]-tau) * tau )
    return y

def sr( t, r ):
    
    tau = 2 * r / c
    
    y = np.zeros( t.shape, dtype=complex )
    
    # up ramp
    y[0:samplesPerUpRamp] = np.exp( -1j * 2  * np.pi * f0 * tau ) \
    * np.exp( 1j * np.pi * alpha * tau**2 ) \
    * np.exp(  -1j * 2 * np.pi * alpha * (t[0:samplesPerUpRamp]) * tau )
    
    # down ramp
    y[samplesPerUpRamp:2*samplesPerUpRamp] = np.exp( -1j * 2  * np.pi * f0 * tau ) \
    * np.exp( 1j * np.pi * (-alpha) * tau**2 ) \
    * np.exp(  -1j * 2 * np.pi * (-alpha) * (t[samplesPerUpRamp:2*samplesPerUpRamp]) * tau )
    
    return y

def srn( t, r, n ):
    
    tau = 2 * r / c
    
    y = np.zeros( t.shape, dtype=complex )
    
    # up ramp
    y[0:samplesPerUpRamp] = np.exp( -1j * 2  * np.pi * f0 * tau ) \
    * np.exp( 1j * np.pi * alpha * tau**2 ) \
    * np.exp(  -1j * 2 * np.pi * alpha * (t[0:samplesPerUpRamp] + n*T + tau) * tau )
    
    # down ramp
    y[samplesPerUpRamp:2*samplesPerUpRamp] = np.exp( -1j * 2  * np.pi * f0 * tau ) \
    * np.exp( 1j * np.pi * (-alpha) * tau**2 ) \
    * np.exp(  -1j * 2 * np.pi * (-alpha) * (t[samplesPerUpRamp:2*samplesPerUpRamp] + n*T + tau) * tau )
    
    return y

def wa( az, rg0, b=1 ):
    Rn = ( rg0**2 + az**2 )**0.5
    wa = np.sinc( b * np.arccos( rg0 / Rn ) )**2
    return wa
