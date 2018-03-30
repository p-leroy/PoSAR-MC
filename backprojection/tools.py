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

def sb1_r( t, r, alpha, fc ):
    tau = 2 * r / c
    y = np.cos( -2  * np.pi * fc * tau \
        + np.pi * alpha * tau**2 \
        -2 * np.pi * alpha * t * tau )
    return y

def sb1_tri( t, r, alpha, fc ):
    
    tau = 2 * r / c
    
    if t.size%2 != 0: 
        print("WARNING t.size is not a multiple of 2")

    y = np.zeros( t.shape, dtype=complex )
    samplesPerUpRamp = int( t.size / 2 )
    
    # up ramp
    y[0:samplesPerUpRamp] = np.exp( -1j * 2  * np.pi * fc * tau ) \
    * np.exp( 1j * np.pi * alpha * tau**2 ) \
    * np.exp(  -1j * 2 * np.pi * alpha * (t[0:samplesPerUpRamp]) * tau )
    
    # down ramp
    y[samplesPerUpRamp:2*samplesPerUpRamp] = np.exp( -1j * 2  * np.pi * fc * tau ) \
    * np.exp( 1j * np.pi * (-alpha) * tau**2 ) \
    * np.exp(  -1j * 2 * np.pi * (-alpha) * (t[samplesPerUpRamp:2*samplesPerUpRamp]) * tau )
    
    return y

def sb1_tri_r( t, r, alpha, fc ):
    
    tau = 2 * r / c

    if t.size%2 != 0: 
        print("WARNING t.size is not a multiple of 2")
    
    y = np.zeros( t.shape )
    samplesPerUpRamp = int( t.size / 2 )
    
    # up ramp
    y[0:samplesPerUpRamp] = np.cos( -2  * np.pi * fc * tau \
        + np.pi * alpha * tau**2 \
        -2 * np.pi * alpha * (t[0:samplesPerUpRamp]) * tau )
    
    # down ramp
    y[samplesPerUpRamp:2*samplesPerUpRamp] = np.cos( -12  * np.pi * fc * tau \
        + np.pi * (-alpha) * tau**2 \
        - 2 * np.pi * (-alpha) * (t[samplesPerUpRamp:2*samplesPerUpRamp]) * tau )
    
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

def wa( az, d0, b=1 ):
    d = ( d0**2 + az**2 )**0.5
    wa = np.sinc( b * np.arccos( d0 / d ) )**2
    return wa

# dx = lambda_c / (4 * sin( phi_a / 2 )) => phi_a = 2 * asin( lambda_c / (4 * dx))

def getApertureAngle( lambda_c, dx ):
    phi =  2 * np.arcsin( lambda_c / (4 * dx) )
    return phi

def getAzimuthResolution( lambda_c, phi ):
    dx = lambda_c / (4 * np.sin( phi / 2 ) )
    return dx