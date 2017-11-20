from scipy.optimize import curve_fit
import numpy as np

R = 6378137  # approximate radius of earth

#startingPoint = ( 48.058403, -2.005964, 0.0 ) # top right hand corner of the runaway
startingPoint = ( 48 + 3 / 60 + 33.14 / 3600, 
	- (2 + 0 / 60 + 15.00 / 3600),
	0.0 )

startingPoint = (48.056830, -2.007717, 0.0)

# crossroad between monterfil north and south
startingPoint = ( 48 + 3 / 60 + 31.27 / 3600,
	- (2 + 0 / 60 + 21.62 / 3600),
	0.0 )

def func(x, a, b):
	return a * x + b

#startingPoint = ( 48.059200, -2.004171, 0.0 )
# (latitude, longitude, altitude)

def getxy( Lat, Long, orig ):

    x = R * ( (Long-orig[1]) * np.pi / 180 ) * np.cos( orig[0] * np.pi / 180 )
    y = R * ( (Lat-orig[0]) * np.pi / 180)

    return x, y

def getlatlong( x, y, orig ):
	pass

def getProjection( coord, proj, startingPoint, trackModel ):
	if len( coord.shape ) == 1:
		x, y = getxy( coord[0], coord[1], startingPoint )
		# project x and y in the track model
		proj[ 0 ] = x * trackModel["ux"][0] + y * trackModel["ux"][1]
		proj[ 1 ] = x * trackModel["uy"][0] + y * trackModel["uy"][1]
	else:
		for k in range( coord.shape[0] ):
			x, y = getxy( coord[k,0], coord[k,1], startingPoint )
			# project x and y in the track model
			proj[ k, 0 ] = x * trackModel["ux"][0] + y * trackModel["ux"][1]
			proj[ k, 1 ] = x * trackModel["uy"][0] + y * trackModel["uy"][1]
