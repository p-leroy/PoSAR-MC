import numpy as np

#EKF NAV
idx_timeStamp = 0
idx_velNorth = 1
idx_velEast = 2
idx_velDown = 3
idx_velNorth_StdDev = 4
idx_velEast_StdDev = 5
idx_velDown_StdDev = 6
idx_Lat = 7
idx_Long = 8
idx_Alt = 9
idx_undulation = 10
idx_Lat_StdDev = 11
idx_Long_StdDev = 12
idx_Alt_StdDev = 13
idx_nav_status = 14

track1_Lat_0  = 48.06069 
track1_Long_0 = -1.99354
track1_Lat_1  = 48.05507 
track1_Long_1 = -2.02359

track2_Lat_0  = 48.06249 
track2_Long_0 = -1.99467
track2_Lat_1  = 48.05687 
track2_Long_1 = -2.02434

track3_Lat_0  = 48.06555 
track3_Long_0 = -1.99619
track3_Lat_1  = 48.06007 
track3_Long_1 = -2.02550

runaway = np.array([ [ 48.057546, -2.010483, 0.0 ],
	[ 48.058403, -2.005964, 0.0 ],
	[ 48.058191, -2.005869, 0.0 ],
	[ 48.057327, -2.010383, 0.0 ],
	[ 48.057546, -2.010483, 0.0 ]
	])

hangar = np.array([ [ 48.056814, -2.007998, 0.0 ],
          [ 48.056822, -2.007793, 0.0 ],
          [ 48.056688, -2.007777, 0.0 ],
          [ 48.056679, -2.007978, 0.0 ],
          [ 48.056814, -2.007998, 0.0 ]
          ])

building = np.array([ [ 48.056830, -2.007717, 0.0 ],
          [ 48.056852, -2.007300, 0.0 ],
          [ 48.056771, -2.007296, 0.0 ],
          [ 48.056754, -2.007706, 0.0 ],
          [ 48.056830, -2.007717, 0.0 ]
          ])

concreteBlocks = np.array([
    [48.058067999622324,-2.0065217857185120,0.0],
    [48.057944804557310,-2.0071700279455684,0.0],
    [48.057820013472934,-2.0078100100003486,0.0],
    [48.057696817814474,-2.0084613762863315,0.0],
    [48.057576395106054,-2.0090977078737673,0.0],
    [48.057453198864690,-2.0097521982186746,0.0]
    ])

church = np.array( [48.066103, -1.978348, 82.27] )

monterfil_top_right = np.array([ 48.069233, -1.975694, 0.0 ])
monterfil_bottom_left = np.array([ 48.066013, -1.986860, 0.0 ])

cornerReflectorLarge = np.array( [48.057226, -2.008023, 0.0] )

cornerReflectorSmall = np.array( [48.057398, -2.008122, 0.0] )

def plotRunaway( ax ):
    ax.plot(runaway[:,1], runaway[:,0], "og", markeredgecolor = 'black')

def plotConcreteBlocks( ax ):
    ax.plot(concreteBlocks[:,1], concreteBlocks[:,0], "ow", markeredgecolor = 'black')    

def plotHangar( ax ):
    ax.plot(hangar[:,1], hangar[:,0], "og", markeredgecolor = 'black')

def plotBuilding( ax ):
    ax.plot(building[:,1], building[:,0], "og", markeredgecolor = 'black')

def plotLongLatAndTrackReferences( ax, Long, Lat ):
    ax.plot( Long, Lat, 'gray' )
    ax.plot( [track1_Long_0, track1_Long_1], [track1_Lat_0, track1_Lat_1], "o-.b", markeredgecolor = 'black' )
    ax.plot( [track2_Long_0, track2_Long_1], [track2_Lat_0, track2_Lat_1], "o-.b", markeredgecolor = 'black' )
    ax.plot( [track3_Long_0, track3_Long_1], [track3_Lat_0, track3_Lat_1], "o-.b", markeredgecolor = 'black' )

def addOnPlot(ax, Long, Lat, idx, color, label=''):
    idxStart = idx[0]
    idxStop = idx[1]
    ax.plot( Long[idxStart:idxStop], Lat[idxStart:idxStop], '--', color=color, label=label )
    ax.plot( Long[idxStart], Lat[idxStart], 'D', color=color, markeredgecolor = 'black' )
    ax.plot( Long[idxStop], Lat[idxStop], 'D', color=color, markeredgecolor = 'black'  )

def addSpot(ax, idx, color):
    ax.plot( Long[idx], Lat[idx], 'D', color=color, markeredgecolor = 'black' )

def getTimeOfDay( seconds ):
    h = np.floor( seconds / 3600 )
    m = np.floor( (seconds - h *3600) / 60 )
    s = seconds - h * 3600 - m * 60
    return (h, m, s)

def getSeconds( timeOfDay ):
    return (timeOfDay[0] * 3600 + timeOfDay[1] * 60 + timeOfDay[2])

def getIndices( start, stop, utc_seconds ):
    idx_start = np.amax( np.where( utc_seconds <= getSeconds(start) ) )
    idx_stop = np.amax( np.where( utc_seconds <= getSeconds(stop) ) )
    return (idx_start, idx_stop)

def getBlockNumberFromUtcSeconds( seconds, t_0 ):
    delta = seconds - getSeconds(t_0)
    return int( np.floor( delta / T_files ) * blocksPerFile )

def getAltitudeVelocity( Alt, Vel, idx0, idx1 ):
    altitude = np.average( Alt[ idx0 : idx1 ] )
    velocity = np.average( Vel[ idx0 : idx1 ] )
    print("(altitude, velocity) = ({:.3f}, {:.3f})".format(altitude, velocity) )

def getSelection(utc_seconds, Alt, Vel, idx):
    (idx0, idx1) = idx
    
    selection0 = getBlockNumberFromUtcSeconds( utc_seconds[idx0], t_0 )
    print( "Data selection starts at file: " + '{}'.format(selection0) )
    
    selection1 = getBlockNumberFromUtcSeconds( utc_seconds[idx1], t_0 )
    nb = int( np.floor( (selection1 - selection0) / blocksPerFile ) )
    print( "Data selection stops at file: " + '{}'.format(selection1) 
          + ", number of files in the selection: " + '{}'.format(nb+1) )
    
    getAltitudeVelocity( Alt, Vel, idx0, idx1 )

def getLastLine( filename ):
    with open(filename, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            last = next(reader)
        print(last)

def printUtc( index, utc ):
    print( str(utc[index][1]) 
          + " : " + str(utc[index][2]) 
          + " : " + "{:.3f}".format(utc[index][3] + utc[index][4] * 1e-9) )

def getDataAtIndex(filename, index):
    with open(filename, newline='') as f:
        reader = csv.reader(f)
        for r in range(index+1):
            next(reader)
        val = next(reader)
    return val

def getBlockNumber( t, t_0 ):
    delta = getSeconds( t ) - getSeconds(t_0)
    return int( np.floor( delta / T_files ) * blocksPerFile )

def getUtcIndex( tod ):
    return np.amax( np.where( utc_seconds <= getSeconds( tod ) ) )

def getUtcData( logUtcData, timeOfDay ):
    with open( logUtcData, newline='') as f:
        reader = csv.reader(f)
        next(reader)
        utc = np.loadtxt( io.StringIO( next(reader)[0] ))
        
        while (utc[5] < timeOfDay[0] - hourOffset):
            utc = np.loadtxt( io.StringIO( next(reader)[0] ))
            
        while (utc[6] < timeOfDay[1]):
            utc = np.loadtxt( io.StringIO( next(reader)[0] ))
            
        while (utc[7] < timeOfDay[2]):
            utc = np.loadtxt( io.StringIO( next(reader)[0] ))
            
    return utc
