import numpy as np
import matplotlib.pyplot as plt

import posarmctools.epsgtools as epsg

class Scene(object):
    def __init__(self, conf, proj, track=None, target_ll=None, xyz=None):
        self.conf = conf
        self.track = track
        self.xyz = xyz

        if target_ll is None:
            self.targetX = track.refX
            self.targetY = track.refY
        else:
            self.targetLat = target_ll[0]
            self.targetLong = target_ll[1]
            self.targetX, self.targetY = epsg.wgs84LongLatToEpsg((self.targetLong, self.targetLat), conf.proj)
            if xyz is not None:
                self.xa = xyz[:,2]
                self.ya = xyz[:,3]
                self.za = xyz[:,4]

        if conf.groundRange == 1:
            self.buildSceneGroundRange(conf, track, target_ll)
        else:
            self.hFocusing = self.buildSceneSlantRange(conf, track)

    def buildSceneGroundRange(self, conf, track, target_ll):
        print(f"buildSceneGroundRange")
        if track is not None:
            if target_ll is None:
                self.baseX = track.refX + np.cos(track.theta) * np.arange(conf.xMax, conf.xMin, -conf.d_x)
                self.baseY = track.refY + np.sin(track.theta) * np.arange(conf.xMax, conf.xMin, -conf.d_x)
            else:
                self.baseX = self.targetX + np.cos(track.theta) * np.arange(conf.xMax, conf.xMin, -conf.d_x)
                self.baseY = self.targetY + np.sin(track.theta) * np.arange(conf.xMax, conf.xMin, -conf.d_x)
            self.nX = self.baseX.size
            y = np.arange(conf.yMin, conf.yMax, conf.d_y)
            self.nY = y.size
            self.X = np.zeros((self.nY, self.nX))
            self.Y = np.zeros((self.nY, self.nX))
            for line in range(self.nY):
                self.X[line, :] = self.baseX + track.uy[0] * y[line]
                self.Y[line, :] = self.baseY + track.uy[1] * y[line]
        else:
            if target_ll is None:
                print("/!\\ ERROR unexpected configuration for Scene definition /!\\")
            else:
                vec_x = self.targetX + np.arange(conf.xMin, conf.xMax, conf.d_x)
                vec_y = self.targetY + np.arange(conf.yMin, conf.yMax, conf.d_y)
                self.nX = vec_x.size
                self.nY = vec_y.size
                self.X, self.Y = np.meshgrid(vec_x, vec_y)
                self.long, self.lat = epsg.epsgToWgs84LongLat((self.X, self.Y), conf.proj)

        self.long, self.lat = epsg.epsgToWgs84LongLat((self.X, self.Y), conf.proj)
        self.X_mean = np.mean(self.X)
        self.Y_mean = np.mean(self.Y)

    def buildSceneSlantRange(self, conf, track):
        print(f"buildSceneSlantRange")
        A = (track.refX, track.refY)
        C = (self.targetX, self.targetY)
        AC = (C[0] - A[0], C[1] - A[1])
        absAC = (AC[0]**2 + AC[1]**2)**0.5
        AC_dot_Ux = (AC[0] * track.ux[0] + AC[1] * track.ux[1])
        self.refX = A[0] + AC_dot_Ux * track.ux[0]
        self.refY = A[1] + AC_dot_Ux * track.ux[1]
        self.closestApproachXY = (absAC**2 - AC_dot_Ux**2)**0.5
        distanceRefTrackXY = ((self.refX - self.xa)**2 + (self.refY - self.ya)**2)**0.5
        idx = np.where(distanceRefTrackXY == np.amin(distanceRefTrackXY))[0][0]
        hFocusing = self.za[idx] - conf.hScene
        self.closestApproach = (self.xa[idx], self.ya[idx], self.za[idx])
        self.closestApproachXYZ = ((self.xa[idx] - self.targetX)**2 + (self.ya[idx] - self.targetY)**2 + hFocusing**2)**0.5
        self.closestApproachAngle = np.arctan(hFocusing / self.closestApproachXY) * 180 / np.pi
        print(f"closestApproachXY {self.closestApproachXY:.1f}, closestApproachXYZ {self.closestApproachXYZ:.1f}")
        print(f"closestApproachAngle {self.closestApproachAngle:.1f}, hFocusing {hFocusing:.1f}")
        self.x = np.arange(conf.xMax, conf.xMin, -conf.d_x)
        self.y = np.arange(conf.yMin, conf.yMax, conf.d_y) + self.closestApproachXYZ
        self.baseX = self.refX + np.cos(track.theta) * self.x
        self.baseY = self.refY + np.sin(track.theta) * self.x
        self.nX = self.baseX.size
        self.nY = self.y.size
        self.X = np.zeros((self.nY, self.nX))
        self.Y = np.zeros((self.nY, self.nX))
        for line in range(self.nY):
            self.X[line, :] = self.baseX + track.uy[0] * self.y[line]
            self.Y[line, :] = self.baseY + track.uy[1] * self.y[line]

        self.long, self.lat = epsg.epsgToWgs84LongLat((self.X, self.Y), conf.proj)
        self.X_mean = np.mean(self.X)
        self.Y_mean = np.mean(self.Y)
        return hFocusing

    def plotScene(self, sceneZ):
        fig, ax = plt.subplots(1,1)
        cs = ax.contourf(
            self.X, 
            self.Y, 
            sceneZ,
            vmin=0,
            vmax=130,
            cmap="terrain")
        ax.contour(cs, colors='k')
        ax.grid()
        ax.set_title("scene elevation " + self.conf.data_date)
        #dia.addColorBar(cs, ax)
