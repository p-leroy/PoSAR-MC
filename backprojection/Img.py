import configparser, os

import numpy as np
import scipy

class Img(object):
    def __init__(self, img):
        self.img = img
        self.dB = 20 * np.log10(np.abs(img))
        self.min_dB = np.amin(self.dB)
        self.max_dB = np.amax(self.dB)
        self.med_dB = np.median(self.dB)
        print("min_dB = {:.2f}, max_dB = {:.2f}, med_dB = {:.2f}".format(self.min_dB, self.max_dB, self.med_dB))

    def dynamicRange(self, dynRange=30, db=1, size=3):
        # box = 3 gives the same result as Nsmooth=2 with the box_filter in posarutils.process.filtering
        if db:
            idx = np.where(self.img != 0)
            imgMin = np.amin(20 * np.log10(np.abs(self.img[idx])))
            idx = np.where(self.img == 0)
            self.img[idx] = imgMin
            if size:
                z = self.box_dB(size)
            else:
                z = self.img_dB
        else:
            if size:
                z = self.box(size)
            else:
                z = np.abs(self.img)
        
        med = np.median(z)
        toPlot = np.minimum(np.maximum(z, med - dynRange / 2), med + dynRange / 2)
        print(f"min {np.amin(toPlot):.1f}, max {np.amax(toPlot):.1f}")

        return toPlot

    def box(self, n):
        y = np.zeros(self.img.shape)
        scipy.ndimage.filters.uniform_filter(np.abs(self.img), n, y, mode='constant')
        return y

    def box_dB(self, n):
        return 20 * np.log10(self.box(n))

    def save(self, focusing):
        np.save(os.path.join(focusing.out_dir, focusing.conf.name), self.img)

    def saveAlt(self, focusing, dtype='complex64'):
        filename = os.path.join(focusing.out_dir, focusing.conf.name + '.data')
        print("save")
        with open(filename, 'w') as f:
            print(filename)
            self.img.astype(dtype).tofile(f)

        config = configparser.ConfigParser()
        # data
        config['data'] = {}
        config['data']['shape'] = str(self.img.shape)
        config['data']['dtype'] = dtype
        # parameters
        config['parameters'] = {}
        config['parameters']['phi_a_deg'] = str(focusing.conf.phi_a_deg)
        config['parameters']['azimuth_resolution'] = f"{focusing.azimuthResolution:.3f}"
        config['parameters']['range_resolution'] = f"{focusing.conf.c / (2 * focusing.conf.B0)}"
        # focusing
        config['focusing'] = {}
        if focusing.conf.groundRange:
            config['focusing']['type'] = "ground range"
        else:
            config['focusing']['type'] = "slant range"
        config['focusing']['window'] = focusing.conf.window
        config['focusing']['fmcw_ramp'] = focusing.conf.upOrDown
        # scene
        config['scene'] = {}
        config['scene']["projection"] = focusing.conf.projStr
        if focusing.conf.groundRange:
            config['scene']['hScene'] = f"{focusing.conf.hScene:.3f}"
            config['scene']['scalex'] = f"{focusing.conf.d_x * (-np.cos(focusing.track.theta)):.3f}"
            config['scene']['skewx'] = f"{focusing.conf.d_x * (-np.sin(focusing.track.theta)):.3f}"
            config['scene']['skewy'] = f"{focusing.conf.d_y * (-np.sin(focusing.track.theta)):.3f}"
            config['scene']['scaley'] = f"{focusing.conf.d_y * (+np.cos(focusing.track.theta)):.3f}"
        else:
            config['scene']['focusing_height'] = f"{focusing.scene.hFocusing:.3f}"
            config['scene']['d_x'] = str(focusing.conf.d_x)
            config['scene']['d_y'] = str(focusing.conf.d_y)
            config['scene']['xMin'] = f"{np.amin(focusing.scene.x):.3f}"
            config['scene']['xMax'] = f"{np.amax(focusing.scene.x):.3f}"
            config['scene']['yMin'] = f"{np.amin(focusing.scene.y):.3f}"
            config['scene']['yMax'] = f"{np.amax(focusing.scene.y):.3f}"
            config['scene']['origin_x'] = f"{focusing.scene.refX:.3f}"
            config['scene']['origin_y'] = f"{focusing.scene.refY:.3f}"
            config['scene']['theta_rad'] = f"{focusing.track.theta:.3f}"
        filename = os.path.join(focusing.out_dir, focusing.conf.name + '.ini')
        with open(filename, 'w') as configfile:
            print(filename)
            config.write(configfile)