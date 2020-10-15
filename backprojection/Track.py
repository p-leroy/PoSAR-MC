import json, os

import numpy as np

class Track(object):
    def __init__(self, conf, filename=None):
        self.conf = conf
        self.filename = filename
        self.load()

    def load(self):
        if self.filename is None:
            # use conf to determine the track filename automatically
            with open(os.path.join(self.conf.out_dir, "track_model.json")) as json_file:  
                data = json.load(json_file)
                self.ux = data['ux']
                self.uy = data['uy']
                self.refX = data['origX']
                self.refY = data['origY']
                self.theta = np.arctan2(self.ux[1], self.ux[0])
        else:
            with open(os.path.join(self.filename)) as json_file:  
                data = json.load(json_file)
                self.ux = data['ux']
                self.uy = data['uy']
                self.refX = data['origX']
                self.refY = data['origY']
                self.theta = np.arctan2(self.ux[1], self.ux[0])
