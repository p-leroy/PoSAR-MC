import xml.etree.ElementTree as etree
import numpy as np

class PosarMCParameters():
    def __init__(self, filename):
        
        self.filename = filename
        
        tree = etree.parse( filename )
        
        # adf4158
        self.requestedTRamp = float( tree.find( './/requestedTRamp' ).attrib["value"])
        self.configuredTRamp = float( tree.find( './/configuredTRamp' ).attrib["value"] )
        self.startFrequency = float( tree.find( './/startFrequency' ).attrib["value"] )
        self.stopFrequency = float( tree.find( './/stopFrequency' ).attrib["value"] )
        self.frequencyBand = float( tree.find( './/frequencyBand' ).attrib["value"] )
        self.frequencyDevReq = float( tree.find( './/frequencyDevReq' ).attrib["value"] )
        self.frequencyDevConf = float( tree.find( './/frequencyDevConf' ).attrib["value"] )
        self.numberOfSteps = int( tree.find( './/numberOfSteps' ).attrib["value"] )
        self.waveformType = int( tree.find( './/waveformType' ).attrib["value"] )
        
        #posarV1
        self.rampsPerBlock = int( tree.find( './/rampsPerBlock' ).attrib["value"] )
        self.blockSize = float( tree.find( './/blockSize' ).attrib["value"] )
        self.blocksPerFile = int( tree.find( './/blocksPerFile' ).attrib["value"] )
        self.rampsPerFile = int( tree.find( './/rampsPerFile' ).attrib["value"] )
        self.fileSize = float( tree.find( './/fileSize' ).attrib["value"] )
        self.samplingFrequency = float( tree.find( './/samplingFrequency' ).attrib["value"] )
        self.samplesPerRamp = int( tree.find( './/samplesPerRamp' ).attrib["value"] )
        self.skipNSamples = int( tree.find( './/skipNSamples' ).attrib["value"] )

class PosarMCParameters_v2():
    def __init__(self, filename):
        
        self.filename = filename
        
        tree = etree.parse( filename )
        
        # adf4158
        self.requestedTRamp = float( tree.find( './/requestedTRamp' ).attrib["value"])
        self.configuredTRamp = float( tree.find( './/configuredTRamp' ).attrib["value"] )
        self.startFrequency = float( tree.find( './/startFrequency' ).attrib["value"] )
        self.stopFrequency = float( tree.find( './/stopFrequency' ).attrib["value"] )
        self.frequencyBand = float( tree.find( './/frequencyBand' ).attrib["value"] )
        self.frequencyDevReq = float( tree.find( './/frequencyDevReq' ).attrib["value"] )
        self.frequencyDevConf = float( tree.find( './/frequencyDevConf' ).attrib["value"] )
        self.numberOfSteps = int( tree.find( './/numberOfSteps' ).attrib["value"] )
        self.waveformType = int( tree.find( './/waveformType' ).attrib["value"] )
        
        #posarV1
        self.rampsPerBuffer = int( tree.find( './/rampsPerBuffer' ).attrib["value"] )
        self.bufferSize = float( tree.find( './/bufferSize' ).attrib["value"] )
        self.buffersPerFile = int( tree.find( './/buffersPerFile' ).attrib["value"] )
        self.rampsPerFile = int( tree.find( './/rampsPerFile' ).attrib["value"] )
        self.fileSize = float( tree.find( './/fileSize' ).attrib["value"] )
        self.samplingFrequency = float( tree.find( './/samplingFrequency' ).attrib["value"] )
        self.samplesPerRamp = int( tree.find( './/samplesPerRamp' ).attrib["value"] )
        self.skipNSamples = int( tree.find( './/skipNSamples' ).attrib["value"] )
        #
        self.kc = 4 * np.pi / 3e8 * ( self.startFrequency + self.stopFrequency ) / 2

    def print(self):
        # adf4158
        print( "requestedTRamp {}".format(self.requestedTRamp) )
        print( "configuredTRamp {}".format(self.configuredTRamp) )
        print( "startFrequency {}".format(self.startFrequency) )
        print( "stopFrequency {}".format(self.stopFrequency) )
        print( "frequencyBand {}".format(self.frequencyBand) )
        print( "frequencyDevReq {}".format(self.frequencyDevReq) )
        print( "frequencyDevConf {}".format(self.frequencyDevConf) )
        print( "numberOfSteps {}".format(self.numberOfSteps) )
        print( "waveformType {}".format(self.waveformType) )
        
        #posarV1
        print( "rampsPerBuffer {}".format(self.rampsPerBuffer) )
        print( "bufferSize {}".format(self.bufferSize) )
        print( "buffersPerFile {}".format(self.buffersPerFile) )
        print( "rampsPerFile {}".format(self.rampsPerFile) )
        print( "fileSize {}".format(self.fileSize) )
        print( "samplingFrequency {}".format(self.samplingFrequency) )
        print( "samplesPerRamp {}".format(self.samplesPerRamp) )
        print( "skipNSamples {}".format(self.skipNSamples) )
        print( f"kc {kc:.3f}" )

class PosarXParameters():
    def __init__(self, filename):
        
        self.filename = filename
        
        tree = etree.parse( filename )
        
        #posarX
        self.rampsPerBuffer = int( tree.find( './/rampsPerBuffer' ).attrib["value"] )
        self.bufferSize = float( tree.find( './/bufferSize' ).attrib["value"] )
        self.buffersPerFile = int( tree.find( './/buffersPerFile' ).attrib["value"] )
        self.rampsPerFile = int( tree.find( './/rampsPerFile' ).attrib["value"] )
        self.fileSize = float( tree.find( './/fileSize' ).attrib["value"] )
        self.samplingFrequency = float( tree.find( './/samplingFrequency' ).attrib["value"] )
        self.samplesPerRamp = int( tree.find( './/samplesPerRamp' ).attrib["value"] )
        self.skipNSamples = int( tree.find( './/skipNSamples' ).attrib["value"] )
        #
        self.rampPeriod = int( tree.find( './/rampPeriod' ).attrib["value"] )
        self.startFrequency = int( tree.find( './/startFrequency' ).attrib["value"] )
        self.stopFrequency = int( tree.find( './/stopFrequency' ).attrib["value"] )
        #
        self.kc = 4 * np.pi / 3e8 * ( self.startFrequency + self.stopFrequency ) / 2 * 1e6

    def print(self):
        #posarX
        print( "rampsPerBuffer {}".format(self.rampsPerBuffer) )
        print( "bufferSize {}".format(self.bufferSize) )
        print( "buffersPerFile {}".format(self.buffersPerFile) )
        print( "rampsPerFile {}".format(self.rampsPerFile) )
        print( "fileSize {}".format(self.fileSize) )
        print( "samplingFrequency {}".format(self.samplingFrequency) )
        print( "samplesPerRamp {}".format(self.samplesPerRamp) )
        print( "skipNSamples {}".format(self.skipNSamples) )
        #
        print( "rampPeriod {}".format(self.rampPeriod) )
        print( "startFrequency {}".format(self.startFrequency) )
        print( "stopFrequency {}".format(self.stopFrequency) )
        #
        print( f"kc {self.kc:.3f}" )
