import _ctypes
from ctypes import *
from numpy import ctypeslib

class MyComplex(Structure):
    _fields_ = [
    ('real', c_double),
    ('imag', c_double) ]

class MyParameters(Structure):
    _fields_ = [
    ('Nx', c_int),
    ('Ny', c_int),
    ('Nover', c_int),
    ('dx', c_double),
    ('Naz', c_int),
    ('Nf', c_int),
    ('hScene', c_double),
    ('phi_a_deg', c_double) ]

class MyParametersPoSAR_GB(Structure):
    _fields_ = [
    ('Nx', c_int),
    ('Ny', c_int),
    ('Nz', c_int),
    ('Nover', c_int),
    ('dx', c_double),
    ('Naz', c_int),
    ('Nf', c_int),
    ('phi_a_deg', c_double),
    ('kc', c_double) ]

class MyParameters_LETG(Structure):
    _fields_ = [
    ('Nx', c_int),
    ('Ny', c_int),
    ('Nover', c_int),
    ('dx', c_double),
    ('Naz', c_int),
    ('Nf', c_int),
    ('hScene', c_double),
    ('phi_a_deg', c_double),
    ('uxx', c_double),
    ('uxy', c_double),
    ('meanX', c_double),
    ('meanY', c_double),
    ('kc', c_double)]

class Lib(object):
  """docstring for LibBackProjection"""
  def __init__(self, filename):
    super(Lib, self).__init__()
    self.filename = filename
    self.so = self.load()
    self.setArgTypesResType()

  def setArgTypesResType(self):

    # backProjection2
    self.so.backProjection2.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int, c_double,
    ctypeslib.ndpointer(complex,  ndim=1, flags='C'), c_int, c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex,  ndim=1, flags='C') ]

    # backProjectionOmp
    self.so.backProjectionOmp.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int, c_double,
    ctypeslib.ndpointer(complex,  ndim=1, flags='C'), c_int, c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex,  ndim=1, flags='C') ]

    # backProjectionOmpSlantRange
    self.so.backProjectionOmpSlantRange.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex,  ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    MyParameters ]

    # backProjectionOmpGroundRange
    self.so.backProjectionOmpGroundRange.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex,  ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    MyParameters ]

    # backProjectionOmpGroundRange_LETG
    self.so.backProjectionOmpGroundRange_LETG.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # x
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # y
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # z
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex,  ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    MyParameters_LETG ]

    # backProjectionOmpGroundRange_corr
    self.so.backProjectionOmpGroundRange_corr.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # x
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # y
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # z
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # r_over
    ctypeslib.ndpointer(complex,  ndim=1, flags='C'), # sr
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex, ndim=1, flags='C'),  # img
    MyParameters_LETG ]

    # backProjectionOmpGroundRange_PoSAR_GB
    self.so.backProjectionOmpGroundRange_PoSAR_GB.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # x
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # y
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # z
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex,  ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    MyParametersPoSAR_GB ]

    # backProjectionOmpGroundRange_PoSAR_GBalt
    self.so.backProjectionOmpGroundRange_PoSAR_GBalt.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # x
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # y
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # z
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex,  ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    MyParametersPoSAR_GB ]

    # backProjectionOmpGroundRange_PoSAR_GB_a
    self.so.backProjectionOmpGroundRange_PoSAR_GB_a.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # x
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # y
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # z
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex,  ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    MyParametersPoSAR_GB ]

    # backProjectionOmpGroundRange_PoSAR_GB_lha
    self.so.backProjectionOmpGroundRange_PoSAR_GB_lha.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # x
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # y
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), # z
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex,  ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    MyParametersPoSAR_GB ]

    # backProjectionOmpGroundRange_NED
    self.so.backProjectionOmpGroundRange_NED.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex,  ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    MyParameters,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C') ]

    # backProjectionOmpGroundRangeb
    self.so.backProjectionOmpGroundRangeb.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex,  ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    MyParameters ]

    # measureAndSavePlans
    self.so.measureAndSavePlans.argtypes = [
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    c_int,
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    c_int ]

    # measurePlans
    self.so.measurePlans.argtypes = [
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    c_int,
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    c_int ]
  
    # resample
    self.so.resample.argtypes = [
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    c_int,
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    c_int ]

    # zeroPaddingAndIfft4
    self.so.zeroPaddingAndIfft4.argtypes = [
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    c_int,
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    c_int ]

    # zeroPaddingAndIfft4b
    self.so.zeroPaddingAndIfft4b.argtypes = [
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    c_int,
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    c_int ]

    # pulse
    self.so.pulse.argtypes = [c_double]

    # interp
    self.so.myInterp.argtypes = [
    c_double, 
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    c_double ]
    self.so.myInterp.restype = MyComplex

  def reload(self):
    _ctypes.dlclose(self.so._handle)
    self.so = cdll.LoadLibrary(self.filename)
    self.so.importPlans()
    self.setArgTypesResType()

  def load(self):
    so = cdll.LoadLibrary(self.filename)
    so.importPlans()
    return so
