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
    ('Nr', c_int),
    ('Nover', c_int),
    ('dx', c_double),
    ('Naz', c_int),
    ('Nf', c_int),
    ('hScene', c_double),
    ('phi_a_deg', c_double) ]

class LibBackProjection(object):
  """docstring for LibBackProjection"""
  def __init__(self, filename):
    super(LibBackProjection, self).__init__()
    self.filename = filename
    self.so = self.load()
    self.setArgTypesResType()

  def setArgTypesResType(self):

    # backProjection2
    self.so.backProjection2.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int, c_double,
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), c_int, c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex, ndim=1, flags='C') ]

    # backProjectionOmp
    self.so.backProjectionOmp.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int, c_double,
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), c_int, c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex, ndim=1, flags='C') ]

    # backProjectionOmpGroundRange
    self.so.backProjectionOmpGroundRange.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int, c_double,
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), c_int, c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    c_double ]

    # call_backProjectionOmpGroundRange
    self.so.call_backProjectionOmpGroundRange.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    MyParameters ]

    # backProjectionOmpGroundRangeb
    self.so.backProjectionOmpGroundRangeb.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int, c_double,
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), c_int, c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    c_double ]

    # backProjectionOmpGroundRange_v2
    self.so.backProjectionOmpGroundRange_v2.argtypes = [
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), c_int, c_double,
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), c_int, c_int,
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'), ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    c_double ]

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

    # resample4
    self.so.resample4.argtypes = [
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    c_int,
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    c_int ]

    # resample4b
    self.so.resample4b.argtypes = [
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    c_int,
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    ctypeslib.ndpointer(complex, ndim=1, flags='C'),
    c_int ]

    # pulse
    self.so.pulse.argtypes = [c_double]
    self.so.interp.restype = c_double

    # interp
    self.so.interp.argtypes = [
    c_double, 
    ctypeslib.ndpointer(c_double, ndim=1, flags='C'),
    ctypeslib.ndpointer(complex, ndim=1, flags='C'), 
    c_double ]
    
    self.so.interp.restype = MyComplex

  def reload(self):
    _ctypes.dlclose(self.so._handle)
    self.so = cdll.LoadLibrary(self.filename)
    self.so.importPlans()
    self.setArgTypesResType()

  def load(self):
    so = cdll.LoadLibrary(self.filename)
    so.importPlans()
    return so
