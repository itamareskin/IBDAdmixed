#cython: profile=True
#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

from InnerModel import InnerModel
from InnerModel cimport InnerModel
from libcpp cimport bool

cdef class WindowedModel(object):
    
