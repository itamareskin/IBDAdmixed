#cython: profile=True
#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

from InnerModel import InnerModel
from InnerModel cimport InnerModel
from TestSet import GenotypePair
from TestSet cimport GenotypePair
from LDModel import LDModel
from LDModel cimport LDModel
from libcpp cimport bool

cdef class GenotypePairModel(InnerModel):
    
    cdef bool _phased
    # the ibd state of the model (0 - no IBD, 1 - IBD between 1st and 3rd chromosomes)
    cdef int _ibd
    
    cdef LDModel _m1
    cdef LDModel _m2
    cdef LDModel _m3
    cdef LDModel _m4
    
    cdef double *****_forward_prob
    
    cdef double *****_backward_prob
    
    cdef double *****_emission_prob
    
    cdef double *_scale_factor
    cdef double *_backward_scale_factor
    
    cpdef calc_emission_probs(self, GenotypePair p)
    cpdef calc_forward_probs(self)
    cpdef calc_backward_probs(self)
    
    cdef rescale_forward(self, int snp_idx)
    cdef rescale_backward(self, int snp_idx)