#cython: profile=True
#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

from InnerModel import InnerModel
from InnerModel cimport InnerModel
from TestSet import TestSet 
from TestSet cimport TestSet 
from TestSet import Genotype
from TestSet cimport Genotype
from LDModel import LDModel
from LDModel cimport LDModel
from libcpp cimport bool

cdef class GenotypeModel(InnerModel):
    
    cdef public bool _phased
    # the ibd state of the model (0 - no IBD, 1 - IBD between 1st and 3rd chromosomes)
    cdef public int _g
    
    cdef public LDModel _m1
    cdef public LDModel _m2
    
    cdef double ***_forward_prob
    
    cdef double ***_backward_prob
    
    cdef double ***_emission_prob
    
    cdef double *_scale_factor
    cdef double *_backward_scale_factor
    
    cpdef alloc_mem(self)
    cpdef free_mem(self)
    cpdef calc_emission_probs(self, Genotype p)
    cpdef calc_forward_probs(self)
    cpdef calc_backward_probs(self)
    
    cdef rescale_forward(self, int snp_idx)
    cdef rescale_backward(self, int snp_idx)

    cdef double anc_trans_prob(self, GenotypeModel other)
    cpdef double trans_prob(self, InnerModel other)
    
    cpdef print_inner_prob(self)
