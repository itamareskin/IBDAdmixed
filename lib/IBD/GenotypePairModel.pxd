from InnerModel import InnerModel
from InnerModel cimport InnerModel
from TestSet import TestSet
from TestSet cimport TestSet
from TestSet import GenotypePair
from TestSet cimport GenotypePair
from LDModel import LDModel
from LDModel cimport LDModel
from libcpp cimport bool

cdef inline double single_anc_trans_prob(LDModel first, LDModel second):
    if first._anc == second._anc:
        return 1
    else:
        return first._alpha * second._alpha

cdef class GenotypePairModel(InnerModel):
    
    cdef public bool _phased
    # the ibd state of the model (0 - no IBD, 1 - IBD between 1st and 3rd chromosomes)
    cdef public int _ibd
    
    cdef public int _g
    
    cdef public LDModel _m1
    cdef public LDModel _m2
    cdef public LDModel _m3
    cdef public LDModel _m4
    
    cdef double *****_forward_prob
    
    cdef double *****_backward_prob
    
    cdef double *****_emission_prob
    
    cdef double *_scale_factor
    cdef double *_backward_scale_factor
    
    cpdef alloc_mem(self)
    cpdef free_mem(self)
    cpdef calc_emission_probs(self, GenotypePair p)
    cpdef calc_forward_probs(self)
    cpdef calc_backward_probs(self)
    
    cdef rescale_forward(self, int snp_idx)
    cdef rescale_backward(self, int snp_idx)
    
    cpdef double ibd_trans_prob(self, GenotypePairModel other)
    cpdef double anc_trans_prob(self, GenotypePairModel other)
    cpdef double trans_prob(self, InnerModel other)
    
    cpdef print_inner_prob(self)
