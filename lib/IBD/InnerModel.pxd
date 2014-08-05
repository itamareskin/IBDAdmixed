from TestSet import TestSet
from TestSet cimport TestSet
from libcpp cimport bool

cdef class InnerModel(object):

    # initial probability of this model (use in Windowed Model)
    cdef double _prior
    
    cpdef InnerModel slice_from_model(self, int start_snp, int snp_num)
    
    cpdef double trans_prob(self, InnerModel other)
    
    cpdef double calc_likelihood(self, TestSet obs_data, bool free_mem=?)