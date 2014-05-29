#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

from TestSet import TestSet
from TestSet cimport TestSet

cdef class InnerModel(object):

    def __cinit__(self):
        self._prior = 1
    
    cpdef InnerModel slice_from_model(self, int start_snp, int snp_num):
        pass
    
    cpdef double trans_prob(self, InnerModel other):
        pass
    
    cpdef double calc_likelihood(self, TestSet g):
        pass