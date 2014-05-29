#cython: profile=True
#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

from InnerModel import InnerModel
from InnerModel cimport InnerModel
from TestSet import TestSet
from TestSet cimport TestSet
from libcpp cimport bool

cdef class WindowedModel(object):
    
    cdef int _model_num
    cdef list _models
    
    cdef int _snp_num
    cdef int _win_size
    
    cdef double **_ems_prob
    
    cdef double **_forward_prob
    cdef double **_backward_prob
    
    cdef double *_scale_factor
    
    cdef double *_posterior_prob
    
    cdef double **_viterbi_prob
    cdef int **_viterbi_back_track
    cdef int *_viterbi_path
 
    cpdef alloc_mem(self)
    
    cpdef free_mem(self)
    
    cpdef alloc_mem_viterbi(self)
    
    cpdef free_mem_viterbi(self)
    
    cpdef calc_ems_probs(self, TestSet obs_data)
    
    cpdef calc_forward_probs(self)
    
    cpdef calc_backward_probs(self)
    
    cpdef posterior_decoding(self)
    
    cpdef viterbi_decoding(self)
    
    cdef int get_num_windows(self)
    
    cdef int start_snp(self, int win_idx)
    
    cdef int end_snp(self, int win_idx)
    
    cpdef compare(self, WindowedModel other)
