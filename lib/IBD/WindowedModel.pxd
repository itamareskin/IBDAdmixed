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
    
    cdef public int _model_num
    cdef public int _model_partitions_num
    cdef public list _models
    cdef public list _models_by_partition
    cdef public list _inner_models
    
    cdef public int _snp_num
    cdef public int _win_size
    
    cdef double **_ems_prob
    
    cdef double **_forward_prob
    cdef double **_backward_prob
    
    cdef double *_scale_factor
    
    cdef double **_posterior_prob
    
    cdef double **_viterbi_prob
    cdef int **_viterbi_back_track
    cdef int *_viterbi_path

    cdef alloc_mem(self)
    
    cdef free_mem(self)
    
    cdef alloc_mem_viterbi(self)
    
    cdef free_mem_viterbi(self)
    
    cdef calc_ems_probs(self, TestSet obs_data, bool keep_inner=?)
    
    cdef calc_forward_probs(self)
    
    cdef calc_backward_probs(self)
    
    cpdef posterior_decoding(self, TestSet obs_data, bool keep_inner=?)
    
    cpdef viterbi_decoding(self, TestSet obs_data)
    
    cdef int get_num_windows(self)
    
    cdef int start_snp(self, int win_idx)
    
    cdef int end_snp(self, int win_idx)
    
    cpdef compare(self, WindowedModel other)

    cpdef compare_partitions(self, int idx1, int idx2)
