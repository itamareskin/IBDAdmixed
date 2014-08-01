#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

from libcpp cimport bool

cdef class GeneticMap(object):
    '''
    represents the genetic map data for a chromosome 
    '''

    cdef bool _is_slice

    cdef public int _snp_num
    
    # physical positions
    cdef long *_position
    
    # genetic distances (cM from 5' end)
    cdef double *_genetic_dist
    
    cpdef GeneticMap get_slice(self, int start_snp, int snp_num)
    
    cpdef dict get_position_dict(self)

    cpdef list get_genetic_dist_list(self)

    cpdef float get_length(self, int start_snp_idx, int end_snp_idx)
    