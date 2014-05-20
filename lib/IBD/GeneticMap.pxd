#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

cdef class GeneticMap(object):
    '''
    represents the genetic map data for a chromosome 
    '''
    
    cdef public int _snp_num
    
    # physical positions
    cdef int *_position
    
    # genetic distances (cM from 5' end)
    cdef double *_genetic_dist