#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

cdef class InnerModel(object):

    def __cinit__(self):
        pass
    
    cpdef double calc_likelihood(self):
        pass