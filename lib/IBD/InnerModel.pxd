#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

cdef class InnerModel(object):

    cpdef double calc_likelihood(self)