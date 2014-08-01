#cython: cdivision=True

import operator
import string

cdef extern from "stdlib.h":
    int ceil(float f)
    float log(float f)
    int RAND_MAX
    int rand()
    int strlen(char *)
    int iabs(int)

cdef class IntervalNode:
    """
    A single node of an `IntervalTree`.

    """
    cdef float priority
    cdef public object interval
    cdef public int start, end
    cdef int minend, maxend, minstart
    cdef IntervalNode cleft, cright, croot

    cpdef IntervalNode insert(IntervalNode self, int start, int end, object interval)

    cdef IntervalNode rotate_right(IntervalNode self)

    cdef IntervalNode rotate_left(IntervalNode self)

    cdef inline void set_ends(IntervalNode self)

    cdef void _intersect( IntervalNode self, int start, int end, list results)

    cdef void _seek_left(IntervalNode self, int position, list results, int n, int max_dist)

    cdef void _seek_right(IntervalNode self, int position, list results, int n, int max_dist)

    cpdef left(self, position, int n=?, int max_dist=?)

    cpdef right(self, position, int n=?, int max_dist=?)

    cdef void _traverse(IntervalNode self, object func)

cdef IntervalNode EmptyNode

cdef class Interval:
    """
    Basic feature, with required integer start and end properties.
    Also accepts optional value

    """
    cdef public int start, end
    cdef public float value

cdef class IntervalTree:
    """
    Data structure for performing window intersect queries on a set of
    of possibly overlapping 1d intervals.

    """

    cdef IntervalNode root