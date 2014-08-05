#cython: profile=True
#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

import operator
import string

cdef extern from "stdlib.h":
    int ceil(float f)
    float log(float f)
    int RAND_MAX
    int rand()
    int strlen(char *)
    int iabs(int)

cdef inline int imax2(int a, int b):
    if b > a: return b
    return a

cdef inline int imax3(int a, int b, int c):
    if b > a:
        if c > b:
            return c
        return b
    if a > c:
        return a
    return c

cdef inline int imin3(int a, int b, int c):
    if b < a:
        if c < b:
            return c
        return b
    if a < c:
        return a
    return c

cdef inline int imin2(int a, int b):
    if b < a: return b
    return a

cdef float nlog = -1.0 / log(0.5)

def string2IntervalTree(s):
    t = IntervalTree()
    l = [string.split(x,",") for x in string.split(s,sep=";")]
    l = [Interval(int(x[0]),int(x[1])) for x in l]
    for interval in l:
        t.add_interval(interval)
    return t

cdef class IntervalNode:
    """
    A single node of an `IntervalTree`.

    """

    property left_node:
        def __get__(self):
            return self.cleft if self.cleft is not EmptyNode else None
    property right_node:
        def __get__(self):
            return self.cright if self.cright is not EmptyNode else None
    property root_node:
        def __get__(self):
            return self.croot if self.croot is not EmptyNode else None
    
    def __repr__(self):
        return "IntervalNode(%i, %i)" % (self.start, self.end)

    def __cinit__(IntervalNode self, int start, int end, object interval):
        # Python lacks the binomial distribution, so we convert a
        # uniform into a binomial because it naturally scales with
        # tree size.  Also, python's uniform is perfect since the
        # upper limit is not inclusive, which gives us undefined here.
        self.priority = ceil(nlog * log(-1.0/(1.0 * rand()/(RAND_MAX+1) - 1)))
        self.start    = start
        self.end      = end
        self.interval = interval
        self.maxend   = end
        self.minstart = start
        self.minend   = end
        self.cleft    = EmptyNode
        self.cright   = EmptyNode
        self.croot    = EmptyNode
        
    cpdef IntervalNode insert(IntervalNode self, int start, int end, object interval):
        """
        Insert a new IntervalNode into the tree of which this node is
        currently the root. The return value is the new root of the tree (which
        may or may not be this node!)
        """
        cdef IntervalNode croot = self
        # If starts are the same, decide which to add interval to based on
        # end, thus maintaining sortedness relative to start/end
        cdef int decision_endpoint = start
        if start == self.start:
            decision_endpoint = end
        
        if decision_endpoint > self.start:
            # insert to cright tree
            if self.cright is not EmptyNode:
                self.cright = self.cright.insert( start, end, interval )
            else:
                self.cright = IntervalNode( start, end, interval )
            # rebalance tree
            if self.priority < self.cright.priority:
                croot = self.rotate_left()
        else:
            # insert to cleft tree
            if self.cleft is not EmptyNode:
                self.cleft = self.cleft.insert( start, end, interval)
            else:
                self.cleft = IntervalNode( start, end, interval)
            # rebalance tree
            if self.priority < self.cleft.priority:
                croot = self.rotate_right()
    
        croot.set_ends()
        self.cleft.croot  = croot
        self.cright.croot = croot
        return croot

    cdef IntervalNode rotate_right(IntervalNode self):
        cdef IntervalNode croot = self.cleft
        self.cleft  = self.cleft.cright
        croot.cright = self
        self.set_ends()
        return croot

    cdef IntervalNode rotate_left(IntervalNode self):
        cdef IntervalNode croot = self.cright
        self.cright = self.cright.cleft
        croot.cleft  = self
        self.set_ends()
        return croot

    cdef inline void set_ends(IntervalNode self):
        if self.cright is not EmptyNode and self.cleft is not EmptyNode: 
            self.maxend = imax3(self.end, self.cright.maxend, self.cleft.maxend)
            self.minend = imin3(self.end, self.cright.minend, self.cleft.minend)
            self.minstart = imin3(self.start, self.cright.minstart, self.cleft.minstart)
        elif self.cright is not EmptyNode:
            self.maxend = imax2(self.end, self.cright.maxend)
            self.minend = imin2(self.end, self.cright.minend)
            self.minstart = imin2(self.start, self.cright.minstart)
        elif self.cleft is not EmptyNode:
            self.maxend = imax2(self.end, self.cleft.maxend)
            self.minend = imin2(self.end, self.cleft.minend)
            self.minstart = imin2(self.start, self.cleft.minstart)
        

    def intersect( self, int start, int end, sort=True ):
        """
        given a start and a end, return a list of features
        falling within that range
        """
        cdef list results = []
        self._intersect( start, end, results )
        return results

    find = intersect
        
    cdef void _intersect( IntervalNode self, int start, int end, list results):
        # Left subtree
        if self.cleft is not EmptyNode and self.cleft.maxend > start:
            self.cleft._intersect( start, end, results )
        # This interval
        if ( self.end > start ) and ( self.start < end ):
            results.append( self.interval )
        # Right subtree
        if self.cright is not EmptyNode and self.start < end:
            self.cright._intersect( start, end, results )
    

    cdef void _seek_left(IntervalNode self, int position, list results, int n, int max_dist):
        # we know we can bail in these 2 cases.
        if self.maxend + max_dist < position:
            return
        if self.minstart > position:
            return

        # the ordering of these 3 blocks makes it so the results are
        # ordered nearest to farest from the query position
        if self.cright is not EmptyNode:
            self.cright._seek_left(position, results, n, max_dist)

        if -1 < position - self.end < max_dist:
            results.append(self.interval)

        # TODO: can these conditionals be more stringent?
        if self.cleft is not EmptyNode:
                self.cleft._seek_left(position, results, n, max_dist)


    
    cdef void _seek_right(IntervalNode self, int position, list results, int n, int max_dist):
        # we know we can bail in these 2 cases.
        if self.maxend < position: return
        if self.minstart - max_dist > position: return

        #print "SEEK_RIGHT:",self, self.cleft, self.maxend, self.minstart, position

        # the ordering of these 3 blocks makes it so the results are
        # ordered nearest to farest from the query position
        if self.cleft is not EmptyNode: 
                self.cleft._seek_right(position, results, n, max_dist)

        if -1 < self.start - position < max_dist:
            results.append(self.interval)

        if self.cright is not EmptyNode:
                self.cright._seek_right(position, results, n, max_dist)

    
    cpdef left(self, position, int n=1, int max_dist=2500):
        """
        find n features with a start > than `position`
        f: a Interval object (or anything with an `end` attribute)
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        cdef list results = []
        # use start - 1 becuase .left() assumes strictly left-of
        self._seek_left( position - 1, results, n, max_dist )
        if len(results) == n: return results
        r = results
        r.sort(key=operator.attrgetter('end'), reverse=True)
        return r[:n]

    cpdef right(self, position, int n=1, int max_dist=2500):
        """
        find n features with a end < than position
        f: a Interval object (or anything with a `start` attribute)
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        cdef list results = []
        # use end + 1 becuase .right() assumes strictly right-of
        self._seek_right(position + 1, results, n, max_dist)
        if len(results) == n: return results
        r = results
        r.sort(key=operator.attrgetter('start'))
        return r[:n]

    def traverse(self, func):
        self._traverse(func)

    cdef void _traverse(IntervalNode self, object func):
        if self.cleft is not EmptyNode: self.cleft._traverse(func)
        func(self)
        if self.cright is not EmptyNode: self.cright._traverse(func)

cdef IntervalNode EmptyNode = IntervalNode( 0, 0, Interval(0, 0))

cdef class Interval:
    """
    Basic feature, with required integer start and end properties.
    Also accepts optional value

    """

    def __init__(self, int start, int end, float value=0):
        assert start <= end, "start must be less than end"
        self.start  = start
        self.end   = end
        self.value = value

    def __repr__(self):
        fstr = "Interval(%d, %d" % (self.start, self.end)
        if not self.value is None:
            fstr += ", value=" + str(self.value)
        fstr += ")"
        return fstr

    def __richcmp__(self, other, op):
        if op == 0:
            # <
            return self.start < other.start or self.end < other.end
        elif op == 1:
            # <=
            return self == other or self < other
        elif op == 2:
            # ==
            return self.start == other.start and self.end == other.end and self.value == other.value
        elif op == 3:
            # !=
            return self.start != other.start or self.end != other.end or self.value != other.value
        elif op == 4:
            # >
            return self.start > other.start or self.end > other.end
        elif op == 5:
            # >=
            return self == other or self > other
    
    def merge(self,other,max_val=False):
        if max_val:
            new_val = max(self.value,other.value)
        else:
            new_val = min(self.value,other.value)
        return Interval(min(self.start,other.start),max(self.end,other.end),new_val)
    
    def substract(self, other):
        if self.start >= other.start:
            if self.end <= other.end:
                return []
            else:
                if self.start >= other.end:
                    return [self]
                else:
                    return [Interval(other.end,self.end,self.value)]
        else:
            if self.end <= other.start:
                return [self]
            else:
                if self.end <= other.end:
                    return [Interval(self.start,other.start,self.value)]
                else:
                    return [Interval(self.start,other.start,self.value),Interval(other.end,self.end,self.value)]

cdef class IntervalTree:
    """
    Data structure for performing window intersect queries on a set of 
    of possibly overlapping 1d intervals.

    """
    
    def __cinit__( self ):
        root = None
    
    def __reduce__( self ):
        return string2IntervalTree, (self.to_string(),)
    
    # ---- Position based interfaces -----------------------------------------
    
    def insert( self, int start, int end, object value=None ):
        """
        Insert the interval [start,end) associated with value `value`.
        """
        if self.root is None:
            self.root = IntervalNode( start, end, value )
        else:
            self.root = self.root.insert( start, end, value )
        
    def to_list( self):
        """
        Return a sorted list of all intervals.
        """
        if self.root is None:
            return []
        return self.root.find( self.root.minstart, self.root.maxend )

    def size( self):
        return len(self.to_list())

    def to_string( self):
        """
        Return a string of all intervals.
        """
        l = self.to_list()
        return string.join([str(x.start)+","+str(x.end) for x in l],";")
    
    @classmethod
    def from_string(cls, s):
        t = IntervalTree()
        l = [string.split(x,",") for x in string.split(s,sep=";")]
        l = [Interval(int(x[0]),int(x[1])) for x in l]
        for interval in l:
            t.insert_interval(interval)
        return t
    
    def find( self, start, end ):
        """
        Return a sorted list of all intervals overlapping [start,end).
        """
        if self.root is None:
            return []
        return self.root.find( start, end )
    
    def before( self, position, num_intervals=1, max_dist=2500 ):
        """
        Find `num_intervals` intervals that lie before `position` and are no
        further than `max_dist` positions away
        """
        if self.root is None:
            return []
        return self.root.left( position, num_intervals, max_dist )

    def after( self, position, num_intervals=1, max_dist=2500 ):
        """
        Find `num_intervals` intervals that lie after `position` and are no
        further than `max_dist` positions away
        """
        if self.root is None:
            return []
        return self.root.right( position, num_intervals, max_dist )

    # ---- Interval-like object based interfaces -----------------------------

    def insert_interval( self, interval ):
        """
        Insert an "interval" like object (one with at least start and end
        attributes)
        """
        self.insert( interval.start, interval.end, interval )

    def traverse(self, fn):
        """
        call fn for each element in the tree
        """
        if self.root is None:
            return None
        return self.root.traverse(fn)

    def intersect(self, other):
        if not isinstance(other, IntervalTree):
            raise ValueError
        intersections = []
        for interval in other.to_list():
            intersections += [(max(interval.start, i.start), min(interval.end, i.end)) for i in self.find(interval.start,interval.end)]
        return intersections