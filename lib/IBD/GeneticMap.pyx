#cython: profile=True
#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

from libc.stdlib cimport malloc, free
from libc.stdlib cimport bsearch
import os
from itertools import islice, combinations_with_replacement

#cdef int _qsortcmp(void const *v, void const *u) nogil:
#    return (v[0] > u[0]) - (v[0] < u[0])

cdef class GeneticMap(object):
    '''
    represents the genetic map data for a chromosome 
    '''
    
    def __cinit__(self, map_file_name=None, int max_snp_num=1000000000):
        
        # total number of SNPs to be analyzed
        if map_file_name is not None:
            with open(map_file_name) as map_file:
                data = map_file.readlines()
                self._snp_num = len(data) if len(data) <= max_snp_num else max_snp_num
            
            self._position = <long *> malloc(self._snp_num * sizeof(long))
            self._genetic_dist = <double *> malloc(self._snp_num * sizeof(double))
            
            self.read_map(map_file_name)
        
            print "genetic map created"

        self._is_slice = False
            
    def __dealloc__(self):
        if not self._is_slice:
            free(self._position)
            free(self._genetic_dist)
    
    cpdef GeneticMap get_slice(self, int start_snp, int snp_num):
        cdef GeneticMap other = GeneticMap()
        other._snp_num = min(snp_num, self._snp_num - start_snp)
        other._position = self._position + start_snp
        other._genetic_dist = self._genetic_dist + start_snp
        other._is_slice = True
        return other
        
    def read_map(self, map_file_name):
        '''
        read the (PLINK) map file
        '''
        if not os.path.exists(map_file_name):
            print "the file: " + map_file_name + " does not exist!"
            exit(-1)
        
        print "reading from map file: " + map_file_name
        
        with open(map_file_name) as map_file:
            
            done = False
            buffer_size = 100
            snp_idx = 0
            
            while True:
                
                if done:
                    break
                # read next buffer_size lines from the file
                lines = list(islice(map_file, buffer_size))
                
                if len(lines) == 0:
                    done = True
                
                for line in lines:
                    #print "line: " + line
                    line = line.split(" ")
                    #print "pos: " + line[0]
                    self._genetic_dist[snp_idx] = float(line[2])
                    self._position[snp_idx] = long(line[3])
                    snp_idx += 1
                    if snp_idx >= self._snp_num:
                        done = True;
                        break
    
    cpdef dict get_position_dict(self):
        cdef snp_idx
        cdef dict pos_dict = {}
        for snp_idx in range(self._snp_num):
            pos_dict[self._position[snp_idx]] = snp_idx
        return pos_dict

    cpdef list get_genetic_dist_list(self):
        cdef snp_idx
        cdef list dists = []
        for snp_idx in range(self._snp_num):
            dists.append(self._genetic_dist[snp_idx])
        return dists


    cpdef float get_length(self, int start_snp_idx, int end_snp_idx):
        return self._genetic_dist[min(end_snp_idx,self._snp_num-1)] - self._genetic_dist[start_snp_idx]