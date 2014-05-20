#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

from libc.stdlib cimport malloc, free
import os
from itertools import islice, combinations_with_replacement

cdef class GeneticMap(object):
    '''
    represents the genetic map data for a chromosome 
    '''
    
    def __cinit__(self, map_file_name, int max_snp_num=1000000000):
        
        # total number of SNPs to be analyzed
        with open(map_file_name) as map_file:
            data = map_file.readlines()
            self._snp_num = len(data) if len(data) <= max_snp_num else max_snp_num
        
        self._position = <int *> malloc(self._snp_num * sizeof(int))
        self._genetic_dist = <double *> malloc(self._snp_num * sizeof(double))
        
        self.read_map(map_file_name)
        
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
                    self._position[snp_idx] = int(line[3])
                    snp_idx += 1
                    if snp_idx >= self._snp_num:
                        done = True;
                        break