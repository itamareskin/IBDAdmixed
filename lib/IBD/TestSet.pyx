#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

import random
import os
from libc.stdlib cimport malloc, free
from itertools import islice, combinations_with_replacement
from sys import stdout
 
cdef extern from "string.h":
    char *strncpy(char *dest, char *src, size_t n)
    int strlen(char *s) 

cdef class TestSet(object):

    def __cinit__(self):
        pass  
    
    def read_haplos(self, file_name, int max_snp_num=1000000000, ids = None, scramble=False):
    
        if not os.path.exists(file_name):
            print "the file: " + file_name + " does not exist!"
            exit(-1)
        
        count = 0
        with open(file_name) as haplos_file:
            for line in haplos_file.xreadlines(  ):
                if count == 0: 
                    self._snp_num = len(line) if len(line) <= max_snp_num else max_snp_num 
                count += 1
                
        if count == 0 or count % 2 == 1:
            print "bad number of haplotypes. quitting..."
            exit(-1)
        
        if ids == None:
            self._nr_haplos = count
            ids = range(self._nr_haplos)
        else:
            self._nr_haplos = len(ids)
        
        print "reading from haplos file: " + file_name
        self._haplos = <bool **> malloc(self._nr_haplos * sizeof(bool *)) 
        
        with open(file_name) as haplos_file:
            
            first_read = True
            done = False
            buffer_size = 100
            hap_idx = 0
            
            while True:
                
                if done:
                    break
                # read next buffer_size lines from the file
                lines = list(islice(haplos_file, buffer_size))
                
                if len(lines) == 0:
                    done = True
                
                for line in lines:
                    if hap_idx in ids:
                        line_trunc = line[:self._snp_num] 
                        self._haplos[hap_idx] = <bool *> malloc((self._snp_num) * sizeof(bool))
                        for snp_idx in range(self._snp_num):
                            if int(chr(line_trunc[snp_idx])) == 0:
                                self._haplos[hap_idx][snp_idx] = 0
                            else:
                                if int(chr(line_trunc[snp_idx])) == 1:
                                    self._haplos[hap_idx][snp_idx] = 1
                                else:
                                    print "unidentified allele: " + str(int(chr(line_trunc[snp_idx])))
                                    break
                        #self._haplos[hap_idx] = <char *> malloc((self._snp_num) * sizeof(char))
                        #strncpy(self._haplos[hap_idx], line_trunc,self._snp_num) 
                    hap_idx += 1
                    if hap_idx >= count:
                        done = True;
                        break
        
        if scramble:
            for hap_idx in range(0,self._nr_haplos,2):
                for snp_idx in range(self._snp_num):
                    switch = random.randint(0, 1)
                    if switch == 1:
                        temp_allele = self._haplos[hap_idx][snp_idx]
                        self._haplos[hap_idx][snp_idx] = self._haplos[hap_idx+1][snp_idx]
                        self._haplos[hap_idx+1][snp_idx] = temp_allele
        
        return self._nr_haplos
            
    cpdef generate_composite_individuals(self, LDModel m, num_inds):
        
        self._nr_haplos = 2 * num_inds
        new_haplos = <bool **> malloc(self._nr_haplos * sizeof(bool *))
        for hap_idx in range(self._nr_haplos):
            new_haplos[hap_idx] = <bool *> malloc(m._snp_num * sizeof(bool))
    
        for i in range(self._nr_haplos):
            print "composing individual " + str(i)
            print "first position: " + str(m._gm._position[0])
            print "last position: " + str(m._gm._position[self._snp_num-1])
            start = 0
            end = 0
            length = int(200 * (self._position[m._snp_num-1] - m._gm._position[0]) / m._snp_num) 
            same_j = 0
            last_j = -1
            while end < m._snp_num - 1:
                stdout.write("%d," % start)
                stdout.flush()
                start = end 
                end = min(m._snp_num-1,start+length) 
                #dists_c = [abs(x-self._genetic_dist[start]-length) for x in self._genetic_dist]
                #end = dists_c.index(min(dists_c))
                j = random.randint(0, self._nr_haplos)
                if last_j == j:
                    same_j+=1
                    if same_j >= 8:
                        while True:
                            j = random.randint(0, self._nr_haplos)
                            if last_j != j:
                                same_j = 0
                                break
                else:
                    same_j = 0 
                    
                last_j = j
                #strncpy(self._haplos[hap_idx], line_trunc,self._snp_num)
                #new_haplos[i][start:end] = self._haplos[j][start:end]
        
cdef class GenotypePair(TestSet):
 
    def __cinit__(self):
        self._nr_haplos = 4
    
    cpdef generate_random_haps_inplace(self, LDModel m, int snp_num=1000000000):
        self._haplos = <bool **> malloc(self._nr_haplos * sizeof(bool *))
        self._snp_num = m._snp_num if m._snp_num <= snp_num else snp_num
        for hap_idx in range(4):
            self._haplos[hap_idx] = m.generate_random_hap(snp_num)
            
cdef class Genotype(TestSet):
 
    def __cinit__(self):
        self._nr_haplos = 2
    
    cpdef generate_random_haps_inplace(self, LDModel m, int snp_num=1000000000):
        self._haplos = <bool **> malloc(self._nr_haplos * sizeof(bool *))
        self._snp_num = m._snp_num if m._snp_num <= snp_num else snp_num
        for hap_idx in range(2):
            self._haplos[hap_idx] = m.generate_random_hap(snp_num)  