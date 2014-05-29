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
    
    cpdef TestSet get_slice(self, int start_snp, int snp_num):
        this_type = type(self)
        cdef TestSet other = this_type()
        other._nr_haplos = self._nr_haplos
        other._snp_num = min(snp_num, self._snp_num - start_snp)
        other._haplos = <bool **> malloc(self._nr_haplos * sizeof(bool *))
        cdef int hap_idx
        cdef int snp_idx
        for hap_idx in range(self._nr_haplos):
            other._haplos[hap_idx] = <bool *> malloc((other._snp_num) * sizeof(bool))
            #TODO: change this to pointer arithmetics instead of copy
            for snp_idx in range(other._snp_num):
                other._haplos[hap_idx][snp_idx] = self._haplos[hap_idx][snp_idx+start_snp]
        return other
    
    def __dealloc__(self):
        cdef int hap_idx
        for hap_idx in range(self._nr_haplos):
            free(self._haplos[hap_idx])
        free(self._haplos)
                         
    cpdef GenotypePair get_genotype_pair(self, int ind1, int ind2):
        '''
        ind1 and ind2 should be zero based
        '''
        this_type = type(self)
        cdef GenotypePair other = GenotypePair()
        other._snp_num = self._snp_num
        other._haplos = <bool **> malloc(other._nr_haplos * sizeof(bool *))
        cdef int hap_idx
        cdef int snp_idx
        for hap_idx in range(self._nr_haplos):
            other._haplos[hap_idx] = <bool *> malloc((self._snp_num) * sizeof(bool))
        for snp_idx in range(self._snp_num):
            other._haplos[0][snp_idx] = self._haplos[ind1*2][snp_idx]
            other._haplos[1][snp_idx] = self._haplos[ind1*2+1][snp_idx]
            other._haplos[2][snp_idx] = self._haplos[ind2*2][snp_idx]
            other._haplos[3][snp_idx] = self._haplos[ind2*2+1][snp_idx]
        return other
    
    def read_haplos(self, file_name, int max_snp_num=1000000000, ids = None, scramble=False):
        
        if not os.path.exists(file_name):
            print "the file: " + file_name + " does not exist!"
            exit(-1)
        
        count = 0
        with open(file_name, 'rU') as haplos_file:
            for line in haplos_file.xreadlines(  ):
                if count == 0: 
                    self._snp_num = len(line.rstrip('\r\n')) if len(line.rstrip('\r\n')) <= max_snp_num else max_snp_num
                if len(line.rstrip('\r\n')) < self._snp_num:
                    raise ValueError  
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
        
        with open(file_name, 'rU') as haplos_file:
            
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
                        line_trunc = line.rstrip('\r\n')[:self._snp_num] 
                        self._haplos[hap_idx] = <bool *> malloc((self._snp_num) * sizeof(bool))
                        for snp_idx in range(self._snp_num):
#                             print "snp: " + str(snp_idx)
                            if int(line_trunc[snp_idx]) == 0:
                                self._haplos[hap_idx][snp_idx] = 0
                            else:
                                if int(line_trunc[snp_idx]) == 1:
                                    self._haplos[hap_idx][snp_idx] = 1
                                else:
                                    print "unidentified allele: " + str(int(line_trunc[snp_idx]))
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
        
        print "finished reading haplos."
        return self._nr_haplos
    
    def write_haplos(self, file_name):
        with open(file_name, "w") as output_file:
            for hap_idx in range(self._nr_haplos):
                for snp_idx in range(self._snp_num):
                    output_file.write(str(int(self._haplos[hap_idx][snp_idx])))
                output_file.write("\n")
        
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
            
    cpdef set_ibd_segment(self, int start_snp, int snp_num):
        cdef int snp_idx
        for snp_idx in range(start_snp,start_snp+snp_num):
            self._haplos[2][snp_idx] = self._haplos[0][snp_idx]
            
    cpdef Genotype get_genotype(self, int ind):
        if ind != 0 and ind != 1:
            raise ValueError
        cdef Genotype other = Genotype()
        other._snp_num = self._snp_num
        other._haplos = <bool **> malloc(other._nr_haplos * sizeof(bool *))
        cdef int hap_idx
        cdef int snp_idx
        for hap_idx in range(self._nr_haplos):
            other._haplos[hap_idx] = <bool *> malloc((self._snp_num) * sizeof(bool))
        for snp_idx in range(self._snp_num):
            other._haplos[0][snp_idx] = self._haplos[ind*2][snp_idx]
            other._haplos[1][snp_idx] = self._haplos[ind*2+1][snp_idx]
        return other
            
cdef class Genotype(TestSet):
 
    def __cinit__(self):
        self._nr_haplos = 2
    
    cpdef generate_random_haps_inplace(self, LDModel m, int snp_num=1000000000):
        self._haplos = <bool **> malloc(self._nr_haplos * sizeof(bool *))
        self._snp_num = m._snp_num if m._snp_num <= snp_num else snp_num
        for hap_idx in range(2):
            self._haplos[hap_idx] = m.generate_random_hap(snp_num)  