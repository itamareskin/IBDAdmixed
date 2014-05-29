#cython: profile=True
#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

from __future__ import division
from InnerModel import InnerModel
from InnerModel cimport InnerModel
from TestSet import TestSet 
from TestSet cimport TestSet
from TestSet import Genotype
from TestSet cimport Genotype
from LDModel import LDModel
from LDModel cimport LDModel
import os
import math
from libc.math cimport exp, log
from libc.float cimport DBL_MIN, DBL_MAX
from libc.limits cimport ULONG_MAX,LONG_MIN
from libc.stdlib cimport malloc, free
from itertools import islice, combinations_with_replacement
from libcpp cimport bool
import logging
#inner_forward_probs_logger = logging.getLogger('innerforwardprobs')

cdef class GenotypeModel(InnerModel):
    
    def __cinit__(self, LDModel m1, LDModel m2, bool phased=False, int g=8):        
        self._phased = phased
        self._g = g
        self._m1 = m1
        self._m2 = m2
    
    cpdef InnerModel slice_from_model(self, int start_snp, int snp_num):
        return GenotypeModel(self._m1.get_slice_model(start_snp,snp_num),
                                 self._m2.get_slice_model(start_snp,snp_num),
                                 self._phased,
                                 self._g)
    
    cpdef alloc_mem(self):
        cdef int snp_idx
        cdef int node_idx1
        cdef int node_idx2
        
        self._forward_prob = <double ***> malloc((self._m1._snp_num) * sizeof(double**))
        self._backward_prob = <double ***> malloc((self._m1._snp_num) * sizeof(double**))
        self._emission_prob = <double ***> malloc((self._m1._snp_num) * sizeof(double**))
        self._scale_factor = <double *> malloc((self._m1._snp_num) * sizeof(double))
        self._backward_scale_factor = <double *> malloc((self._m1._snp_num) * sizeof(double))
        for snp_idx in range(self._m1._snp_num):
            self._forward_prob[snp_idx]= <double **> malloc((self._m1._layer_state_nums[snp_idx]) * sizeof(double*))
            self._backward_prob[snp_idx]= <double **> malloc((self._m1._layer_state_nums[snp_idx]) * sizeof(double*))
            self._emission_prob[snp_idx]= <double **> malloc((self._m1._layer_state_nums[snp_idx]) * sizeof(double*))
            self._scale_factor[snp_idx] = 0
            self._backward_scale_factor[snp_idx] = 0
            for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
                self._forward_prob[snp_idx][node_idx1] = <double *> malloc((self._m2._layer_state_nums[snp_idx]) * sizeof(double))
                self._backward_prob[snp_idx][node_idx1] = <double *> malloc((self._m2._layer_state_nums[snp_idx]) * sizeof(double))
                self._emission_prob[snp_idx][node_idx1] = <double *> malloc((self._m2._layer_state_nums[snp_idx]) * sizeof(double))
                for node_idx2 in range(self._m2._layer_state_nums[snp_idx]):
                    self._forward_prob[snp_idx][node_idx1][node_idx2] = 0
                    self._backward_prob[snp_idx][node_idx1][node_idx2] = 0
                    self._emission_prob[snp_idx][node_idx1][node_idx2] = 0
    
    cpdef free_mem(self):
        cdef int snp_idx
        cdef int node_idx1
        cdef int node_idx2
        
        free(self._scale_factor)
        free(self._backward_scale_factor)
        for snp_idx in range(self._m1._snp_num):
            for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
                    free(self._forward_prob[snp_idx][node_idx1])
                    free(self._backward_prob[snp_idx][node_idx1])
                    free(self._emission_prob[snp_idx][node_idx1])
        free(self._forward_prob)
        free(self._backward_prob)
        free(self._emission_prob)
        
    cpdef calc_emission_probs(self, Genotype p):
        
        cdef int snp_idx 
        cdef int node_idx1
        cdef int node_idx2

        for snp_idx in range(self._m1._snp_num):
            print "calculating emission probs, snp: " + str(snp_idx)
            for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
                for node_idx2 in range(self._m2._layer_state_nums[snp_idx]):
                    print "nodes: "+ str(node_idx1) + " " + str(node_idx2)
                    if not self._phased:
                        print "1"
                        if p.chr1(snp_idx) == p.chr2(snp_idx):
                            print "2"
                            self._emission_prob[snp_idx][node_idx1][node_idx2] = \
                            self._m1._states[snp_idx][node_idx1].prob_em[p.chr1(snp_idx)] * self._m2._states[snp_idx][node_idx2].prob_em[p.chr2(snp_idx)]
                        else:
                            print "3"
                            self._emission_prob[snp_idx][node_idx1][node_idx2] = \
                            self._m1._states[snp_idx][node_idx1].prob_em[p.chr1(snp_idx)] * self._m2._states[snp_idx][node_idx2].prob_em[p.chr2(snp_idx)] + \
                            self._m1._states[snp_idx][node_idx1].prob_em[p.chr2(snp_idx)] * self._m2._states[snp_idx][node_idx2].prob_em[p.chr1(snp_idx)]
                    else:
                        print "4"
                        self._emission_prob[snp_idx][node_idx1][node_idx2] = \
                        self._m1._states[snp_idx][node_idx1].prob_em[p.chr1(snp_idx)] * \
                        self._m2._states[snp_idx][node_idx2].prob_em[p.chr2(snp_idx)]
                                                       
    cpdef calc_forward_probs(self):
        
        cdef int snp_idx 
        cdef int node_idx1
        cdef int node_idx2
        cdef int prev_node_idx1
        cdef int prev_node_idx2
        cdef int prev_node1
        cdef int prev_node2
        cdef double eps_or_1_eps

        # first layer
        for node_idx1 in range(self._m1._layer_state_nums[0]):
            for node_idx2 in range(self._m2._layer_state_nums[0]):
                self._forward_prob[0][node_idx1][node_idx2] = \
                self._m1._pi[node_idx1] * self._m2._pi[node_idx2] * \
                self._emission_prob[0][node_idx1][node_idx2]
        # rescaling to avoid underflow
        self._scale_factor[0] = 1
    
        # all other layers
        for snp_idx in range(self._m1._snp_num-1):
            # calculate forward probabilities
            for node_idx1 in range(self._m1._layer_state_nums[snp_idx+1]):
                for node_idx2 in range(self._m2._layer_state_nums[snp_idx+1]):
                    for prev_node_idx1 in range(self._m1._states[snp_idx+1][node_idx1].in_trans_num):
                        for prev_node_idx2 in range(self._m2._states[snp_idx+1][node_idx2].in_trans_num):
                            prev_node1 = self._m1._back_trans_idx[snp_idx+1][node_idx1][prev_node_idx1]
                            prev_node2 = self._m2._back_trans_idx[snp_idx+1][node_idx2][prev_node_idx2]
                            self._forward_prob[snp_idx+1][node_idx1][node_idx2] += \
                            self._forward_prob[snp_idx][prev_node1][prev_node2] * \
                            self._emission_prob[snp_idx+1][node_idx1][node_idx2] * \
                            self._m1._back_trans[snp_idx+1][node_idx1][prev_node_idx1] * \
                            self._m2._back_trans[snp_idx+1][node_idx2][prev_node_idx2]
            
            # rescaling to avoid underflow
            self.rescale_forward(snp_idx+1)

    cdef rescale_forward(self, int snp_idx):
        cdef int node_idx1
        cdef int node_idx2
        
        for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
            for node_idx2 in range(self._m2._layer_state_nums[snp_idx]):
                self._scale_factor[snp_idx] += self._forward_prob[snp_idx][node_idx1][node_idx2]
        if self._scale_factor[snp_idx] > 0:                 
            self._scale_factor[snp_idx] = 1.0 / self._scale_factor[snp_idx]
        else: 
            self._scale_factor[snp_idx] = DBL_MAX
                                
        for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
            for node_idx2 in range(self._m2._layer_state_nums[snp_idx]):
                self._forward_prob[snp_idx][node_idx1][node_idx2] = \
                self._forward_prob[snp_idx][node_idx1][node_idx2] * self._scale_factor[snp_idx]
    
    cpdef calc_backward_probs(self):
         
        cdef int node_idx1
        cdef int node_idx2
        cdef int nxt_node_idx1
        cdef int nxt_node_idx2
        cdef int nxt_node1
        cdef int nxt_node2
        
        # last layer
        for node_idx1 in range(self._m1._layer_state_nums[self._m1._snp_num - 1]):
            for node_idx2 in range(self._m2._layer_state_nums[self._m1._snp_num - 1]):
                self._backward_prob[self._m1._snp_num - 1][node_idx1][node_idx2] = 1
                                                
        # rescaling to avoid underflow
        for node_idx1 in range(self._m1._layer_state_nums[self._m1._snp_num - 1]):
            for node_idx2 in range(self._m2._layer_state_nums[self._m1._snp_num - 1]):
                self._backward_scale_factor[self._m1._snp_num - 1] += self._backward_prob[self._m1._snp_num - 1][node_idx1][node_idx2]
        self._backward_scale_factor[self._m1._snp_num - 1] = 1.0
        
        # all other layers
        snp_idx = self._m1._snp_num - 2 
        for snp_idx in reversed(range(0, self._m1._snp_num - 1)):
            # calculate forward probabilities
            for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
                for node_idx2 in range(self._m2._layer_state_num[snp_idx]):
                    for nxt_node_idx1 in range(self._m1._states[snp_idx][node_idx1].out_trans_num):
                        for nxt_node_idx2 in range(self._m2._states[snp_idx][node_idx2].out_trans_num):
                                    #for nxt_ibd in range(2):
                            nxt_node1 = self._m1._trans_idx[snp_idx][node_idx1][nxt_node_idx1]
                            nxt_node2 = self._m2._trans_idx[snp_idx][node_idx2][nxt_node_idx2]
                            self._backward_prob[snp_idx][node_idx1][node_idx2] += \
                            self._backward_prob[snp_idx+1][nxt_node1][nxt_node2] * \
                            self._m1._trans[snp_idx][node_idx1][nxt_node_idx1] * \
                            self._m2._trans[snp_idx][node_idx2][nxt_node_idx2] * \
                            self._emission_prob[snp_idx+1][nxt_node1][nxt_node2] #* \
            # rescaling to avoid underflow
            self.rescale_backward(snp_idx)

    cdef rescale_backward(self, int snp_idx):
        for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
            for node_idx2 in range(self._m2._layer_state_num[snp_idx]):
                if snp_idx > 0:
                    self._backward_scale_factor[snp_idx] += self._backward_prob[snp_idx][node_idx1][node_idx2]
                else:
                    self._backward_scale_factor[snp_idx] += self._backward_prob[snp_idx][node_idx1][node_idx2] * \
                    self._emission_prob[0][node_idx1][node_idx2] * \
                    self._m1._pi[node_idx1] * self._m2._pi[node_idx2]
                                        
                    if self._backward_scale_factor[snp_idx] > 0:                 
                        self._backward_scale_factor[snp_idx] = 1.0 / self._backward_scale_factor[snp_idx]
                    else: 
                        self._backward_scale_factor[snp_idx] = DBL_MAX
                            
        for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
            for node_idx2 in range(self._m2._layer_state_nums[snp_idx]):
                self._backward_prob[snp_idx][node_idx1][node_idx2] = \
                self._backward_prob[snp_idx][node_idx1][node_idx2] * self._backward_scale_factor[snp_idx]
    
    cpdef double calc_likelihood(self, TestSet obs_data):
        cdef Genotype g = <Genotype?>obs_data
        cdef int snp_idx       
        cdef likelihood = 0
        
        self.alloc_mem()
        self.calc_emission_probs(g)
        self.calc_forward_probs()
        #self.print_inner_prob() 
        for snp_idx in range(self._m1._snp_num):
            likelihood = likelihood - log(self._scale_factor[snp_idx])
        self.free_mem()
        return likelihood
    
    cdef double ibd_trans_prob(self, GenotypeModel other):
        cdef double d = self._m1._gm._genetic_dist[self._m1._snp_num - 1] - self._m1._gm._genetic_dist[0]
        cdef double s_0_0, s_0_1, s_1_0, s_1_1
        s_1_1 = exp(-self._m1._t_1_0 * d)
        s_0_0 = exp(-self._m1._t_0_1 * d)
        s_1_0 = (1 - s_1_1)
        s_0_1 = (1 - s_0_0)
        if self._ibd == 1 and other._ibd == 1:
            return s_1_1
        if self._ibd == 0 and other._ibd == 0:
            return s_0_0
        if self._ibd == 1 and other._ibd == 0:
            return s_1_0
        if self._ibd == 0 and other._ibd == 1:
            return s_0_1
    
    cdef double anc_trans_prob(self, GenotypeModel other):
        cdef double win_recomb = (self._m1._gm._genetic_dist[self._m1._snp_num - 1] - self._m1._gm._genetic_dist[0])
        if win_recomb > 0:
            return ((self._g - 1) * win_recomb)
        else:
            return 1
    
    cpdef double trans_prob(self, InnerModel other):
        cdef GenotypeModel other_gpm = <GenotypeModel?>other
        return self.anc_trans_prob(other_gpm) * self.ibd_trans_prob(other_gpm)
    
    cpdef print_inner_prob(self):
        for snp_idx in range(self._m1._snp_num):
            for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
                for node_idx2 in range(self._m2._layer_state_nums[snp_idx]):
                    print str(snp_idx) + " " + \
                    str(node_idx1) + " " + \
                    str(node_idx2) + " " + \
                    str(self._forward_prob[snp_idx][node_idx1][node_idx2]) + " " + \
                    str(self._backward_prob[snp_idx][node_idx1][node_idx2]) + " " + \
                    str(self._emission_prob[snp_idx][node_idx1][node_idx2])
            
