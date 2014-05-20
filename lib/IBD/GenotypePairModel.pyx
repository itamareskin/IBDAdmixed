#cython: profile=True
#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

from __future__ import division
from InnerModel import InnerModel
from InnerModel cimport InnerModel
from TestSet import GenotypePair
from TestSet cimport GenotypePair
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

cdef class GenotypePairModel(InnerModel):
    
    def __cinit__(self, LDModel m1, LDModel m2, LDModel m3, LDModel m4, bool phased, int ibd):
        cdef int snp_idx
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        
        self._phased = phased
        self._ibd = ibd
        self._m1 = m1
        self._m2 = m2
        self._m3 = m3
        self._m4 = m4
         
        self._forward_prob = <double *****> malloc((self._m1._snp_num) * sizeof(double****))
        self._backward_prob = <double *****> malloc((self._m1._snp_num) * sizeof(double****))
        self._emission_prob = <double *****> malloc((self._m1._snp_num) * sizeof(double****))
        self._scale_factor = <double *> malloc((self._m1._snp_num) * sizeof(double))
        self._backward_scale_factor = <double *> malloc((self._m1._snp_num) * sizeof(double))
        for snp_idx in range(self._m1._snp_num):
            self._forward_prob[snp_idx]= <double ****> malloc((self._m1._layer_state_nums[snp_idx]) * sizeof(double***))
            self._backward_prob[snp_idx]= <double ****> malloc((self._m1._layer_state_nums[snp_idx]) * sizeof(double***))
            self._emission_prob[snp_idx]= <double ****> malloc((self._m1._layer_state_nums[snp_idx]) * sizeof(double***))
            self._scale_factor[snp_idx] = 0
            self._backward_scale_factor[snp_idx] = 0
            for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
                self._forward_prob[snp_idx][node_idx1] = <double ***> malloc((self._m2._layer_state_nums[snp_idx]) * sizeof(double**))
                self._backward_prob[snp_idx][node_idx1] = <double ***> malloc((self._m2._layer_state_nums[snp_idx]) * sizeof(double**))
                self._emission_prob[snp_idx][node_idx1] = <double ***> malloc((self._m2._layer_state_nums[snp_idx]) * sizeof(double**))
                for node_idx2 in range(self._m2._layer_state_nums[snp_idx]):
                    self._forward_prob[snp_idx][node_idx1][node_idx2] = <double **> malloc((self._m3._layer_state_nums[snp_idx]) * sizeof(double*))
                    self._backward_prob[snp_idx][node_idx1][node_idx2] = <double **> malloc((self._m3._layer_state_nums[snp_idx]) * sizeof(double*))
                    self._emission_prob[snp_idx][node_idx1][node_idx2] = <double **> malloc((self._m3._layer_state_nums[snp_idx]) * sizeof(double*))
                    for node_idx3 in range(self._m3._layer_state_nums[snp_idx]):
                        self._forward_prob[snp_idx][node_idx1][node_idx2][node_idx3] = <double *> malloc((self._m4._layer_state_nums[snp_idx]) * sizeof(double))
                        self._backward_prob[snp_idx][node_idx1][node_idx2][node_idx3] = <double *> malloc((self._m4._layer_state_nums[snp_idx]) * sizeof(double))
                        self._emission_prob[snp_idx][node_idx1][node_idx2][node_idx3] = <double *> malloc((self._m4._layer_state_nums[snp_idx]) * sizeof(double))
                        for node_idx4 in range(self._m4._layer_state_nums[snp_idx]):
                            self._forward_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4] = 0                        
                            self._backward_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4] = 0                        
                            self._emission_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4] = 0                        

    cpdef calc_emission_probs(self, GenotypePair p):
        
        cdef int snp_idx 
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4

        for snp_idx in range(self._m1._snp_num):
            for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
                for node_idx2 in range(self._m2._layer_state_nums[snp_idx]):
                    for node_idx3 in range(self._m3._layer_state_nums[snp_idx]):
                        for node_idx4 in range(self._m4._layer_state_nums[snp_idx]):
                            if not self._phased:
                                if p.chr1(snp_idx) == p.chr2(snp_idx):
                                    self._emission_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4] = \
                                    self._m1._states[snp_idx][node_idx1].prob_em[p.chr1(snp_idx)] * self._m2._states[snp_idx][node_idx2].prob_em[p.chr2(snp_idx)]
                                else:
                                    self._emission_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4] = \
                                    self._m1._states[snp_idx][node_idx1].prob_em[p.chr1(snp_idx)] * self._m2._states[snp_idx][node_idx2].prob_em[p.chr2(snp_idx)] + \
                                    self._m1._states[snp_idx][node_idx1].prob_em[p.chr2(snp_idx)] * self._m2._states[snp_idx][node_idx2].prob_em[p.chr1(snp_idx)]
                                
                                if p.chr3(snp_idx) == p.chr4(snp_idx):
                                    self._emission_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4] *= \
                                    self._m3._states[snp_idx][node_idx3].prob_em[p.chr3(snp_idx)] * self._m4._states[snp_idx][node_idx4].prob_em[p.chr4(snp_idx)]
                                else:
                                    self._emission_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4] *= \
                                    self._m3._states[snp_idx][node_idx3].prob_em[p.chr3(snp_idx)] * self._m4._states[snp_idx][node_idx4].prob_em[p.chr4(snp_idx)] + \
                                    self._m3._states[snp_idx][node_idx3].prob_em[p.chr4(snp_idx)] * self._m4._states[snp_idx][node_idx4].prob_em[p.chr3(snp_idx)]
                            else:
                                self._emission_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4] = \
                                self._m1._states[snp_idx][node_idx1].prob_em[p.chr1(snp_idx)] * \
                                self._m2._states[snp_idx][node_idx2].prob_em[p.chr2(snp_idx)] * \
                                self._m3._states[snp_idx][node_idx3].prob_em[p.chr3(snp_idx)] * \
                                self._m4._states[snp_idx][node_idx4].prob_em[p.chr4(snp_idx)]
                                                       
    cpdef calc_forward_probs(self):
        
        cdef int snp_idx 
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        cdef int prev_node_idx1
        cdef int prev_node_idx2
        cdef int prev_node_idx3
        cdef int prev_node_idx4
        cdef int prev_node1
        cdef int prev_node2
        cdef int prev_node3
        cdef int prev_node4
        cdef double eps_or_1_eps

        # first layer
        for node_idx1 in range(self._m1._layer_state_nums[0]):
            for node_idx2 in range(self._m2._layer_state_nums[0]):
                for node_idx3 in range(self._m3._layer_state_nums[0]):
                    for node_idx4 in range(self._m4._layer_state_nums[0]):
                        if self._ibd == 0:
                            self._forward_prob[0][node_idx1][node_idx2][node_idx3][node_idx4] = \
                            self._m1._pi[node_idx1] * self._m2._pi[node_idx2] * self._m3._pi[node_idx3] * self._m4._pi[node_idx4] * \
                            self._emission_prob[0][node_idx1][node_idx2][node_idx3][node_idx4]
                        else:
                            if self._m1._anc == self._mself._m1.anc:
                                self._forward_prob[0][node_idx1][node_idx2][node_idx3][node_idx4] = \
                                self._m1._pi[node_idx1] * self._m2._pi[node_idx2] * self._m4._pi[node_idx4] * \
                                self._emission_prob[0][node_idx1][node_idx2][node_idx3][node_idx4]
        # rescaling to avoid underflow
        self._scale_factor[0] = 1
    
        # all other layers
        for snp_idx in range(self._m1._snp_num):
            # calculate forward probabilities
            for node_idx1 in range(self._m1._layer_state_nums[snp_idx+1]):
                for node_idx2 in range(self._m2._layer_state_nums[snp_idx+1]):
                    for node_idx3 in range(self._m3._layer_state_nums[snp_idx+1]):
                        for node_idx4 in range(self._m4._layer_state_nums[snp_idx+1]):
                            for prev_node_idx1 in range(self._m1._states[snp_idx+1][node_idx1].in_trans_num):
                                for prev_node_idx2 in range(self._m2._states[snp_idx+1][node_idx2].in_trans_num):
                                    for prev_node_idx3 in range(self._m3._states[snp_idx+1][node_idx3].in_trans_num):
                                        for prev_node_idx4 in range(self._m4._states[snp_idx+1][node_idx4].in_trans_num):
                                            #for prev_ibd in range(2):
                                                prev_node1 = self._m1._back_trans_idx[snp_idx+1][node_idx1][prev_node_idx1]
                                                prev_node2 = self._m2._back_trans_idx[snp_idx+1][node_idx2][prev_node_idx2]
                                                prev_node3 = self._m3._back_trans_idx[snp_idx+1][node_idx3][prev_node_idx3]
                                                prev_node4 = self._m4._back_trans_idx[snp_idx+1][node_idx4][prev_node_idx4]
                                                if self._ibd == 0:
                                                    self._forward_prob[snp_idx+1][node_idx1][node_idx2][node_idx3][node_idx4] += \
                                                    self._forward_prob[snp_idx][prev_node1][prev_node2][prev_node3][prev_node4] * \
                                                    self._emission_prob[snp_idx+1][node_idx1][node_idx2][node_idx3][node_idx4] * \
                                                    self._m1._back_trans[snp_idx+1][node_idx1][prev_node_idx1] * \
                                                    self._m2._back_trans[snp_idx+1][node_idx2][prev_node_idx2] * \
                                                    self._m3._back_trans[snp_idx+1][node_idx3][prev_node_idx3] * \
                                                    self._m4._back_trans[snp_idx+1][node_idx4][prev_node_idx4]
                                                else:
                                                    if self._m1._anc == self._m3._anc and node_idx1 == node_idx3 and prev_node1 == prev_node3:
                                                        self._forward_prob[snp_idx+1][node_idx1][node_idx2][node_idx3][node_idx4] += \
                                                        self._forward_prob[snp_idx][prev_node1][prev_node2][prev_node3][prev_node4] * \
                                                        self._emission_prob[snp_idx+1][node_idx1][node_idx2][node_idx3][node_idx4] * \
                                                        self._m1._back_trans[snp_idx+1][node_idx1][prev_node_idx1] * \
                                                        self._m2._back_trans[snp_idx+1][node_idx2][prev_node_idx2] * \
                                                        self._m4._back_trans[snp_idx+1][node_idx4][prev_node_idx4]
            
            # rescaling to avoid underflow
            self.rescale_forward(snp_idx+1)

    cdef rescale_forward(self, int snp_idx):
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        
        for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
            for node_idx2 in range(self._m2._layer_state_nums[snp_idx]):
                for node_idx3 in range(self._m3._layer_state_nums[snp_idx]):
                    for node_idx4 in range(self._m4._layer_state_nums[snp_idx]):
                        self._scale_factor[snp_idx] += self._forward_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4]
        if self._scale_factor[snp_idx] > 0:                 
            self._scale_factor[snp_idx] = 1.0 / self._scale_factor[snp_idx]
        else: 
            self._scale_factor[snp_idx+1] = DBL_MAX
                                
        for node_idx1 in range(self._m1._layer_state_nums[snp_idx+1]):
            for node_idx2 in range(self._m2._layer_state_nums[snp_idx+1]):
                for node_idx3 in range(self._m3._layer_state_nums[snp_idx+1]):
                    for node_idx4 in range(self._m4._layer_state_nums[snp_idx+1]):
                        self._forward_prob[snp_idx+1][node_idx1][node_idx2][node_idx3][node_idx4] = \
                        self._forward_prob[snp_idx+1][node_idx1][node_idx2][node_idx3][node_idx4] * self._scale_factor[snp_idx+1]
    
    cpdef calc_backward_probs(self):
         
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        cdef int nxt_node_idx1
        cdef int nxt_node_idx2
        cdef int nxt_node_idx3
        cdef int nxt_node_idx4
        cdef int nxt_node1
        cdef int nxt_node2
        cdef int nxt_node3
        cdef int nxt_node4
        
        # last layer
        for node_idx1 in range(self._m1._layer_state_nums[self._m1._snp_num - 1]):
            for node_idx2 in range(self._m2._layer_state_nums[self._m1._snp_num - 1]):
                for node_idx3 in range(self._m3._layer_state_nums[self._m1._snp_num - 1]):
                    for node_idx4 in range(self._m4._layer_state_nums[self._m1._snp_num - 1]):
                        if self._ibd == 0 or self._m1._anc == self._m3._anc:
                            self._backward_prob[self._m1._snp_num - 1][node_idx1][node_idx2][node_idx3][node_idx4] = 1
                                                
        # rescaling to avoid underflow
        for node_idx1 in range(self._m1._layer_state_nums[self._m1._snp_num - 1]):
            for node_idx2 in range(self._m2._layer_state_nums[self._m1._snp_num - 1]):
                for node_idx3 in range(self._m3._layer_state_nums[self._m1._snp_num - 1]):
                    for node_idx4 in range(self._m4._layer_state_nums[self._m1._snp_num - 1]):
                        self._backward_scale_factor[self._m1._snp_num - 1] += self._backward_prob[self._m1._snp_num - 1][node_idx1][node_idx2][node_idx3][node_idx4]
        self._backward_scale_factor[self._m1._snp_num - 1] = 1.0
        
        # all other layers
        snp_idx = self._m1._snp_num - 2 
        for snp_idx in reversed(range(0, self._m1._snp_num - 1)):
            # calculate forward probabilities
            for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
                for node_idx2 in range(self._m2._layer_state_num[snp_idx]):
                    for node_idx3 in range(self._m3._layer_state_nums[snp_idx]):
                        for node_idx4 in range(self._m4._layer_state_nums[snp_idx]):
                            for nxt_node_idx1 in range(self._m1._states[snp_idx][node_idx1].out_trans_num):
                                for nxt_node_idx2 in range(self._m2._states[snp_idx][node_idx2].out_trans_num):
                                    for nxt_node_idx3 in range(self._m3._states[snp_idx][node_idx3].out_trans_num):
                                        for nxt_node_idx4 in range(self._m4._states[snp_idx][node_idx4].out_trans_num):
                                            #for nxt_ibd in range(2):
                                                nxt_node1 = self._m1._trans_idx[snp_idx][node_idx1][nxt_node_idx1]
                                                nxt_node2 = self._m2._trans_idx[snp_idx][node_idx2][nxt_node_idx2]
                                                nxt_node3 = self._m3._trans_idx[snp_idx][node_idx3][nxt_node_idx3]
                                                nxt_node4 = self._m4._trans_idx[snp_idx][node_idx4][nxt_node_idx4]
                                                if self._ibd == 0:
                                                    self._backward_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4] += \
                                                    self._backward_prob[snp_idx+1][nxt_node1][nxt_node2][nxt_node3][nxt_node4] * \
                                                    self._m1._trans[snp_idx][node_idx1][nxt_node_idx1] * \
                                                    self._m2._trans[snp_idx][node_idx2][nxt_node_idx2] * \
                                                    self._m3._trans[snp_idx][node_idx3][nxt_node_idx3] * \
                                                    self._m4._trans[snp_idx][node_idx4][nxt_node_idx4] * \
                                                    self._emission_prob[snp_idx+1][nxt_node1][nxt_node2][nxt_node3][nxt_node4] #* \
                                                else:
                                                    if self._m1._anc == self._m3._anc and node_idx1 == node_idx3 and nxt_node1 == nxt_node3:
                                                        self._backward_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4] += \
                                                        self._backward_prob[snp_idx+1][nxt_node1][nxt_node2][nxt_node3][nxt_node4] * \
                                                        self._m1._trans[snp_idx][node_idx1][nxt_node_idx1] * \
                                                        self._m2._trans[snp_idx][node_idx2][nxt_node_idx2] * \
                                                        self._m3._trans[snp_idx][node_idx4][nxt_node_idx4] * \
                                                        self._m4._emission_prob[snp_idx+1][nxt_node1][nxt_node2][nxt_node3][nxt_node4]
            # rescaling to avoid underflow
            self.rescale_backward(snp_idx)

    cdef rescale_backward(self, int snp_idx):
        for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
            for node_idx2 in range(self._m2._layer_state_num[snp_idx]):
                for node_idx3 in range(self._m3._layer_state_nums[snp_idx]):
                    for node_idx4 in range(self._m4._layer_state_nums[snp_idx]):
                        if snp_idx > 0:
                            self._backward_scale_factor[snp_idx] += self._backward_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4]
                        else:
                            if self._ibd == 0:
                                self._backward_scale_factor[snp_idx] += self._backward_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4] * \
                                self._emission_prob[0][node_idx1][node_idx2][node_idx3][node_idx4] * \
                                self._m1._pi[node_idx1] * self._m2._pi[node_idx2] * self._m3._pi[node_idx3] * self._m4._pi[node_idx4]
                            else:
                                if self._m1._anc == self._m3._anc:
                                    self._backward_scale_factor[snp_idx] += self._backward_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4] * \
                                    self._emission_prob[0][node_idx1][node_idx2][node_idx3][node_idx4] * \
                                    self._m1._pi[node_idx1] * self._m2._pi[node_idx2] * self._m4._pi[node_idx4]
                                                
                            if self._backward_scale_factor[snp_idx] > 0:                 
                                self._backward_scale_factor[snp_idx] = 1.0 / self._backward_scale_factor[snp_idx]
                            else: 
                                self._backward_scale_factor[snp_idx] = DBL_MAX
                            
        for node_idx1 in range(self._m1._layer_state_nums[snp_idx]):
            for node_idx2 in range(self._m2._layer_state_nums[snp_idx]):
                for node_idx3 in range(self._m3._layer_state_nums[snp_idx]):
                    for node_idx4 in range(self._m4._layer_state_nums[snp_idx]):
                        self._backward_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4] = \
                        self._backward_prob[snp_idx][node_idx1][node_idx2][node_idx3][node_idx4] * self._backward_scale_factor[snp_idx]
    
    cpdef double calc_likelihood(self):
        cdef int snp_idx       
        cdef likelihood = 0
        for snp_idx in range(self._m1._snp_num):
            likelihood = likelihood - log(self._scale_factor[snp_idx])
        return likelihood
            
