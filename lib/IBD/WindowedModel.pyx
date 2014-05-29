# cython: profile=True
# cython: boundscheck=False
# cython: cdivision=True
# cython.wraparound=False
# cython.nonecheck=False

from __future__ import division
import os
from libc.stdlib cimport malloc, free
from InnerModel import InnerModel
from InnerModel cimport InnerModel
from TestSet import TestSet
from TestSet cimport TestSet
from libcpp cimport bool
from libc.float cimport DBL_MIN, DBL_MAX
from libc.math cimport exp, log
import logging
posterior_probs_logger = logging.getLogger('posteriorprobs')

cdef extern from "structs.h":
    double logsumexp(double first, double second)

cdef class WindowedModel(object):
    
    def __cinit__(self, models, snp_num, win_size):
        self._model_num = len(models)
        self._models = models 
        for model_idx in range(self._model_num):
            if not isinstance(models[model_idx], InnerModel):
                raise TypeError
        self._snp_num = snp_num
        self._win_size = win_size
    
    def __dealloc__(self):
        cdef int win_idx                 
        free(self._scale_factor)
        free(self._posterior_prob)
        for win_idx in range(self.get_num_windows()):
            free(self._ems_prob[win_idx])
            free(self._forward_prob[win_idx])
            free(self._backward_prob[win_idx])
        free(self._ems_prob)
        free(self._forward_prob)
        free(self._backward_prob)
        
    cpdef alloc_mem(self):
        cdef int win_idx
        cdef int model_idx 
         
        self._ems_prob = < double **> malloc(self.get_num_windows() * sizeof(double *))
        self._forward_prob = < double **> malloc(self.get_num_windows() * sizeof(double *))
        self._backward_prob = < double **> malloc(self.get_num_windows() * sizeof(double *))
        self._scale_factor = < double *> malloc(self.get_num_windows() * sizeof(double))
        self._posterior_prob = < double *> malloc(self.get_num_windows() * sizeof(double))
        for win_idx in range(self.get_num_windows()):
#             print "win_idx " + str(win_idx) + " out of " + str(self.get_num_windows())
            self._scale_factor[win_idx] = 0
            self._posterior_prob[win_idx] = 0
            self._ems_prob[win_idx] = < double *> malloc(self._model_num * sizeof(double))
            self._forward_prob[win_idx] = < double *> malloc(self._model_num * sizeof(double))
            self._backward_prob[win_idx] = < double *> malloc(self._model_num * sizeof(double))
            for model_idx in range(self._model_num):
#                 print "model_idx " + str(model_idx) + " out of " + str(self._model_num)
                self._ems_prob[win_idx][model_idx] = 0
                self._forward_prob[win_idx][model_idx] = -DBL_MAX
                self._backward_prob[win_idx][model_idx] = -DBL_MAX
    
    cpdef free_mem(self):
        cdef int win_idx
        free(self._scale_factor)
        free(self._posterior_prob)
        for win_idx in range(self.get_num_windows()):
            free(self._ems_prob[win_idx])
            free(self._forward_prob[win_idx])
            free(self._backward_prob[win_idx])
        free(self._ems_prob)
        free(self._forward_prob)
        free(self._backward_prob)
        
    cpdef alloc_mem_viterbi(self):
        cdef int win_idx
        cdef int model_idx 
         
        self._ems_prob = < double **> malloc(self.get_num_windows() * sizeof(double *))
        self._viterbi_prob = < double **> malloc(self.get_num_windows() * sizeof(double *))
        self._viterbi_back_track = < int **> malloc(self.get_num_windows() * sizeof(int *))
        self._viterbi_path = < int *> malloc(self.get_num_windows() * sizeof(int))
        for win_idx in range(self.get_num_windows()):
#             print "win_idx " + str(win_idx) + " out of " + str(self.get_num_windows())
            self._viterbi_path[win_idx] = 0
            self._ems_prob[win_idx] = < double *> malloc(self._model_num * sizeof(double))
            self._viterbi_prob[win_idx] = < double *> malloc(self._model_num * sizeof(double))
            self._viterbi_back_track[win_idx] = < int *> malloc(self._model_num * sizeof(int))
            for model_idx in range(self._model_num):
#                 print "model_idx " + str(model_idx) + " out of " + str(self._model_num)
                self._ems_prob[win_idx][model_idx] = 0
                self._viterbi_prob[win_idx][model_idx] = -DBL_MAX
                self._viterbi_back_track[win_idx][model_idx] = -1
    
    cpdef free_mem_viterbi(self):
        cdef int win_idx
        free(self._viterbi_path)
        for win_idx in range(self.get_num_windows()):
            free(self._ems_prob[win_idx])
            free(self._viterbi_prob[win_idx])
            free(self._viterbi_back_track[win_idx])
        free(self._ems_prob)
        free(self._viterbi_prob)
        free(self._viterbi_back_track)
    
    cpdef calc_ems_probs(self, TestSet obs_data):
        cdef int win_idx
        cdef int model_idx
        cdef double tmp
        cdef InnerModel inner_model
         
        for win_idx in range(self.get_num_windows()):
            for model_idx in range(self._model_num):
                inner_model = self._models[model_idx].slice_from_model(self.start_snp(win_idx), self._win_size)
                print "calculating likelihood in window " + str(win_idx)
                print "snp num in slice: " + str(obs_data.get_slice(self.start_snp(win_idx), self._win_size)._snp_num) 
                print "haplos num in slice: " + str(obs_data.get_slice(self.start_snp(win_idx), self._win_size)._nr_haplos)
                self._ems_prob[win_idx][model_idx] = inner_model.calc_likelihood(obs_data.get_slice(self.start_snp(win_idx), self._win_size))
                print "finished calculating likelihood in window " + str(win_idx)
    
    cpdef calc_forward_probs(self):
        cdef int win_idx 
        cdef int model_idx
        cdef int prev_model_idx
        cdef InnerModel inner_model
        cdef InnerModel prev_inner_model
       
        for model_idx in range(self._model_num):
            self._forward_prob[0][model_idx] = \
            logsumexp(self._forward_prob[0][model_idx], self._ems_prob[0][model_idx]) # + self._models[model_idx]._prior)
                                  
        for win_idx in range(self.get_num_windows() - 1):
            for model_idx in range(self._model_num):
                for prev_model_idx in range(self._model_num):
                    inner_model = self._models[model_idx].slice_from_model(self.start_snp(win_idx + 1), self._win_size)
                    prev_inner_model = self._models[prev_model_idx].slice_from_model(self.start_snp(win_idx), self._win_size)
                    trans = prev_inner_model.trans_prob(inner_model)
                    self._forward_prob[win_idx + 1][model_idx] = \
                    logsumexp(self._forward_prob[win_idx + 1][model_idx], \
                              self._forward_prob[win_idx][prev_model_idx] + \
                              self._ems_prob[win_idx + 1][model_idx] + \
                              log(trans))
                    
    cpdef calc_backward_probs(self):
        cdef int win_idx
        cdef int model_idx
        cdef int nxt_model_idx
        cdef InnerModel inner_model
        cdef InnerModel nxt_inner_model
         
        for model_idx in range(self._model_num):
            self._backward_prob[self.get_num_windows() - 1][model_idx] = 0

        for win_idx in reversed(range(self.get_num_windows() - 1)):                
            for model_idx in range(self._model_num):
                for nxt_model_idx in range(self._model_num):
                    inner_model = self._models[model_idx].slice_from_model(self.start_snp(win_idx), self._win_size)
                    nxt_inner_model = self._models[nxt_model_idx].slice_from_model(self.start_snp(win_idx + 1), self._win_size)
                    trans = inner_model.trans_prob(nxt_inner_model)
                    self._backward_prob[win_idx][model_idx] = \
                    logsumexp(self._backward_prob[win_idx][model_idx], \
                    self._backward_prob[win_idx + 1][nxt_model_idx] + \
                    self._ems_prob[win_idx + 1][nxt_model_idx] + \
                    log(trans))
                    
    cpdef posterior_decoding(self):
        cdef int win_idx
        cdef int model_idx
        cdef double curr_gamma
        
        for win_idx in range(self.get_num_windows()):
            self._posterior_prob[win_idx] = -DBL_MAX
            for model_idx in range(self._model_num):
                curr_gamma = self._forward_prob[win_idx][model_idx] + \
                self._backward_prob[win_idx][model_idx]
                self._posterior_prob[win_idx] = logsumexp(self._posterior_prob[win_idx], curr_gamma)
            posterior_probs_logger.debug("win: " + str(win_idx) + " ibd: " + str(self._models[0]._ibd) + " posterior prob: " + str(self._posterior_prob[win_idx]))
            
    cpdef viterbi_decoding(self):
        cdef int win_idx 
        cdef int model_idx
        cdef int prev_model_idx
        cdef InnerModel inner_model
        cdef InnerModel prev_inner_model
        cdef double curr_prob
        cdef double max_prob
        cdef int max_prev
        cdef int max_idx
        cdef int curr_model_
       
        for model_idx in range(self._model_num):
            self._viterbi_prob[0][model_idx] = self._ems_prob[0][model_idx] + self._models[model_idx]._prior
                                  
        for win_idx in range(self.get_num_windows() - 1):
            print "viterbi probs, win: " + str(win_idx)
            for model_idx in range(self._model_num):
                max_prob = -DBL_MAX
                max_prev = -1
                for prev_model_idx in range(self._model_num):
                    inner_model = self._models[model_idx].slice_from_model(self.start_snp(win_idx + 1), self._win_size)
                    prev_inner_model = self._models[prev_model_idx].slice_from_model(self.start_snp(win_idx), self._win_size)
                    trans = prev_inner_model.trans_prob(inner_model)
                    curr_prob = self._forward_prob[win_idx][prev_model_idx] + \
                                self._ems_prob[win_idx + 1][model_idx] + \
                                log(trans)
                    if curr_prob > max_prob:
                        max_prob = curr_prob
                        max_prev = prev_model_idx
                self._viterbi_prob[win_idx + 1][model_idx] = max_prob
                self._viterbi_back_track[win_idx + 1][model_idx] = max_prev
        
        max_prob = -DBL_MAX
        max_prev = -1
        for model_idx in range(self._model_num):
            curr_prob = self._viterbi_prob[self.get_num_windows() - 1][model_idx]
            if curr_prob > max_prob:
                max_prob = curr_prob
                max_idx = model_idx
                
        self._viterbi_path[self.get_num_windows() - 1] = max_idx
        for win_idx in reversed(range(self.get_num_windows() - 1)):
            print "viterbi path, win: " + str(win_idx)
            self._viterbi_path[win_idx] = self._viterbi_back_track[win_idx+1][self._viterbi_path[win_idx+1]]
                
    def get_posterior_probs(self):
        cdef int win_idx
        cdef list res = [] 
        for win_idx in range(self.get_num_windows()):
            res.append(self._posterior_prob[win_idx])
        return res 
    
    def get_viterbi_path_models(self):
        cdef int win_idx
        cdef list path = []
        for win_idx in range(self.get_num_windows()):
            path.append(self._models[self._viterbi_path[win_idx]])
        return path
    
    cdef int get_num_windows(self):
        cdef int num_win
        num_win = int(self._snp_num / self._win_size)
        if self._snp_num % self._win_size > 0:
            num_win += 1
        return num_win
    
    cdef int start_snp(self, int win_idx):
        return max(0, win_idx * self._win_size)
    
    cdef int end_snp(self, int win_idx):
        return min((win_idx + 1) * self._win_size, self._snp_num)
    
    cpdef compare(self, WindowedModel other):
        lod_scores = []
        if self.get_num_windows() != other.get_num_windows():
            raise ValueError 
        for win_idx in range(self.get_num_windows()):
            lod_scores.append(2*(self._posterior_prob[win_idx] - other._posterior_prob[win_idx]))
        return lod_scores
 
