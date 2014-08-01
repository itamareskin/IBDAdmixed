#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

from __future__ import division
import random
import os
import math
from libc.stdlib cimport malloc, free
from itertools import islice, combinations_with_replacement
from libc.math cimport exp, log
from libc.float cimport DBL_MIN, DBL_MAX
from libc.limits cimport ULONG_MAX,LONG_MIN
from libcpp cimport bool
from sys import stdout
import numpy as np
cimport numpy as np
 
cdef extern from "string.h":
    char *strncpy(char *dest, char *src, size_t n)
    int strlen(char *s) 

cdef extern from "structs.h":
    int c_isfinite(double x)
    
    cdef struct struct_state:
        double prob_em[2]
        int out_trans_num
        int in_trans_num
        bool likely_allele

    ctypedef struct_state state
    
    state create_state(double allele_0_prob, double allele_1_prob, int out_trans_num, int in_trans_num)
    
    cdef struct struct_transition:
        double prob
        int next_node

    ctypedef struct_transition transition
    
    transition create_transition(double prob, int next_node)
    
    cdef struct struct_gen_map_entry:
        double recomb_rate
        double genetic_dist

    ctypedef struct_gen_map_entry gen_map_entry
    
    gen_map_entry create_gen_map_entry(int position, double recomb_rate, double genetic_dist)
    
    bool get_likely_allele(state s)
    
    double logsumexp(double first, double second)

cdef read_until_blank_line(file_name, first_line = None):
        buffer_size = 500
        blank_line = '\n'
        # read a buffer_size of lines from the file
        if first_line == None:
            lines = list(islice(file_name, buffer_size))
        else:
            lines = list(islice(file_name, first_line, buffer_size))
        # read some more lines until a blank line is reached
        while True:
            line = list(islice(file_name, 1))
            if not line:
                lines += [blank_line]
                break
            lines += line
            if line == [blank_line]:
                break
        return lines

cdef class LDModel(object):
    '''
    Hidden Markov Model for a single ancestral population
    '''
     
    def __cinit__(self, map_file_name=None, beagle_file_name=None, int anc=0, double alpha=1, double t_0_1 = 1e-5, double t_1_0 = 1, int max_snp_num=1000000000, double eps = 1e-4):
        self._anc = anc
        self._alpha = alpha
        self._t_0_1 = t_0_1
        self._t_1_0 = t_1_0
        self.eps = eps
        self._snp_num = 0
        if map_file_name is not None:
            self._gm = GeneticMap(map_file_name, max_snp_num)
            self._snp_num = self._gm._snp_num
        if beagle_file_name is not None:
            self.read_from_bgl_file(beagle_file_name)
        self._is_slice = False
        
    cpdef LDModel get_slice_model(self, int start_snp, int snp_num):
        cdef LDModel other = LDModel()
        other._anc = self._anc
        other.eps = self.eps
        other._gm = self._gm.get_slice(start_snp,snp_num)
        other._snp_num = min(snp_num, self._snp_num - start_snp)
        other._snp_ids = self._snp_ids + start_snp
        other._states = self._states + start_snp
        other._layer_state_nums = self._layer_state_nums + start_snp
        other._pi = <double *> malloc(other._layer_state_nums[0] * sizeof(double))
        if other._layer_state_nums[0] == 0:
            exit(-1)
        for i in range(other._layer_state_nums[0]):
            other._pi[i] = 1
        other._trans = self._trans + start_snp
        other._trans_idx = self._trans_idx + start_snp
        other._back_trans = self._back_trans + start_snp
        other._back_trans_idx = self._back_trans_idx + start_snp
        other._allele_0 = self._allele_0
        other._allele_1 = self._allele_1
        other._is_slice = True
        return other 
    
    def __dealloc__(self):
        cdef int snp_idx
        cdef int node_idx
        free(self._pi)
        if not self._is_slice:
            for snp_idx in range(self._snp_num):
                free(self._states[snp_idx])
                for node_idx in range(self._layer_state_nums[snp_idx]):
                    if snp_idx < self._snp_num - 1:
                        free(self._trans[snp_idx][node_idx])
                        free(self._trans_idx[snp_idx][node_idx])
                    if snp_idx > 0:
                        free(self._back_trans[snp_idx][node_idx])
                        free(self._back_trans_idx[snp_idx][node_idx])
                if snp_idx < self._snp_num - 1:
                    free(self._trans[snp_idx])
                    free(self._trans_idx[snp_idx])
                if snp_idx > 0:
                    free(self._back_trans[snp_idx])
                    free(self._back_trans_idx[snp_idx])
            free(self._states)
            free(self._layer_state_nums)
            free(self._trans)
            free(self._trans_idx)
            free(self._back_trans)
            free(self._back_trans_idx)
    
    def read_from_bgl_file(self, file_name):
        ''' 
        read the model from a beagle model file (.dag file)
        '''
        
        # allocate memory
        self._snp_ids = <char **> malloc(self._snp_num * sizeof(char *))
        self._states = <state **> malloc(self._snp_num * sizeof(state *))
        self._layer_state_nums = <int *> malloc(self._snp_num * sizeof(int))
        self._trans = <double ***> malloc((self._snp_num-1) * sizeof(double **))
        self._trans_idx = <int ***> malloc((self._snp_num-1) * sizeof(int **))
        self._back_trans = <double ***> malloc((self._snp_num) * sizeof(double **))
        self._back_trans_idx = <int ***> malloc((self._snp_num) * sizeof(int **))
        
        if not os.path.exists(file_name):
            print "the file: " + file_name + " does not exist!"  
            return
        
        print "reading from bgl model file: " + file_name
        
        # identify allele coding
        with open(file_name) as model_file:
            model_file.readline()
            model_file.readline()
            line = model_file.readline()
            node = line.split("\t")
            self._allele_0 = int(node[4])
            while True:
                line = model_file.readline()
                if len(line) > 0:
                    node = line.split("\t")
                    if len(node) > 4:
                        curr_allele = node[4]
                        if int(curr_allele) != self._allele_0:
                            if int(curr_allele) > self._allele_0:
                                self._allele_1 = int(curr_allele)
                            else:
                                self._allele_1 = self._allele_0
                                self._allele_0 = int(curr_allele)
                            break;
        
        print "allele 0 symbol: " + str(self._allele_0) + "  allele 1 symbol: " + str(self._allele_1)
                
        with open(file_name) as model_file:
            
            blank_line = '\n'
            layer = 0
            nodes = []
            nodes_prev = []
            nodes_curr = []
            first_read = True
            done = False
            
            while True:
                
                if done:
                    break
                # read next buffer_size lines from the file
                if first_read: 
                    lines = read_until_blank_line(model_file,2)
                    first_read = False
                else:
                    lines = read_until_blank_line(model_file)
                if lines == [blank_line]:   
                    break
                
                for line in lines:
                    if line != blank_line: # split next line
                        node = line.split("\t")
                        nodes.append(node[:6])
                    else: # an entire level has been read
                        nodes_prev = nodes_curr
                        nodes_curr = nodes
                        nodes = []
                        
                        # allocate memory for transition probabilities
                        if layer == 0:
                            self._pi = <double *> malloc(len(nodes_curr) * sizeof(double))
                        if layer > 0:
                            self._trans[layer-1] = <double **> malloc(len(nodes_prev) * sizeof(double*))
                            self._trans_idx[layer-1] = <int **> malloc(len(nodes_prev) * sizeof(int*))
                        self._back_trans[layer] = <double **> malloc(len(nodes_curr) * sizeof(double*))
                        self._back_trans_idx[layer] = <int **> malloc(len(nodes_curr) * sizeof(int*))
                        
                        # set layer state nums and snp ids
                        self._layer_state_nums[layer] = len(nodes_curr)
                        self._snp_ids[layer] = nodes_curr[0][1]
                        
                        # allocate memory for the states
                        self._states[layer] = <state *> malloc(len(nodes_curr) * sizeof(state))
                        
                        for j in range(len(nodes_curr)):
                            if layer == 0:
                                total_parent_count = sum([int(x[5]) for x in nodes_curr])
                                self._pi[j] = int(nodes_curr[j][5]) / total_parent_count
                            
                        for j in range(len(nodes_curr)):
                            # find how many edges go into this node from the previous layer
                            edges_num = 0
                            for k in range(len(nodes_prev)):
                                if nodes_prev[k][3] == nodes_curr[j][2]:
                                    edges_num += 1
                            # create states and set emission probability
                            error_eps = self.eps*1e-1  
                            if int(nodes_curr[j][4]) == self._allele_0:
                                self._states[layer][j] = create_state(1 - error_eps, error_eps, 0, edges_num)
                            else: 
                                self._states[layer][j] = create_state(error_eps, 1 - error_eps, 0, edges_num)
                            
                            if layer > 0:
                                # allocate memory for the edges
                                self._back_trans[layer][j] = <double *> malloc(edges_num * sizeof(double))
                                self._back_trans_idx[layer][j] = <int *> malloc(edges_num * sizeof(int))
                                # set the edges transition probabilities
                                edge_idx = 0
                                sum_back_trans = 0
                                for k in range(len(nodes_prev)):
                                    if nodes_prev[k][3] == nodes_curr[j][2]:
                                        total_parent_count = sum([int(x[5]) for x in nodes_curr if x[2] == nodes_curr[j][2]])
                                        self._back_trans[layer][j][edge_idx] = int(nodes_curr[j][5]) / total_parent_count
                                        self._back_trans_idx[layer][j][edge_idx] = k
                                        sum_back_trans += self._back_trans[layer][j][edge_idx]
                                        edge_idx += 1
                        
                        if layer > 0:
                            for j in range(len(nodes_prev)):
                                # find how many edges go into this node from the previous layer
                                edges_num = 0
                                for k in range(len(nodes_curr)):
                                    if nodes_prev[j][3] == nodes_curr[k][2]:
                                        edges_num += 1
                                        
                                # set states
                                self._states[layer-1][j].out_trans_num = edges_num
                                # allocate memory for the edges
                                self._trans[layer-1][j] = <double *> malloc(edges_num * sizeof(double))
                                self._trans_idx[layer-1][j] = <int *> malloc(edges_num * sizeof(int))
                                # set the edges transition probabilities
                                edge_idx = 0
                                for k in range(len(nodes_curr)):
                                    if nodes_prev[j][3] == nodes_curr[k][2]:
                                        total_parent_count = sum([int(x[5]) for x in nodes_curr if x[2] == nodes_curr[k][2]])
                                        self._trans[layer-1][j][edge_idx] = int(nodes_curr[k][5]) / total_parent_count
                                        self._trans_idx[layer-1][j][edge_idx] = k
                                        edge_idx += 1
                        
                        layer += 1
                        if layer >= self._snp_num:
                            done = True;
                            break
                    
        print "Finished reading beagle file."
    
    def get_layer_node_nums(self, anc):
        node_nums = []
        for layer in range(self._snp_num):
            node_nums.append(self._layer_state_nums[layer])
        return node_nums 
    
    def get_node_edges(self, anc):
        edges = []
        edge_weights = []
        for layer in range(self._snp_num-1):
            layer_edges = []
            layer_edge_weights = []
            for node in range(self._layer_state_nums[layer]):
                layer_edges.append([self._trans_idx[layer][node][i] for i in xrange(self._states[layer][node].out_trans_num)])
                layer_edge_weights.append([self._trans[layer][node][i] for i in xrange(self._states[layer][node].out_trans_num)])
            edges.append(layer_edges)
            edge_weights.append(layer_edge_weights)
        return (edges, edge_weights)
    
    def get_node_edges_back(self):
        edges = []
        edge_weights = []
        for layer in range(1,self._snp_num):
            layer_edges = []
            layer_edge_weights = []
            for node in range(self._layer_state_nums[layer]):
                layer_edges.append([self._back_trans_idx[layer][node][i] for i in xrange(self._states[layer][node].in_trans_num)])
                layer_edge_weights.append([self._back_trans[layer][node][i] for i in xrange(self._states[layer][node].in_trans_num)])
            edges.append(layer_edges)
            edge_weights.append(layer_edge_weights)
        return (edges, edge_weights)
    
    def get_node_ems_probs(self):
        ems_probs = []
        for layer in range(self._snp_num):
            layer_ems_probs = []
            for node in range(self._layer_state_nums[layer]):
                layer_ems_probs.append((self._states[layer][node].prob_em[0],self._states[layer][node].prob_em[1]))
            ems_probs.append(layer_ems_probs)
        return ems_probs
    
    def get_trans_probs(self):
        probs = []
        for layer in range(self._snp_num-1):
            layer_probs = []
            for node in range(self._layer_state_nums[layer]):
                layer_probs.append([self._trans[layer][node][i] for i in xrange(self._states[layer][node].out_trans_num)])
            probs.append(layer_probs)
        return probs
    
    def print_emissions(self):
        self._ems_file.write("snp node ems_prob0 ems_prob1\n")
        for snp_idx in range(self._snp_num):
            for node_idx in range(self._layer_state_nums[snp_idx]):
                self._ems_file.write(str(snp_idx) + " " + 
                              str(node_idx) + " " + 
                              str(self._states[snp_idx][node_idx].prob_em[0]) + " " + 
                              str(self._states[snp_idx][node_idx].prob_em[1]) + "\n") 
    
    def print_transitions(self):
        self._trans_file.write("snp node next_node transition\n")
        for snp_idx in range(self._snp_num):
            for node_idx in range(self._layer_state_nums[snp_idx]):
                print str(snp_idx) + " " + str(node_idx) + " " + str(self._states[snp_idx][node_idx].out_trans_num)
                for nxt_node_idx in range(self._states[snp_idx][node_idx].out_trans_num):
                    self._trans_file.write(str(snp_idx) + " " + 
                              str(node_idx) + " " +
                              str(nxt_node_idx) + " " + 
                              str(self._trans[snp_idx][node_idx][nxt_node_idx]) + "\n") 
        self._trans_file.flush()

    cpdef int start_position(self):
        return self._gm._position[0] 
    
    cpdef int end_position(self):
        return self._gm._position[self._snp_num-1]
    
    cpdef int get_position(self, int snp_num):
        return self._gm._position[snp_num]
    
    cdef bool* generate_random_hap(self, int max_snp_num):
        cdef int snp_num = self._snp_num if self._snp_num <= max_snp_num else max_snp_num
        cdef bool* hap = <bool *> malloc(self._snp_num * sizeof(bool))
        p = [self._pi[i] for i in xrange(self._layer_state_nums[0])]
        cdef int next_node = np.array(p).cumsum().searchsorted(np.random.sample(1))[0]
        for layer in range(self._snp_num):
            if self._states[layer][next_node].prob_em[0] > 0.5:
                hap[layer] = False
            else:
                hap[layer] = True
            if layer < self._snp_num - 1:
                t = [self._trans[layer][next_node][i] for i in xrange(self._states[layer][next_node].out_trans_num)]
                next_node = self._trans_idx[layer][next_node][np.array(t).cumsum().searchsorted(np.random.sample(1))]
        return hap
        
