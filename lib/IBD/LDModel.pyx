'''
Created on May 15, 2012
 
@author: eskin3
'''

from __future__ import division
import random
import numpy as np
cimport numpy as np
import os
import math
from libc.stdlib cimport malloc, free
from itertools import islice
from libc.math cimport exp, log
from libc.float cimport DBL_MIN, DBL_MAX
from IBD.cIBD import cPairIBD

cdef double eps = 1e-4

cdef extern from "string.h":
    char *strncpy(char *dest, char *src, size_t n)

cdef extern from "structs.h":
    int c_isfinite(double x)
    
    cdef struct struct_state:
        double prob_em[2]
        int out_trans_num
        int in_trans_num

    ctypedef struct_state state
    
    state create_state(double allele_0_prob, double allele_1_prob, int out_trans_num, int in_trans_num)
    
    cdef struct struct_transition:
        double prob
        int next_node

    ctypedef struct_transition transition
    
    transition create_transition(double prob, int next_node)
    
    cdef struct struct_gen_map_entry:
        int position
        double recomb_rate
        double genetic_dist

    ctypedef struct_gen_map_entry gen_map_entry
    
    gen_map_entry create_gen_map_entry(int position, double recomb_rate, double genetic_dist)
    
    char get_likely_allele(state s)

cdef read_until_blank_line(file_name, first_line = None):
        buffer_size = 100
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
     
    def __cinit__(self, int snp_num, int k, int g, int win_size):
        
        # total number of SNPs to be analyzed
        self._snp_num = snp_num
        
        # number of haplotypes to analyze
        self._nr_haplos = 0
        
        # number of ancestral populations
        self.K = k 
        
        # number of generations since admixture began
        self.g = g 
        
        # window size (in number of consecutive SNPs). 
        # Inside each window, the ancestral origin of each haplotype is assumed to be unique, 
        # that is - transitions from one ancestry to the other are only allowd between consecutive windows
        self._win_size = win_size 
        
        # probabilities of transition from IBD to No-IBD and vice-versa (notation from the Browning paper)
        self._t_0_1 = 1e-4
        self._t_1_0 = 1
        
        # prior probability of IBD
        self._ibd_prior = <double *> malloc(2 * sizeof(double))
        self._ibd_prior[1] = self._t_0_1 / (self._t_0_1 + self._t_1_0)
        self._ibd_prior[0] = 1-self._ibd_prior[1]
        
        # proportions of ancestry from each ancestral populations
        self._alphas = <double *> malloc(self.K * sizeof(double))
        for i in range(self.K):
            self._alphas[i] = (1.0/self.K)
            
        # allocate memory
        self._snp_ids = <char **> malloc(self._snp_num * sizeof(char *))
        
        self._states = <state ***> malloc(self.K * sizeof(state **))
        self._pi = <double ***> malloc(self.K * sizeof(double **))
        self._layer_state_nums = <int **> malloc(self.K * sizeof(int*))
        self._trans = <double ****> malloc(self.K * sizeof(double ***))
        self._trans_idx = <int ****> malloc(self.K * sizeof(int ***))
        self._back_trans = <double ****> malloc(self.K * sizeof(double ***))
        self._back_trans_idx = <int ****> malloc(self.K * sizeof(int ***))
        for i in range(self.K):
            self._states[i] = <state **> malloc(self._snp_num * sizeof(state **))
            self._layer_state_nums[i] = <int *> malloc(self._snp_num * sizeof(int))
            self._trans[i] = <double ***> malloc((self._snp_num-1) * sizeof(double **))
            self._trans_idx[i] = <int ***> malloc((self._snp_num-1) * sizeof(int **))
            self._back_trans[i] = <double ***> malloc((self._snp_num) * sizeof(double **))
            self._back_trans_idx[i] = <int ***> malloc((self._snp_num) * sizeof(int **))
            self._pi[i] = <double **> malloc(self.get_num_windows() * sizeof(double *))
    
    
    cpdef set_alphas(self, alphas):
        if len(alphas) != self.K:
            print "wrong number of alphas"
            return 
        for i in range(self.K):
            self._alphas[i] = alphas[i]
    
    
    def read_genetic_map(self, genetic_map_file_name):
        '''
        read the genetic map from hapmap format file
        '''
        if not os.path.exists(genetic_map_file_name):
            print "the file: " + genetic_map_file_name + " does not exist!"
            return
        
        print "reading from genetic map file: " + genetic_map_file_name
        
        # allocate memory
        self._genetic_map = <gen_map_entry *> malloc(self._snp_num * sizeof(gen_map_entry))
        print "allocated mem"
        
        with open(genetic_map_file_name) as genetic_map_file:
            
            first_read = True
            done = False
            buffer_size = 100
            snp_idx = 0
            
            while True:
                
                if done:
                    break
                # read next buffer_size lines from the file
                if not first_read:
                    lines = list(islice(genetic_map_file, buffer_size))
                else:
                    lines = list(islice(genetic_map_file, 1, buffer_size))
                    first_read = False
                
                if len(lines) == 0:
                    done = True
                
                for line in lines:
                    #print "line: " + line
                    line = line.split(" ")
                    #print "pos: " + line[0]
                    self._genetic_map[snp_idx] = create_gen_map_entry(int(line[0]),float(line[1]),float(line[2]))
                    snp_idx += 1
                    if snp_idx >= self._snp_num:
                        done = True;
                        break
        
    
    def read_haplos(self, file_name, nr_haplos):
        if not os.path.exists(file_name):
            print "the file: " + file_name + " does not exist!"
            return
        
        self._nr_haplos = nr_haplos
        
        print "reading from haplos file: " + file_name
        self._haplos = <char **> malloc(nr_haplos * sizeof(char *)) 
        
        with open(file_name) as bgl_file:
            
            first_read = True
            done = False
            buffer_size = 100
            hap_idx = 0
            
            while True:
                
                if done:
                    break
                # read next buffer_size lines from the file
                lines = list(islice(bgl_file, buffer_size))
                
                if len(lines) == 0:
                    done = True
                
                for line in lines:
                    #print "line: " + line
                    #print "pos: " + line[0]
                    line_trunc = line[:self._snp_num] 
                    self._haplos[hap_idx] = <char *> malloc((self._snp_num) * sizeof(char))
                    strncpy(self._haplos[hap_idx], line_trunc,self._snp_num) 
                    #py_bytes = line_trunc.encode('UTF-8')
                    #self._haplos[hap_idx] = py_bytes
                    #print self._haplos[hap_idx]
                    hap_idx += 1
                    if hap_idx >= nr_haplos:
                        done = True;
                        break
    
    def read_from_bgl_file(self, file_name, anc):
        '''
        read the model from a beagle model file (.dag file)
        '''
        
        if not os.path.exists(file_name):
            print "the file: " + file_name + " does not exist!"  
            return
        
        print "reading from bgl model file: " + file_name
        
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
                        print "reading layer: " + str(layer) + ", snp id: " + nodes_curr[0][1]
                        
                        # allocate memory for transition probabilities
                        if layer % self._win_size == 0:
                            self._pi[anc][int(layer / self._win_size)] = <double *> malloc(len(nodes_curr) * sizeof(double))
                        if layer > 0:
                            self._trans[anc][layer-1] = <double **> malloc(len(nodes_prev) * sizeof(double*))
                            self._trans_idx[anc][layer-1] = <int **> malloc(len(nodes_prev) * sizeof(int*))
                        self._back_trans[anc][layer] = <double **> malloc(len(nodes_curr) * sizeof(double*))
                        self._back_trans_idx[anc][layer] = <int **> malloc(len(nodes_curr) * sizeof(int*))
                        
                        # set layer state nums and snp ids
                        self._layer_state_nums[anc][layer] = len(nodes_curr)
                        self._snp_ids[layer] = nodes_curr[0][1]
                        
                        # allocate memory for the states
                        self._states[anc][layer] = <state *> malloc(len(nodes_curr) * sizeof(state))
                        
                        for j in range(len(nodes_curr)):
                            if layer == 0:
                                total_parent_count = sum([int(x[5]) for x in nodes_curr])
                                self._pi[anc][0][j] = int(nodes_curr[j][5]) / total_parent_count
                                #print "layer: " + str(layer) + " pi in node " + str(j) + " is: " + str(self._pi[anc][0][j])
                                #if j > 0:
                                #    print "  pi in previous node " + str(j-1) + " is: " + str(self._pi[anc][0][j-1])
                            
                        for j in range(len(nodes_curr)):
                            # find how many edges go into this node from the previous layer
                            edges_num = 0
                            for k in range(len(nodes_prev)):
                                if nodes_prev[k][3] == nodes_curr[j][2]:
                                    edges_num += 1
                            # create states and set emission probability
                            error_eps = 0
                            if int(nodes_curr[j][4]) == 0:
                                self._states[anc][layer][j] = create_state(1 - error_eps, error_eps, 0, edges_num)
                                #self._states[anc][layer][j] = create_state(1, 0, 0, edges_num)
                            else: 
                                self._states[anc][layer][j] = create_state(error_eps, 1 - error_eps, 0, edges_num)
                                #self._states[anc][layer][j] = create_state(0, 1, 0, edges_num)
                            
                            if layer > 0:
                                # allocate memory for the edges
                                self._back_trans[anc][layer][j] = <double *> malloc(edges_num * sizeof(double))
                                self._back_trans_idx[anc][layer][j] = <int *> malloc(edges_num * sizeof(int))
                                # set the edges transition probabilities
                                edge_idx = 0
                                sum_back_trans = 0
                                for k in range(len(nodes_prev)):
                                    if nodes_prev[k][3] == nodes_curr[j][2]:
                                        total_parent_count = sum([int(x[5]) for x in nodes_curr if x[2] == nodes_curr[j][2]])
                                        self._back_trans[anc][layer][j][edge_idx] = int(nodes_curr[j][5]) / total_parent_count
                                        self._back_trans_idx[anc][layer][j][edge_idx] = k
                                        sum_back_trans += self._back_trans[anc][layer][j][edge_idx]
                                        edge_idx += 1
                                
                                # set initial probabilities for every window
                                if layer % self._win_size == 0:
                                    self._pi[anc][int(layer / self._win_size)][j] = sum_back_trans
                        
                        if layer > 0:
                            for j in range(len(nodes_prev)):
                                # find how many edges go into this node from the previous layer
                                edges_num = 0
                                for k in range(len(nodes_curr)):
                                    if nodes_prev[j][3] == nodes_curr[k][2]:
                                        edges_num += 1
                                        
                                # set states
                                self._states[anc][layer-1][j].out_trans_num = edges_num
                                # allocate memory for the edges
                                self._trans[anc][layer-1][j] = <double *> malloc(edges_num * sizeof(double))
                                self._trans_idx[anc][layer-1][j] = <int *> malloc(edges_num * sizeof(int))
                                # set the edges transition probabilities
                                edge_idx = 0
                                for k in range(len(nodes_curr)):
                                    if nodes_prev[j][3] == nodes_curr[k][2]:
                                        total_parent_count = sum([int(x[5]) for x in nodes_curr if x[2] == nodes_curr[k][2]])
                                        self._trans[anc][layer-1][j][edge_idx] = int(nodes_curr[k][5]) / total_parent_count
                                        self._trans_idx[anc][layer-1][j][edge_idx] = k
                                        edge_idx += 1
                        
                        layer += 1
                        if layer >= self._snp_num:
                            done = True;
                            break
                #p = [self._pi[anc][i] for i in xrange(self._layer_state_nums[anc][0])]
                #print "pi: " + str(p)  
    
    def get_haplo_num(self):
        return self._nr_haplos
    
    def get_snp_num(self):
        return self._snp_num
    
    def get_haplo(self,hap_idx):
        #pystring = self._haplos[hap_idx][:self._snp_num].decode('UTF-8')
        pystring = self._haplos[hap_idx][:self._snp_num]
        return pystring
    
    def get_layer_node_nums(self, anc):
        node_nums = []
        for layer in range(self._snp_num):
            node_nums.append(self._layer_state_nums[anc][layer])
        return node_nums 
    
    def get_node_edges(self, anc):
        edges = []
        edge_weights = []
        for layer in range(self._snp_num-1):
            print "layer: " + str(layer)
            layer_edges = []
            layer_edge_weights = []
            for node in range(self._layer_state_nums[anc][layer]):
                print "node: " + str(node)
                layer_edges.append([self._trans_idx[anc][layer][node][i] for i in xrange(self._states[anc][layer][node].out_trans_num)])
                layer_edge_weights.append([self._trans[anc][layer][node][i] for i in xrange(self._states[anc][layer][node].out_trans_num)])
            edges.append(layer_edges)
            edge_weights.append(layer_edge_weights)
        return (edges, edge_weights)
    
    def get_node_edges_back(self, anc):
        edges = []
        edge_weights = []
        for layer in range(1,self._snp_num):
            print "layer: " + str(layer)
            layer_edges = []
            layer_edge_weights = []
            for node in range(self._layer_state_nums[anc][layer]):
                print "node: " + str(node)
                layer_edges.append([self._back_trans_idx[anc][layer][node][i] for i in xrange(self._states[anc][layer][node].in_trans_num)])
                layer_edge_weights.append([self._back_trans[anc][layer][node][i] for i in xrange(self._states[anc][layer][node].in_trans_num)])
            edges.append(layer_edges)
            edge_weights.append(layer_edge_weights)
        return (edges, edge_weights)
    
    def get_node_ems_probs(self, anc):
        ems_probs = []
        for layer in range(self._snp_num):
            print "layer: " + str(layer)
            layer_ems_probs = []
            for node in range(self._layer_state_nums[anc][layer]):
                print "node: " + str(node)
                layer_ems_probs.append((self._states[anc][layer][node].prob_em[0],self._states[anc][layer][node].prob_em[1]))
            ems_probs.append(layer_ems_probs)
        return ems_probs
    
    def get_trans_probs(self, anc):
        probs = []
        for layer in range(self._snp_num-1):
            layer_probs = []
            for node in range(self._layer_state_nums[anc][layer]):
                layer_probs.append([self._trans[anc][layer][node][i] for i in xrange(self._states[anc][layer][node].out_trans_num)])
            probs.append(layer_probs)
        return probs
    
    def get_forward_probs(self, anc):
        probs = []
        for layer in range(self._snp_num):
            layer_probs = []
            for node in range(self._layer_state_nums[anc][layer]):
                layer_probs.append(self._forward_probs[layer][node])
            probs.append(layer_probs)
        return probs
    
    cpdef generate_random_hap(self, int anc):
        hap = ""
        p = [self._pi[anc][0][i] for i in xrange(self._layer_state_nums[anc][0])]
        #print "layer state nums: " + str(self._layer_state_nums[anc][0])
        #print "pi: " + str(p)  
        cdef int next_node = np.array(p).cumsum().searchsorted(np.random.sample(1))[0]  
        #print "next node: " + str(next_node)
        for layer in range(self._snp_num):
            print "generating haplotype for layer: " + str(layer)
            if self._states[anc][layer][next_node].prob_em[0] > 0.5:
                hap += '0'
            else:
                hap += '1'
            if layer < self._snp_num - 1:
                #print "states in layer: " + str(layer)
                t = [self._trans[anc][layer][next_node][i] for i in xrange(self._states[anc][layer][next_node].out_trans_num)]
                #print "transition probs: " + str(t)
                next_node = self._trans_idx[anc][layer][next_node][np.array(t).cumsum().searchsorted(np.random.sample(1))]
                #print "next node: " + str(next_node)
        return hap
        
    cpdef calc_ibd_prior(self):
        cdef int win_idx
        cdef double d
        self._s = <double ***> malloc((self.get_num_windows()) * sizeof(double **))
        for win_idx in range(self.get_num_windows()):
            self._s[win_idx] = <double **> malloc(2 * sizeof(double *))
            for j in range(2):
                self._s[win_idx][j] = <double *> malloc(2 * sizeof(double))
        
        for win_idx in range(self.get_num_windows()):
            d = self._genetic_map[self.end_snp(win_idx)].genetic_dist - self._genetic_map[self.start_snp(win_idx)].genetic_dist
            self._s[win_idx][1][1] = exp(-self._t_1_0 * d)  
            self._s[win_idx][0][0] = exp(-self._t_0_1 * d)
            self._s[win_idx][1][0] = 1 - self._s[win_idx][1][1]
            self._s[win_idx][0][1] = 1 - self._s[win_idx][0][0]
            #for i in range(2):
            #    for j in range(2):
            #        print "ibd trans: " + str(win_idx) + " " + str(i) + " " + str(j) + " : " + str(self._s[win_idx][i][j])
            
        #self._s[self.get_num_windows()-1][0][0] = 1
        #self._s[self.get_num_windows()-1][0][1] = 1
        #self._s[self.get_num_windows()-1][1][0] = 1
        #self._s[self.get_num_windows()-1][1][1] = 1
    
    cpdef calc_anc_trans(self):
        cdef int win_idx
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int nxt_admx_idx1
        cdef int nxt_admx_idx2
        cdef int nxt_admx_idx3
        cdef int nxt_admx_idx4
        
        cdef double *sum_recomb
        sum_recomb = <double *> malloc((self.get_num_windows()) * sizeof(double))
        for win_idx in range(self.get_num_windows()):
            sum_recomb[win_idx] = 0
        for snp_idx in range(self._snp_num):
            sum_recomb[int(snp_idx/self._win_size)] += self._genetic_map[snp_idx].recomb_rate #* 1e-6
        
        self._anc_trans = <double *********> malloc((self.get_num_windows()) * sizeof(double ********))
        for win_idx in range(self.get_num_windows()):
            sum = 0
            #print "calculating ancestry transitions for window: " + str(win_idx)
            self._anc_trans[win_idx] = <double ********> malloc(self.K * sizeof(double *******))
            for admx_idx1 in range(self.K):
                self._anc_trans[win_idx][admx_idx1] = <double *******> malloc(self.K * sizeof(double ******))
                for admx_idx2 in range(admx_idx1+1):
                    self._anc_trans[win_idx][admx_idx1][admx_idx2] = <double ******> malloc(self.K * sizeof(double *****))
                    for admx_idx3 in range(self.K):
                        self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3] = <double *****> malloc(self.K * sizeof(double ****))
                        for admx_idx4 in range(admx_idx3+1):
                            self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4] = <double ****> malloc(self.K * sizeof(double ***))
                            for nxt_admx_idx1 in range(self.K):
                                self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1] = <double ***> malloc(self.K * sizeof(double **))
                                for nxt_admx_idx2 in range(nxt_admx_idx1+1):
                                    self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2] = <double **> malloc(self.K * sizeof(double *))
                                    for nxt_admx_idx3 in range(self.K):
                                        self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3] = <double *> malloc(self.K * sizeof(double))
                                        for nxt_admx_idx4 in range(nxt_admx_idx3+1):
                                            if sum_recomb[win_idx] > 0:
                                                self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4] = (self.g - 1) * sum_recomb[win_idx]
                                            else:
                                                self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4] = 1
                                            sum += self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4]
            scale = 1.0                                        
            if sum == 0:
                scale = 1.0
            else:
                scale = 1.0/sum
                if not c_isfinite(scale):
                    scale = 1.0                     
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for nxt_admx_idx1 in range(self.K):
                                for nxt_admx_idx2 in range(nxt_admx_idx1+1):
                                    for nxt_admx_idx3 in range(self.K):
                                        for nxt_admx_idx4 in range(nxt_admx_idx3+1):
                                            self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4] = \
                                            self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4] * scale
                                            if not c_isfinite(self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4]):
                                                self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4] = 0.0                        
    
    cpdef emission_prob_ibd_admx_mem_alloc(self, int win_idx):
        cdef int snp_idx 
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef snp_idx_win
        #cdef double sum
        
        snp_idx_win = 0
        self._emission_prob_ibd_admx = <double *********> malloc(self._win_size * sizeof(double ********))
        for snp_idx in range(self.start_snp(win_idx), self.end_snp(win_idx)):
            self._emission_prob_ibd_admx[snp_idx_win] = <double ********> malloc(self.K * sizeof(double*******))
            for admx_idx1 in range(self.K):
                self._emission_prob_ibd_admx[snp_idx_win][admx_idx1] = <double *******> malloc(self.K * sizeof(double******))
                for admx_idx2 in range(admx_idx1+1):
                    self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2] = <double ******> malloc(self.K * sizeof(double*****))
                    for admx_idx3 in range(self.K):
                        self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3] = <double *****> malloc(self.K * sizeof(double****))
                        for admx_idx4 in range(admx_idx3+1):
                            self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4] = <double ****> malloc((self._layer_state_nums[admx_idx1][snp_idx]) * sizeof(double***))
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1] = <double ***> malloc((self._layer_state_nums[admx_idx2][snp_idx]) * sizeof(double**))
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2] = <double **> malloc((self._layer_state_nums[admx_idx3][snp_idx]) * sizeof(double*))
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3] = <double *> malloc((self._layer_state_nums[admx_idx4][snp_idx]) * sizeof(double))
            snp_idx_win += 1
    
#    cpdef emission_prob_ibd_admx_mem_free(self, int win_idx):
#        cdef int snp_idx 
#        cdef int node_idx1
#        cdef int node_idx2
#        cdef int node_idx3
#        cdef int node_idx4
#        cdef int admx_idx1
#        cdef int admx_idx2
#        cdef int admx_idx3
#        cdef int admx_idx4
#        cdef snp_idx_win
#        #cdef double sum
#        
#        snp_idx_win = 0
#        self._emission_prob_ibd_admx = <double *********> malloc(self._win_size * sizeof(double ********))
#        for snp_idx in range(self.start_snp(win_idx), self.end_snp(win_idx)):
#            self._emission_prob_ibd_admx[snp_idx_win] = <double ********> malloc(self.K * sizeof(double*******))
#            for admx_idx1 in range(self.K):
#                self._emission_prob_ibd_admx[snp_idx_win][admx_idx1] = <double *******> malloc(self.K * sizeof(double******))
#                for admx_idx2 in range(admx_idx1+1):
#                    self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2] = <double ******> malloc(self.K * sizeof(double*****))
#                    for admx_idx3 in range(self.K):
#                        self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3] = <double *****> malloc(self.K * sizeof(double****))
#                        for admx_idx4 in range(admx_idx3+1):
#                            self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4] = <double ****> malloc((self._layer_state_nums[admx_idx1][snp_idx]) * sizeof(double***))
#                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
#                                self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1] = <double ***> malloc((self._layer_state_nums[admx_idx2][snp_idx]) * sizeof(double**))
#                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
#                                    self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2] = <double **> malloc((self._layer_state_nums[admx_idx3][snp_idx]) * sizeof(double*))
#                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
#                                        self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3] = <double *> malloc((self._layer_state_nums[admx_idx4][snp_idx]) * sizeof(double))
#            snp_idx_win += 1
    
    cpdef calc_emission_probs_ibd_admx(self, char* chr1, char* chr2, char* chr3, char* chr4, int win_idx):
        
        cdef int snp_idx 
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        cdef int ibd
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef frst_ems
        cdef scnd_emd
        cdef snp_idx_win
        #cdef double sum
        
        snp_idx_win = 0
        for snp_idx in range(self.start_snp(win_idx), self.end_snp(win_idx)):
            #print "calculating emission probs in layer: " + str(snp_idx)
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                            #print "admx_idx1: " + str(admx_idx1) + " admx_idx2: " + str(admx_idx2) + " admx_idx3: " + str(admx_idx3) + " admx_idx4: " + str(admx_idx4) + " snp_idx: " + str(snp_idx) + " node_idx1: " + str(node_idx1) + " node_idx2: " + str(node_idx2) + " node_idx3: " + str(node_idx3) + " node_idx4: " + str(node_idx4)
                                            if chr1[snp_idx] == chr2[snp_idx]:
                                                self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] = \
                                                self._states[admx_idx1][snp_idx][node_idx1].prob_em[int(chr(chr1[snp_idx]))] * \
                                                self._states[admx_idx2][snp_idx][node_idx2].prob_em[int(chr(chr2[snp_idx]))]
                                            else:
                                                self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] = \
                                                self._states[admx_idx1][snp_idx][node_idx1].prob_em[int(chr(chr1[snp_idx]))] * self._states[admx_idx2][snp_idx][node_idx2].prob_em[int(chr(chr2[snp_idx]))] + \
                                                self._states[admx_idx1][snp_idx][node_idx1].prob_em[int(chr(chr2[snp_idx]))] * self._states[admx_idx2][snp_idx][node_idx2].prob_em[int(chr(chr1[snp_idx]))] 
                                                
                                            if chr3[snp_idx] == chr4[snp_idx]:
                                                self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] *= \
                                                self._states[admx_idx3][snp_idx][node_idx3].prob_em[int(chr(chr3[snp_idx]))] * \
                                                self._states[admx_idx4][snp_idx][node_idx4].prob_em[int(chr(chr4[snp_idx]))]
                                            else:
                                                self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] *= \
                                                self._states[admx_idx3][snp_idx][node_idx3].prob_em[int(chr(chr3[snp_idx]))] * self._states[admx_idx4][snp_idx][node_idx4].prob_em[int(chr(chr4[snp_idx]))] + \
                                                self._states[admx_idx3][snp_idx][node_idx3].prob_em[int(chr(chr4[snp_idx]))] * self._states[admx_idx4][snp_idx][node_idx4].prob_em[int(chr(chr3[snp_idx]))]
            snp_idx_win += 1                                                     
    
    cpdef forward_probs_mem_alloc(self, int win_idx):
        cdef int snp_idx 
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        cdef int ibd
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef snp_idx_win
        
        snp_idx_win = 0
        self._forward_probs_ibd_admx = <double **********> malloc(self._win_size * sizeof(double *********))
        for snp_idx in range(self.start_snp(win_idx), self.end_snp(win_idx)):
            self._forward_probs_ibd_admx[snp_idx_win] = <double *********> malloc(self.K * sizeof(double********))
            for admx_idx1 in range(self.K):
                self._forward_probs_ibd_admx[snp_idx_win][admx_idx1] = <double ********> malloc(self.K * sizeof(double*******))
                for admx_idx2 in range(admx_idx1+1):
                    self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2] = <double *******> malloc(self.K * sizeof(double******))
                    for admx_idx3 in range(self.K):
                        self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3] = <double ******> malloc(self.K * sizeof(double*****))
                        for admx_idx4 in range(admx_idx3+1):
                            self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4]= <double *****> malloc((self._layer_state_nums[admx_idx1][snp_idx]) * sizeof(double****))
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1] = <double ****> malloc((self._layer_state_nums[admx_idx2][snp_idx]) * sizeof(double***))
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2] = <double ***> malloc((self._layer_state_nums[admx_idx3][snp_idx]) * sizeof(double**))
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3] = <double **> malloc((self._layer_state_nums[admx_idx4][snp_idx]) * sizeof(double*))
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                            self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] = <double *> malloc(2 * sizeof(double))                        
                                            for ibd in range(2):
                                                self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = 0
            snp_idx_win += 1
            
    cpdef forward_probs_mem_free(self, int win_idx):
        cdef int snp_idx 
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        cdef int ibd
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef snp_idx_win
        
        snp_idx_win = 0
        for snp_idx in range(self.start_snp(win_idx), self.end_snp(win_idx)):
            #print "freeing in snp: " + str(snp_idx)
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                            free(self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4])
                                        free(self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3])
                                    free(self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2])
                                free(self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1])
                            free(self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4])
                        free(self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3])
                    free(self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2])
                free(self._forward_probs_ibd_admx[snp_idx_win][admx_idx1])
            free(self._forward_probs_ibd_admx[snp_idx_win])
            snp_idx_win += 1
        free(self._forward_probs_ibd_admx)
                                                
    cdef forward_probs_init(self, int win_idx):
        cdef int snp_idx 
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        cdef int ibd
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef snp_idx_win
    
        snp_idx_win = 0
        for snp_idx in range(self.start_snp(win_idx), self.end_snp(win_idx)):
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                            for ibd in range(2):
                                                self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = 0
            snp_idx_win += 1
    
    cpdef calc_forward_probs_ibd_admx(self, char* chr1, char* chr2, char* chr3, char* chr4, int win_idx):
        
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
        cdef int ibd
        cdef int prev_ibd
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int prev_admx_idx1
        cdef int prev_admx_idx2
        cdef int prev_admx_idx3
        cdef int prev_admx_idx4
        cdef double sum
        cdef double eps_or_1_eps
        cdef snp_idx_win
        cdef int start_snp = self.start_snp(win_idx)

        # first layer
        for admx_idx1 in range(self.K):
            for admx_idx2 in range(admx_idx1+1):
                for admx_idx3 in range(self.K):
                    for admx_idx4 in range(admx_idx3+1):
                        for node_idx1 in range(self._layer_state_nums[admx_idx1][start_snp]):
                            for node_idx2 in range(self._layer_state_nums[admx_idx2][start_snp]):
                                for node_idx3 in range(self._layer_state_nums[admx_idx3][start_snp]):
                                    for node_idx4 in range(self._layer_state_nums[admx_idx4][start_snp]):
                                        for ibd in range(2):
                                            if ibd == 0:
                                                self._forward_probs_ibd_admx[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = \
                                                self._pi[admx_idx1][win_idx][node_idx1] * self._pi[admx_idx2][win_idx][node_idx2] * self._pi[admx_idx3][win_idx][node_idx3] * self._pi[admx_idx4][win_idx][node_idx4] * \
                                                self._emission_prob_ibd_admx[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4]
                                            else:
                                                self._forward_probs_ibd_admx[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = \
                                                min(self._pi[admx_idx1][win_idx][node_idx1],self._pi[admx_idx3][win_idx][node_idx3]) * self._pi[admx_idx2][win_idx][node_idx2] * self._pi[admx_idx4][win_idx][node_idx4] * \
                                                self._emission_prob_ibd_admx[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4]
        
        # all other layers
        snp_idx_win = 0
        for snp_idx in range(self.start_snp(win_idx), self.end_snp(win_idx)-1):
            #print "calculating forward probs in layer: " + str(snp_idx)
            
            sum = 0              
            # calculate forward probabilities
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx+1]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx+1]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx+1]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx+1]):
                                            for ibd in range(2):
                                                
                                                for prev_node_idx1 in range(self._states[admx_idx1][snp_idx+1][node_idx1].in_trans_num):
                                                    for prev_node_idx2 in range(self._states[admx_idx2][snp_idx+1][node_idx2].in_trans_num):
                                                        for prev_node_idx3 in range(self._states[admx_idx3][snp_idx+1][node_idx3].in_trans_num):
                                                            for prev_node_idx4 in range(self._states[admx_idx4][snp_idx+1][node_idx4].in_trans_num):
                                                                #for prev_ibd in range(2):
                                                                    prev_ibd = ibd
                                                                    prev_node1 = self._back_trans_idx[admx_idx1][snp_idx+1][node_idx1][prev_node_idx1]
                                                                    prev_node2 = self._back_trans_idx[admx_idx2][snp_idx+1][node_idx2][prev_node_idx2]
                                                                    prev_node3 = self._back_trans_idx[admx_idx3][snp_idx+1][node_idx3][prev_node_idx3]
                                                                    prev_node4 = self._back_trans_idx[admx_idx4][snp_idx+1][node_idx4][prev_node_idx4]
                                                                    if ibd == 0:
                                                                        a=1
                                                                        
                                                                        self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] += \
                                                                        self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][prev_node1][prev_node2][prev_node3][prev_node4][prev_ibd] * \
                                                                        self._back_trans[admx_idx1][snp_idx+1][node_idx1][prev_node_idx1] * \
                                                                        self._back_trans[admx_idx2][snp_idx+1][node_idx2][prev_node_idx2] * \
                                                                        self._back_trans[admx_idx3][snp_idx+1][node_idx3][prev_node_idx3] * \
                                                                        self._back_trans[admx_idx4][snp_idx+1][node_idx4][prev_node_idx4] * \
                                                                        self._emission_prob_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] #* \
                                                                        #self._s[snp_idx+1][prev_ibd][ibd]
                                                                        
                                                                    else:
                                                                        if chr1[snp_idx+1] == chr3[snp_idx+1]:
                                                                            eps_or_1_eps = 1 - eps
                                                                        else:
                                                                            eps_or_1_eps = eps
                                                                        
                                                                        self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] += \
                                                                        self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][prev_node1][prev_node2][prev_node3][prev_node4][prev_ibd] * \
                                                                        min(self._back_trans[admx_idx1][snp_idx+1][node_idx1][prev_node_idx1], 
                                                                            self._back_trans[admx_idx3][snp_idx+1][node_idx3][prev_node_idx3]) * \
                                                                        self._back_trans[admx_idx2][snp_idx+1][node_idx2][prev_node_idx2] * \
                                                                        self._back_trans[admx_idx4][snp_idx+1][node_idx4][prev_node_idx4] * \
                                                                        self._emission_prob_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] * eps_or_1_eps
                                                                        #self._s[snp_idx+1][prev_ibd][ibd] * \                 
                                                sum += self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd]        
            
            # rescaling to avoid underflow
            if sum == 0:
                sum = 1
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx+1]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx+1]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx+1]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx+1]):
                                            for ibd in range(2):
                                                self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = \
                                                self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] * (1/sum)
                                                
                                                #if snp_idx == (self.end_snp(win_idx) - 1):
                                                #    self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] += \
                                                #    self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd]
            snp_idx_win += 1
    
    cpdef backward_probs_mem_alloc(self, int win_idx):
        cdef int snp_idx 
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        cdef int ibd
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef snp_idx_win
        
        snp_idx_win = 0
        self._backward_probs_ibd_admx = <double **********> malloc(self._win_size * sizeof(double *********))
        for snp_idx in range(self.start_snp(win_idx), self.end_snp(win_idx)):
            self._backward_probs_ibd_admx[snp_idx_win] = <double *********> malloc(self.K * sizeof(double********))
            for admx_idx1 in range(self.K):
                self._backward_probs_ibd_admx[snp_idx_win][admx_idx1] = <double ********> malloc(self.K * sizeof(double*******))
                for admx_idx2 in range(admx_idx1+1):
                    self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2] = <double *******> malloc(self.K * sizeof(double******))
                    for admx_idx3 in range(self.K):
                        self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3] = <double ******> malloc(self.K * sizeof(double*****))
                        for admx_idx4 in range(admx_idx3+1):
                            self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4]= <double *****> malloc((self._layer_state_nums[admx_idx1][snp_idx]) * sizeof(double****))
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1] = <double ****> malloc((self._layer_state_nums[admx_idx2][snp_idx]) * sizeof(double***))
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2] = <double ***> malloc((self._layer_state_nums[admx_idx3][snp_idx]) * sizeof(double**))
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3] = <double **> malloc((self._layer_state_nums[admx_idx4][snp_idx]) * sizeof(double*))
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                            self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] = <double *> malloc(2 * sizeof(double))                        
                                            for ibd in range(2):
                                                self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = 0
            snp_idx_win += 1 
    
    cdef backward_probs_init(self, int win_idx):
        cdef int snp_idx 
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        cdef int ibd
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef snp_idx_win
        
        snp_idx_win = 0
        for snp_idx in range(self.start_snp(win_idx), self.end_snp(win_idx)):
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                            for ibd in range(2):
                                                self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = 0
            snp_idx_win += 1
                             
    cpdef calc_backward_probs_ibd_admx(self, char* chr1, char* chr2, char* chr3, char* chr4, int win_idx):
        
        cdef int snp_idx 
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
        cdef int ibd
        cdef int nxt_ibd
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int nxt_admx_idx1
        cdef int nxt_admx_idx2
        cdef int nxt_admx_idx3
        cdef int nxt_admx_idx4
        cdef double sum
        cdef double eps_or_1_eps
        cdef snp_idx_win
        
        # last layer
        for admx_idx1 in range(self.K):
            for admx_idx2 in range(admx_idx1+1):
                for admx_idx3 in range(self.K):
                    for admx_idx4 in range(admx_idx3+1):
                        for node_idx1 in range(self._layer_state_nums[admx_idx1][self.end_snp(win_idx) - 1]):
                            for node_idx2 in range(self._layer_state_nums[admx_idx2][self.end_snp(win_idx) - 1]):
                                for node_idx3 in range(self._layer_state_nums[admx_idx3][self.end_snp(win_idx) - 1]):
                                    for node_idx4 in range(self._layer_state_nums[admx_idx4][self.end_snp(win_idx) - 1]):
                                        for ibd in range(2):
                                            self._backward_probs_ibd_admx[self._win_size - 1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = 1
        
        # all other layers
        snp_idx_win = self._win_size - 2 
        for snp_idx in reversed(range(self.start_snp(win_idx), self.end_snp(win_idx) - 1)):
            #print "calculating backward probs in layer: " + str(snp_idx)
            
            sum = 0              
            # calculate forward probabilities
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                            for ibd in range(2):
                                                
                                                for nxt_node_idx1 in range(self._states[admx_idx1][snp_idx][node_idx1].out_trans_num):
                                                    for nxt_node_idx2 in range(self._states[admx_idx2][snp_idx][node_idx2].out_trans_num):
                                                        for nxt_node_idx3 in range(self._states[admx_idx3][snp_idx][node_idx3].out_trans_num):
                                                            for nxt_node_idx4 in range(self._states[admx_idx4][snp_idx][node_idx4].out_trans_num):
                                                                #for nxt_ibd in range(2):
                                                                    nxt_ibd = ibd
                                                                    nxt_node1 = self._trans_idx[admx_idx1][snp_idx][node_idx1][nxt_node_idx1]
                                                                    nxt_node2 = self._trans_idx[admx_idx2][snp_idx][node_idx2][nxt_node_idx2]
                                                                    nxt_node3 = self._trans_idx[admx_idx3][snp_idx][node_idx3][nxt_node_idx3]
                                                                    nxt_node4 = self._trans_idx[admx_idx4][snp_idx][node_idx4][nxt_node_idx4]
                                                                    if ibd == 0:
                                                                        a=1            
                                                                        self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] += \
                                                                        self._backward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_node1][nxt_node2][nxt_node3][nxt_node4][nxt_ibd] * \
                                                                        self._trans[admx_idx1][snp_idx][node_idx1][nxt_node_idx1] * \
                                                                        self._trans[admx_idx2][snp_idx][node_idx2][nxt_node_idx2] * \
                                                                        self._trans[admx_idx3][snp_idx][node_idx3][nxt_node_idx3] * \
                                                                        self._trans[admx_idx4][snp_idx][node_idx4][nxt_node_idx4] * \
                                                                        self._emission_prob_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_node1][nxt_node2][nxt_node3][nxt_node4] #* \
                                                                        #self._s[snp_idx][ibd][nxt_ibd] 
                                                                    else:
                                                                        if chr1[snp_idx+1] == chr3[snp_idx+1]:
                                                                            eps_or_1_eps = 1 - eps
                                                                        else:
                                                                            eps_or_1_eps = eps                                
                                                                        self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] += \
                                                                        self._backward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_node1][nxt_node2][nxt_node3][nxt_node4][nxt_ibd] * \
                                                                        min(self._trans[admx_idx1][snp_idx][node_idx1][nxt_node_idx1], 
                                                                            self._trans[admx_idx3][snp_idx][node_idx3][nxt_node_idx3]) * \
                                                                        self._trans[admx_idx2][snp_idx][node_idx2][nxt_node_idx2] * \
                                                                        self._trans[admx_idx4][snp_idx][node_idx4][nxt_node_idx4] * \
                                                                        self._emission_prob_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_node1][nxt_node2][nxt_node3][nxt_node4] * eps_or_1_eps
                                                                        #self._s[snp_idx][ibd][nxt_ibd] * \
                                                sum += self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd]        
            
            # rescaling to avoid underflow
            if sum == 0:
                sum = 1
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                            for ibd in range(2):
                                                self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = \
                                                self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] * (1/sum)
                                                #print " snp_idx: " + str(snp_idx) + " admx_idx1: " + str(admx_idx1) + " admx_idx2: " + str(admx_idx2) + " admx_idx3: " + str(admx_idx3) + " admx_idx4: " + str(admx_idx4) + \
                                                #" node_idx1: " + str(node_idx1) + " node_idx2: " + str(node_idx2) + " node_idx3: " + str(node_idx3) + " node_idx4: " + str(node_idx4) + " ibd: " + str(ibd),
                                                #print " backward prob: " + str(self._backward_probs_ibd_admx[snp_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd])
            snp_idx_win -= 1                                            
    
    cpdef top_level_alloc_mem(self):
        cdef int win_idx
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int ibd
         
        self._top_level_ems_prob = <double ******> malloc(self.get_num_windows() * sizeof(double *****))
        self._top_level_forward_probs = <double ******> malloc(self.get_num_windows() * sizeof(double *****))
        self._top_level_backward_probs = <double ******> malloc(self.get_num_windows() * sizeof(double *****))
        for win_idx in range(self.get_num_windows()):
            self._top_level_ems_prob[win_idx] = <double *****> malloc(self.K * sizeof(double ****))
            self._top_level_forward_probs[win_idx] = <double *****> malloc(self.K * sizeof(double ****))
            self._top_level_backward_probs[win_idx] = <double *****> malloc(self.K * sizeof(double ****))
            for admx_idx1 in range(self.K):
                self._top_level_ems_prob[win_idx][admx_idx1] = <double ****> malloc((admx_idx1+1) * sizeof(double ***))
                self._top_level_forward_probs[win_idx][admx_idx1] = <double ****> malloc((admx_idx1+1) * sizeof(double ***))
                self._top_level_backward_probs[win_idx][admx_idx1] = <double ****> malloc((admx_idx1+1) * sizeof(double ***))
                for admx_idx2 in range(admx_idx1+1):
                    self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2] = <double ***> malloc(self.K * sizeof(double **))
                    self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2] = <double ***> malloc(self.K * sizeof(double **))
                    self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2] = <double ***> malloc(self.K * sizeof(double **))
                    for admx_idx3 in range(self.K):
                        self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3] = <double **> malloc((admx_idx3+1) * sizeof(double *))
                        self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3] = <double **> malloc((admx_idx3+1) * sizeof(double *))
                        self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3] = <double **> malloc((admx_idx3+1) * sizeof(double *))
                        for admx_idx4 in range(admx_idx3+1):
                            self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4] = <double *> malloc(2 * sizeof(double))
                            self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4] = <double *> malloc(2 * sizeof(double))
                            self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4] = <double *> malloc(2 * sizeof(double))
                            for ibd in range(2):
                                self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
                                self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
                                self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
    
    cpdef top_level_init(self):
        cdef int win_idx
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int ibd
         
        for win_idx in range(self.get_num_windows()):
            #print "initializing window %d" % win_idx
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for ibd in range(2):
                                self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
                                self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
                                self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0         
    
    cpdef posterior_decoding_ibd_admx(self, int win_idx):
        cdef int max_node_idx1
        cdef int max_node_idx2
        cdef int max_node_idx3
        cdef int max_node_idx4
        cdef int max_admx_idx1
        cdef int max_admx_idx2
        cdef int max_admx_idx3
        cdef int max_admx_idx4
        cdef int max_ibd
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int ibd
        cdef double max_gamma
        cdef double curr_gamma
        cdef snp_idx_win
        a1 = ""
        a2 = ""
        a3 = ""
        a4 = ""
        d1 = ""
        d2 = ""
        d3 = ""
        d4 = ""  
        i = ""
        snp_idx_win = 0
        for snp_idx in range(self.start_snp(win_idx), self.end_snp(win_idx)):
            #print "calculating posterior decoding in layer: " + str(snp_idx)
            max_gamma = -DBL_MAX
            max_node_idx1 = -1
            max_node_idx2 = -1
            max_node_idx3 = -1
            max_node_idx4 = -1
            max_admx_idx1 = -1
            max_admx_idx2 = -1
            max_admx_idx3 = -1
            max_admx_idx4 = -1
            max_ibd = -1
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                            for ibd in range(2):
                                                curr_gamma = self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] * \
                                                self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd]
                                                if curr_gamma  > max_gamma:
                                                    max_gamma = curr_gamma
                                                    max_node_idx1 = node_idx1
                                                    max_node_idx2 = node_idx2
                                                    max_node_idx3 = node_idx3
                                                    max_node_idx4 = node_idx4
                                                    max_admx_idx1 = admx_idx1
                                                    max_admx_idx2 = admx_idx2
                                                    max_admx_idx3 = admx_idx3
                                                    max_admx_idx4 = admx_idx4
                                                    max_ibd = ibd
            a1 += str(max_admx_idx1)
            a2 += str(max_admx_idx2)
            a3 += str(max_admx_idx3)
            a4 += str(max_admx_idx4)
            d1 += chr(get_likely_allele(self._states[max_admx_idx1][snp_idx][max_node_idx1]))
            d2 += chr(get_likely_allele(self._states[max_admx_idx2][snp_idx][max_node_idx2]))
            d3 += chr(get_likely_allele(self._states[max_admx_idx3][snp_idx][max_node_idx3]))
            d4 += chr(get_likely_allele(self._states[max_admx_idx4][snp_idx][max_node_idx4]))
            i += str(max_ibd)
            snp_idx_win += 1
        return (a1,a2,a3,a4,d1,d2,d3,d4,i)
    
    cpdef calc_top_level_ems_probs_inner(self, char* chr1, char* chr2, char* chr3, char* chr4):
        cdef int win_idx
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int ibd
        cdef last_snp_win
        for win_idx in range(self.get_num_windows()):
            print "calculating top level ems probs for window: " + str(win_idx)
            last_snp_win = self.end_snp(win_idx) - self.start_snp(win_idx)
            self.emission_prob_ibd_admx_mem_alloc(win_idx)
            self.calc_emission_probs_ibd_admx(chr1,chr2,chr3,chr4,win_idx)
            self.forward_probs_mem_alloc(win_idx)
            self.calc_forward_probs_ibd_admx(chr1,chr2,chr3,chr4,win_idx)
            #self.backward_probs_mem_alloc(win_idx)
            #self.calc_backward_probs_ibd_admx(chr1,chr2,chr3,chr4,win_idx)
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for ibd in range(2):
                                sum = 0
                                for node_idx1 in range(self._layer_state_nums[admx_idx1][(win_idx*self._win_size)+last_snp_win-1]):
                                    for node_idx2 in range(self._layer_state_nums[admx_idx2][(win_idx*self._win_size)+last_snp_win-1]):
                                        for node_idx3 in range(self._layer_state_nums[admx_idx3][(win_idx*self._win_size)+last_snp_win-1]):
                                            for node_idx4 in range(self._layer_state_nums[admx_idx4][(win_idx*self._win_size)+last_snp_win-1]):
                                                self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] += \
                                                self._forward_probs_ibd_admx[last_snp_win-1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd]
                                                #print ">top level ems probs (win_idx admx_idx1 admx_idx2 admx_idx3 admx_idx4 ibd node_idx1 node_idx2 node_idx3 node_idx4)" + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + " " + str(node_idx1) + " " +  str(node_idx2) + " " +  str(node_idx3) + " " +  str(node_idx4) + " " + ": " + str(self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]) + " ibd_admx_forward_probs " + str(self._forward_probs_ibd_admx[last_snp_win-1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd])
            
            self.forward_probs_mem_free(win_idx)
    
    cpdef calc_top_level_ems_probs(self, int hap_idx1, int hap_idx2, int hap_idx3, int hap_idx4):
        self.calc_top_level_ems_probs_inner(self._haplos[hap_idx1],self._haplos[hap_idx2],self._haplos[hap_idx3],self._haplos[hap_idx4])
    
    cpdef calc_top_level_forward_probs(self):
        cdef int win_idx
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int ibd
        cdef int prev_admx_idx1
        cdef int prev_admx_idx2
        cdef int prev_admx_idx3
        cdef int prev_admx_idx4
        cdef int prev_ibd
        
        for admx_idx1 in range(self.K):
            for admx_idx2 in range(admx_idx1+1):
                for admx_idx3 in range(self.K):
                    for admx_idx4 in range(admx_idx3+1):
                        for ibd in range(2):
                            self._top_level_forward_probs[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] += \
                            self._top_level_ems_prob[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] * \
                            self._alphas[admx_idx1] * self._alphas[admx_idx2] * self._alphas[admx_idx3] * self._alphas[admx_idx4] 
                            #print "first window ems probs " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_ems_prob[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                            #print "params: " + str(self._ibd_prior[ibd]) + " " + str(self._alphas[admx_idx1]) + " " + str(self._alphas[admx_idx2]) + " " + str(self._alphas[admx_idx3]) + " "  + str(self._alphas[admx_idx4]) 
                            #print "first window forward probs " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_forward_probs[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
        
        for win_idx in range(self.get_num_windows()-1):
            print "calculating top level forward probs for window: " + str(win_idx)
            sum = 0
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for ibd in range(2):
                                for prev_admx_idx1 in range(self.K):
                                    for prev_admx_idx2 in range(prev_admx_idx1+1):
                                        for prev_admx_idx3 in range(self.K):
                                            for prev_admx_idx4 in range(prev_admx_idx3+1):
                                                for prev_ibd in range(2):
                                                    self._top_level_forward_probs[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] += \
                                                    self._top_level_forward_probs[win_idx][prev_admx_idx1][prev_admx_idx2][prev_admx_idx3][prev_admx_idx4][prev_ibd] * \
                                                    self._top_level_ems_prob[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] * \
                                                    self._s[win_idx][prev_ibd][ibd] * \
                                                    self._anc_trans[win_idx][prev_admx_idx1][prev_admx_idx2][prev_admx_idx3][prev_admx_idx4][admx_idx1][admx_idx2][admx_idx3][admx_idx4] * \
                                                    self._alphas[admx_idx1] * self._alphas[admx_idx2] * self._alphas[admx_idx3] * self._alphas[admx_idx4]
                                sum += self._top_level_forward_probs[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]
                                if False: #ibd == 0:
                                    print str(win_idx) + " top level ems probs " + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_ems_prob[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                                    print str(win_idx) + " top level ibd prior and alphas: " + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._ibd_prior[0]) + " " + str(self._alphas[admx_idx1]) + " " + str(self._alphas[admx_idx2]) + " " + str(self._alphas[admx_idx3]) + " "  + str(self._alphas[admx_idx4])
                                    print str(win_idx) + " top level ibd trans: " + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._s[win_idx][0][0]) + " " + str(self._s[win_idx][0][1]) + " " + str(self._s[win_idx][1][0]) + " " + str(self._s[win_idx][1][1])
                                    for prev_admx_idx1 in range(self.K):
                                        for prev_admx_idx2 in range(prev_admx_idx1+1):
                                            for prev_admx_idx3 in range(self.K):
                                                for prev_admx_idx4 in range(prev_admx_idx3+1):
                                                    print str(win_idx) + " top level anc trans: " + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + ": " + str(self._anc_trans[win_idx][prev_admx_idx1][prev_admx_idx2][prev_admx_idx3][prev_admx_idx4][admx_idx1][admx_idx2][admx_idx3][admx_idx4])
                                    print str(win_idx) + " top level forward probs " + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_forward_probs[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
            
            scale = 1.0                                        
            if sum == 0:
                scale = 1.0
            else:
                scale = 1.0/sum
                if not c_isfinite(scale):
                    scale = 1.0 
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for ibd in range(2):
                                self._top_level_forward_probs[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = \
                                self._top_level_forward_probs[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] * scale
                                if not c_isfinite(self._top_level_forward_probs[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]):
                                    self._top_level_forward_probs[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0.0 
                                
    
    cpdef calc_top_level_backward_probs(self):
        cdef int win_idx
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int ibd
        cdef int nxt_admx_idx1
        cdef int nxt_admx_idx2
        cdef int nxt_admx_idx3
        cdef int nxt_admx_idx4
        cdef int nxt_ibd
         
        for admx_idx1 in range(self.K):
            for admx_idx2 in range(admx_idx1+1):
                for admx_idx3 in range(self.K):
                    for admx_idx4 in range(admx_idx3+1):
                        for ibd in range(2):
                            self._top_level_backward_probs[self.get_num_windows()-1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 1
         
        for win_idx in reversed(range(self.get_num_windows()-1)):
            print "calculating top level forward probs for window: " + str(win_idx)
            sum = 0
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for ibd in range(2):
                                for nxt_admx_idx1 in range(self.K):
                                    for nxt_admx_idx2 in range(nxt_admx_idx1+1):
                                        for nxt_admx_idx3 in range(self.K):
                                            for nxt_admx_idx4 in range(nxt_admx_idx3+1):
                                                for nxt_ibd in range(2):
                                                    self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] += \
                                                    self._top_level_backward_probs[win_idx+1][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4][nxt_ibd] * \
                                                    self._top_level_ems_prob[win_idx+1][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4][nxt_ibd] * \
                                                    self._s[win_idx][ibd][nxt_ibd] * \
                                                    self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4] * \
                                                    self._alphas[nxt_admx_idx1] * self._alphas[nxt_admx_idx2] * self._alphas[nxt_admx_idx3] * self._alphas[nxt_admx_idx4]
                                                    #print ">backward " + str(win_idx) + " " + str(admx_idx1) + " " + str(admx_idx2) + " " + str(admx_idx3) + " " + str(admx_idx4) + " " + str(ibd) + " " + str(nxt_admx_idx1) + " " + str(nxt_admx_idx2) + " " + str(nxt_admx_idx3) + " " + str(nxt_admx_idx4) + " " + str(nxt_ibd) + " ems_prob: " + str(self._top_level_ems_prob[win_idx+1][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4][nxt_ibd])
                                                    #print ">backward " + str(win_idx) + " " + str(admx_idx1) + " " + str(admx_idx2) + " " + str(admx_idx3) + " " + str(admx_idx4) + " " + str(ibd) + " " + str(nxt_admx_idx1) + " " + str(nxt_admx_idx2) + " " + str(nxt_admx_idx3) + " " + str(nxt_admx_idx4) + " " + str(nxt_ibd) + " _s: " + str(self._s[win_idx][ibd][nxt_ibd])
                                                    #print ">backward " + str(win_idx) + " " + str(admx_idx1) + " " + str(admx_idx2) + " " + str(admx_idx3) + " " + str(admx_idx4) + " " + str(ibd) + " " + str(nxt_admx_idx1) + " " + str(nxt_admx_idx2) + " " + str(nxt_admx_idx3) + " " + str(nxt_admx_idx4) + " " + str(nxt_ibd) + " _anc_trans: " + str(self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4])
                                                    #print ">backward " + str(win_idx) + " " + str(admx_idx1) + " " + str(admx_idx2) + " " + str(admx_idx3) + " " + str(admx_idx4) + " " + str(ibd) + " " + str(nxt_admx_idx1) + " " + str(nxt_admx_idx2) + " " + str(nxt_admx_idx3) + " " + str(nxt_admx_idx4) + " " + str(nxt_ibd) + " _backward prob i+1: " + str(self._top_level_backward_probs[win_idx+1][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4][nxt_ibd])
                                                    #print "yyyyyyyyyyyyyyyyyyyy " + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + str(self._top_level_backward_probs[win_idx+1][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4][nxt_ibd])
                                sum += self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]
                                #if ibd == 0:
                                #    print "back top level params: " + str(self._ibd_prior[ibd]) + " " + str(self._alphas[admx_idx1]) + " " + str(self._alphas[admx_idx2]) + " " + str(self._alphas[admx_idx3]) + " "  + str(self._alphas[admx_idx4])
                                #    print "back top level forward probs " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                                                    
            scale = 1.0                                        
            if sum == 0:
                scale = 1.0
            else:
                scale = 1.0/sum
                if not c_isfinite(scale):
                    scale = 1.0 
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for ibd in range(2):
                                self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = \
                                self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] * scale
                                if not c_isfinite(self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]):
                                    self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0.0
                                                    
    cpdef posterior_top_level_decoding(self):
        cdef int win_idx
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int ibd
        cdef int max_admx_idx1
        cdef int max_admx_idx2
        cdef int max_admx_idx3
        cdef int max_admx_idx4
        cdef int max_ibd
        cdef double max_gamma
        cdef double curr_gamma
        a1 = ""
        a2 = ""
        a3 = ""
        a4 = ""
        i = ""
        
        for win_idx in range(self.get_num_windows()):
            max_admx_idx1 = -1
            max_admx_idx2 = -1
            max_admx_idx3 = -1
            max_admx_idx4 = -1
            max_ibd = -1
            max_gamma = -DBL_MAX
            print "top level decoding in window %d" % win_idx
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(admx_idx1+1):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(admx_idx3+1):
                            for ibd in range(2):
                                curr_gamma = self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] * \
                                self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]
                                print str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]) + " " + str(self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]) + " " + str(curr_gamma)
                                if curr_gamma  > max_gamma:
                                    max_gamma = curr_gamma
                                    max_admx_idx1 = admx_idx1
                                    max_admx_idx2 = admx_idx2
                                    max_admx_idx3 = admx_idx3
                                    max_admx_idx4 = admx_idx4
                                    max_ibd = ibd
            a1 += str(max_admx_idx1)
            a2 += str(max_admx_idx2)
            a3 += str(max_admx_idx3)
            a4 += str(max_admx_idx4)
            i += str(max_ibd)
            
            a1_filt = ''
            a2_filt = ''
            a3_filt = ''
            a4_filt = ''
            i_filt = ''
            for ind in range(len(a1)):
                start = max(0,ind-3)
                end = min(ind+4,len(a1))
                a1_filt += str(int(np.median([int(x) for x in a1[start:end]])))
                a2_filt += str(int(np.median([int(x) for x in a2[start:end]])))
                a3_filt += str(int(np.median([int(x) for x in a3[start:end]])))
                a4_filt += str(int(np.median([int(x) for x in a4[start:end]])))
                i_filt += str(int(np.median([int(x) for x in i[start:end]])))
                
            pairIBD=cPairIBD()
            for win_i in range(len(i_filt)):
                if i[win_i] == '1':
                    pairIBD.add_interval(win_i*self._win_size,(win_i+1)*self._win_size)
            pairIBD.merge_intervals()
        
        return (a1_filt,a2_filt,a3_filt,a4_filt,i_filt,pairIBD)
    
    cdef int start_snp(self, int win_idx):
        return win_idx * self._win_size
    
    cdef int end_snp(self, int win_idx):
        return min((win_idx + 1) * self._win_size, self._snp_num)
    
    cpdef int get_num_windows(self):
        cdef int num_win
        num_win = int(self._snp_num / self._win_size)
        if self._snp_num % self._win_size > 0:
            num_win += 1
        return num_win
        
        
        
        
