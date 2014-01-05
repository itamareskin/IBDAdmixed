#cython: profile=True

from __future__ import division
import random
import numpy as np
cimport numpy as np
import os
import math
from libc.stdlib cimport malloc, free
from itertools import islice, combinations_with_replacement
from libc.math cimport exp, log
from libc.float cimport DBL_MIN, DBL_MAX
from libc.limits cimport ULONG_MAX,LONG_MIN
from libcpp cimport bool
from IBD.cIBD import cPairIBD
from intersection import Interval, IntervalTree
from sys import stdout
import scipy.misc
 
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
     
    def __cinit__(self, map_file_name, char* log_dir, char* log_prefix, int k=1, int g=8, int win_size=25, int max_snp_num=1000000000, double eps = 1e-4, double min_score = 0, phased=True, bool debug=False):
        
        # total number of SNPs to be analyzed
        with open(map_file_name) as map_file:
            data = map_file.readlines()
            self._snp_num = len(data) if len(data) <= max_snp_num else max_snp_num
        
        self.eps = eps
        
        self._min_score = min_score
        
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
        
        self._log_dir = log_dir
        
        self._log_prefix = log_prefix
        
        self._debug = debug
        
        self._phased = phased
        
        # rate of transition from IBD to No-IBD and vice-versa (notation from the Browning paper)
        self._t_0_1 = <double *> malloc(self.K * sizeof(double))
        self._t_1_0 = <double *> malloc(self.K * sizeof(double))
        for i in range(self.K):
            self._t_0_1[i] = 1e-5
            self._t_1_0[i] = 1
        
        # prior probability of IBD
        self._ibd_prior = <double **> malloc(self.K * sizeof(double*))
        for i in range(self.K):
            self._ibd_prior[i] = <double *> malloc(2 * sizeof(double))
            self._ibd_prior[i][1] = self._t_0_1[i] / (self._t_0_1[i] + self._t_1_0[i])
            self._ibd_prior[i][0] = 1-self._ibd_prior[i][1]
        
        # proportions of ancestry from each ancestral populations
        self._alphas = <double *> malloc(self.K * sizeof(double))
        for i in range(self.K):
            self._alphas[i] = (1.0/self.K)
        
        # allocate memory
        self._snp_ids = <char **> malloc(self._snp_num * sizeof(char *))
        self._position = <int *> malloc(self._snp_num * sizeof(int))
        self._genetic_dist = <double *> malloc(self._snp_num * sizeof(double))
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
            
        #self._genetic_map = <gen_map_entry *> malloc(self._snp_num * sizeof(gen_map_entry))
        #print "allocated mem"
        #for snp_idx in range(self._snp_num):
        #    self._genetic_map[snp_idx] = create_gen_map_entry(snp_idx,1.63936,0.0010840)
            
        if self._debug:
            print "DEBUG!"
            self._inner_probs_file = open(self._log_dir+"/"+self._log_prefix+".inner.probs.txt","w")
            self._probs_file = open(self._log_dir+"/"+self._log_prefix+".probs.txt","w")
            self._trans_file = open(self._log_dir+"/"+self._log_prefix+".trans.txt","w")
            self._ems_file = open(self._log_dir+"/"+self._log_prefix+".ems.txt","w") 
            self._ibs_file = open(self._log_dir+"/"+self._log_prefix+".ibs.txt","w")
            self._breakpoints_file = open(self._log_dir+"/"+self._log_prefix+".breakpoints.txt","w")
        
        self._ibs = <bool *> malloc(self.get_num_windows() * sizeof(bool))
        self._ibd_segment_start = <bool *> malloc(self.get_num_windows() * sizeof(bool))
        self._ibd_segment_end = <bool *> malloc(self.get_num_windows() * sizeof(bool))
        for win_idx in range(self.get_num_windows()):
            self._ibs[win_idx] = True
            self._ibd_segment_start[win_idx] = False
            self._ibd_segment_end[win_idx] = False
            
        self._exact_ibd_starts = <int *> malloc(self.get_num_windows() * sizeof(int))
        self._exact_ibd_ends = <int *> malloc(self.get_num_windows() * sizeof(int))
        
        self._ibd_probs = <double *> malloc(self.get_num_windows() * sizeof(double))
        self._no_ibd_probs = <double *> malloc(self.get_num_windows() * sizeof(double))
        self._lod_scores = <double *> malloc(self.get_num_windows() * sizeof(double))
            
        self._chr1 = <bool *> malloc(self._snp_num * sizeof(bool))
        self._chr2 = <bool *> malloc(self._snp_num * sizeof(bool))
        self._chr3 = <bool *> malloc(self._snp_num * sizeof(bool))
        self._chr4 = <bool *> malloc(self._snp_num * sizeof(bool))
        
        self.read_map(map_file_name)
                
    def __dealloc__(self):
        if self._debug:
            self._inner_probs_file.close()
            self._probs_file.close() 
            self._trans_file.close()
            self._ems_file.close()
            self._ibs_file.close()
            self._breakpoints_file.close()
    
    cpdef set_ibd_trans_rate(self, anc, t_0_1, t_1_0):
        self._t_0_1[anc] = t_0_1
        self._t_1_0[anc] = t_1_0
        
        self._ibd_prior[anc][1] = self._t_0_1[anc] / (self._t_0_1[anc] + self._t_1_0[anc])
        self._ibd_prior[anc][0] = 1-self._ibd_prior[anc][1]
            
    cpdef set_alphas(self, alphas):
        if len(alphas) != self.K:
            print "wrong number of alphas"
            return 
        for i in range(self.K):
            self._alphas[i] = alphas[i]
            
    cpdef set_chrs(self, char* chr1, char* chr2, char* chr3, char* chr4):
    
        for snp_idx in range(self._snp_num):
            if int(chr(chr1[snp_idx])) == 0:
                self._chr1[snp_idx] = 0
            else:
                if int(chr(chr1[snp_idx])) == 1:
                    self._chr1[snp_idx] = 1
                else:
                    print "unidentified allele: " + str(int(chr(chr1[snp_idx])))
                    exit(-1)
            if int(chr(chr2[snp_idx])) == 0:
                self._chr2[snp_idx] = 0
            else:
                if int(chr(chr2[snp_idx])) == 1:
                    self._chr2[snp_idx] = 1
                else:
                    print "unidentified allele: " + str(int(chr(chr2[snp_idx])))
                    exit(-1)
            if int(chr(chr3[snp_idx])) == 0:
                self._chr3[snp_idx] = 0
            else:
                if int(chr(chr3[snp_idx])) == 1:
                    self._chr3[snp_idx] = 1
                else:
                    print "unidentified allele: " + str(int(chr(chr3[snp_idx])))
                    exit(-1)
            if int(chr(chr4[snp_idx])) == 0:
                self._chr4[snp_idx] = 0
            else:
                if int(chr(chr4[snp_idx])) == 1:
                    self._chr4[snp_idx] = 1
                else:
                    print "unidentified allele: " + str(int(chr(chr4[snp_idx])))
                    exit(-1)
            
    
    def set_prefix_string(self, char* new_prefix_string):
        self._prefix_string = <char *> malloc((strlen(new_prefix_string) + 1) * sizeof(char))
        strncpy(self._prefix_string, new_prefix_string, strlen(new_prefix_string)) 
        self._prefix_string[strlen(new_prefix_string)] = '\0'
        
    
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
        
#     def read_genetic_map(self, genetic_map_file_name):
#         '''
#         read the genetic map from hapmap format file
#         '''
#         if not os.path.exists(genetic_map_file_name):
#             print "the file: " + genetic_map_file_name + " does not exist!"
#             exit(-1)
#         
#         print "reading from genetic map file: " + genetic_map_file_name
#         
#         with open(genetic_map_file_name) as genetic_map_file:
#             
#             first_read = True
#             done = False
#             buffer_size = 100
#             snp_idx = 0
#             
#             while True:
#                 
#                 if done:
#                     break
#                 # read next buffer_size lines from the file
#                 if not first_read:
#                     lines = list(islice(genetic_map_file, buffer_size))
#                 else:
#                     lines = list(islice(genetic_map_file, 1, buffer_size))
#                     first_read = False
#                 
#                 if len(lines) == 0:
#                     done = True
#                 
#                 for line in lines:
#                     #print "line: " + line
#                     line = line.split(" ")
#                     #print "pos: " + line[0]
#                     # TODO: this is a bug - need to read only positions that are in _position
#                     self._genetic_map[snp_idx] = create_gen_map_entry(int(line[0]),float(line[1]),float(line[2]))
#                     snp_idx += 1
#                     if snp_idx >= self._snp_num:
#                         done = True;
#                         break
        
    def generate_random_haps_inplace(self, int anc, int count):
        self._nr_haplos = count
        self._haplos = <char **> malloc(self._nr_haplos * sizeof(char *))
        for hap_idx in range(4):    
            self._haplos[hap_idx] = <char *> malloc((self._snp_num) * sizeof(char))
            self.generate_random_hap(anc)
            hap = self.generate_random_hap(anc)
            strncpy(self._haplos[hap_idx], hap, self._snp_num)
    
    def read_haplos(self, file_name, scramble=False):
    
        if not os.path.exists(file_name):
            print "the file: " + file_name + " does not exist!"
            exit(-1)
        
        count = 0
        with open(file_name) as haplos_file:
            for line in haplos_file.xreadlines(  ): 
                count += 1
                
        if count == 0 or count % 2 == 1:
            print "bad number of haplotypes. quitting..."
            exit(-1)
        
        self._nr_haplos = count
        
        print "reading from haplos file: " + file_name
        self._haplos = <char **> malloc(self._nr_haplos * sizeof(char *)) 
        
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
                    #print "line: " + line
                    #print "pos: " + line[0]
                    line_trunc = line[:self._snp_num] 
                    self._haplos[hap_idx] = <char *> malloc((self._snp_num) * sizeof(char))
                    strncpy(self._haplos[hap_idx], line_trunc,self._snp_num) 
                    #py_bytes = line_trunc.encode('UTF-8')
                    #self._haplos[hap_idx] = py_bytes
                    #print self._haplos[hap_idx]
                    hap_idx += 1
                    if hap_idx >= self._nr_haplos:
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
    
    def read_from_bgl_file(self, file_name, anc):
        '''
        read the model from a beagle model file (.dag file)
        '''
        
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
                        #print "reading layer: " + str(layer) + ", snp id: " + nodes_curr[0][1]
                        
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
                            error_eps = self.eps*1e-1  
                            if int(nodes_curr[j][4]) == self._allele_0:
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
                                    self._pi[anc][int(layer / self._win_size)][j] = 1 #sum_back_trans
                        
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
        print "Finished reading beagle file."
#         if self._debug:
#             self.print_transitions()
            
                #p = [self._pi[anc][i] for i in xrange(self._layer_state_nums[anc][0])]
                #print "pi: " + str(p)  
    
    def set_ibs(self, ibs_dic, by_position = True):
        if self._debug: 
            self._ibs_file.write(self._prefix_string + " ")
        for win_idx in range(self.get_num_windows()):
            if by_position:
                inter = ibs_dic.find(self._position[self.start_snp(win_idx)],self._position[self.end_snp(win_idx)-1])
            else:
                inter = ibs_dic.find(self.start_snp(win_idx),self.end_snp(win_idx))
            if len(inter) > 0:
                self._ibs[win_idx] = True
            else:
                self._ibs[win_idx] = False
            #self._ibs[win_idx] = True
            if self._debug:
                self._ibs_file.write(str(int(self._ibs[win_idx])) + " ")
           
        if self._debug:
            self._ibs_file.write("\n")
            self._ibs_file.flush()
    
    def set_ibd_segment_endpoints(self, ibs_dic, by_position = True):
        cdef bool* tmp_ibs = <bool *> malloc(self.get_num_windows() * sizeof(bool))
        for win_idx in range(self.get_num_windows()):
            if by_position:
                inter = ibs_dic.find(self._position[self.start_snp(win_idx)],self._position[self.end_snp(win_idx)-1])
            else:
                inter = ibs_dic.find(self.start_snp(win_idx),self.end_snp(win_idx))
            if len(inter) > 0:
                tmp_ibs[win_idx] = True
            else:
                tmp_ibs[win_idx] = False
        
        self._ibd_segment_start[0] = tmp_ibs[0]
        self._ibd_segment_end[self.get_num_windows()-1] = tmp_ibs[self.get_num_windows()-1]
        for win_idx in range(1,self.get_num_windows()-1):
            # keep only start and end windows as IBS
            if tmp_ibs[win_idx] and not tmp_ibs[win_idx-1]: 
                self._ibd_segment_start[win_idx] = True
            elif tmp_ibs[win_idx] and not tmp_ibs[win_idx+1]:
                self._ibd_segment_end[win_idx] = True
                
        if self._debug:
            for win_idx in range(self.get_num_windows()):
                self._ibs_file.write(str(int(self._ibd_segment_start[win_idx])) + " ")
            self._ibs_file.write("\n")
            for win_idx in range(self.get_num_windows()):
                self._ibs_file.write(str(int(self._ibd_segment_end[win_idx])) + " ")
        
        free(tmp_ibs)
    
    def print_ibs(self):
        ibs_s = "" 
        for win_idx in range(self.get_num_windows()):
            ibs_s += str(int(self._ibs[win_idx])) + " "
        print ibs_s
        
    def print_ibd_segment_start(self):
        ibs_s = "" 
        for win_idx in range(self.get_num_windows()):
            ibs_s += str(int(self._ibd_segment_start[win_idx])) + " "
        print ibs_s
        
    def print_ibd_segment_end(self):
        ibs_s = "" 
        for win_idx in range(self.get_num_windows()):
            ibs_s += str(int(self._ibd_segment_end[win_idx])) + " "
        print ibs_s
        
    def get_haplo_num(self):
        return self._nr_haplos
    
    def get_snp_num(self):
        return self._snp_num
    
    def get_haplo(self,hap_idx):
        #pystring = self._haplos[hap_idx][:self._snp_num].decode('UTF-8')
        pystring = self._haplos[hap_idx][:self._snp_num]
        return pystring
    
    def set_haplo(self,hap_idx,hap):
        strncpy(self._haplos[hap_idx], hap, self._snp_num)
    
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
    
    def print_emissions(self):
        self._ems_file.write("snp ancestry node ems_prob0 ems_prob1\n")
        for snp_idx in range(self._snp_num):
            for admx_idx in range(self.K):
                for node_idx in range(self._layer_state_nums[admx_idx][snp_idx]):
                    self._ems_file.write(str(snp_idx) + " " + 
                                  str(admx_idx) + " " +
                                  str(node_idx) + " " + 
                                  str(self._states[admx_idx][snp_idx][node_idx].prob_em[0]) + " " + 
                                  str(self._states[admx_idx][snp_idx][node_idx].prob_em[1]) + "\n") 
    
    def print_transitions(self):
        self._trans_file.write("snp ancestry node next_node transition\n")
        for snp_idx in range(self._snp_num):
            for admx_idx in range(self.K):
                for node_idx in range(self._layer_state_nums[admx_idx][snp_idx]):
                    print str(snp_idx) + " " + str(admx_idx) + " " + str(node_idx) + " " + str(self._states[admx_idx][snp_idx][node_idx].out_trans_num)
                    for nxt_node_idx in range(self._states[admx_idx][snp_idx][node_idx].out_trans_num):
                        self._trans_file.write(str(snp_idx) + " " + 
                                  str(admx_idx) + " " +
                                  str(node_idx) + " " +
                                  str(nxt_node_idx) + " " + 
                                  str(self._trans[admx_idx][snp_idx][node_idx][nxt_node_idx]) + "\n") 
        self._trans_file.flush()
    
    cpdef generate_random_hap(self, int anc):
        hap = ""
        p = [self._pi[anc][0][i] for i in xrange(self._layer_state_nums[anc][0])]
        #print "layer state nums: " + str(self._layer_state_nums[anc][0])
        #print "pi: " + str(p)  
        cdef int next_node = np.array(p).cumsum().searchsorted(np.random.sample(1))[0]  
        #print "next node: " + str(next_node)
        for layer in range(self._snp_num):
            #print "generating haplotype for layer: " + str(layer)
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
        self._s = <double ****> malloc(self.K * sizeof(double ***))
        for anc in range(self.K):
            self._s[anc] = <double ***> malloc((self.get_num_windows()) * sizeof(double ***))
            for win_idx in range(self.get_num_windows()):
                self._s[anc][win_idx] = <double **> malloc(2 * sizeof(double *))
                for j in range(2):
                    self._s[anc][win_idx][j] = <double *> malloc(2 * sizeof(double))
            
            for win_idx in range(self.get_num_windows()):
                d = self._genetic_dist[self.end_snp(win_idx)] - self._genetic_dist[self.start_snp(win_idx)]
                self._s[anc][win_idx][1][1] = exp(-self._t_1_0[anc] * d)  
                self._s[anc][win_idx][0][0] = exp(-self._t_0_1[anc] * d)
                self._s[anc][win_idx][1][0] = 1 - self._s[anc][win_idx][1][1]
                self._s[anc][win_idx][0][1] = 1 - self._s[anc][win_idx][0][0]
                #for i in range(2):
                #    for j in range(2):
                #        print "ibd trans: " + str(anc) + " " + str(win_idx) + " " + str(i) + " " + str(j) + " : " + str(self._s[anc][win_idx][i][j])
            
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
        
        cdef double *win_recomb
        win_recomb = <double *> malloc((self.get_num_windows()) * sizeof(double))
        for win_idx in range(self.get_num_windows()):
            # the approximated recombination rate accross the window: for smnll distances, r = d/100 (approx.) where d is genetic distance (cM) and r is recombination rate (M/b) 
            win_recomb[win_idx] = (self._genetic_dist[self.end_snp(win_idx)] - self._genetic_dist[self.start_snp(win_idx)])

        self._anc_trans = <double *********> malloc((self.get_num_windows()) * sizeof(double ********))
        for win_idx in range(self.get_num_windows()):
            sum = 0
            #print "calculating ancestry transitions for window: " + str(win_idx)
            self._anc_trans[win_idx] = <double ********> malloc(self.K * sizeof(double *******))
            for admx_idx1 in range(self.K):
                self._anc_trans[win_idx][admx_idx1] = <double *******> malloc(self.K * sizeof(double ******))
                for admx_idx2 in range(self.K):
                    self._anc_trans[win_idx][admx_idx1][admx_idx2] = <double ******> malloc(self.K * sizeof(double *****))
                    for admx_idx3 in range(self.K):
                        self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3] = <double *****> malloc(self.K * sizeof(double ****))
                        for admx_idx4 in range(self.K):
                            self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4] = <double ****> malloc(self.K * sizeof(double ***))
                            for nxt_admx_idx1 in range(self.K):
                                self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1] = <double ***> malloc(self.K * sizeof(double **))
                                for nxt_admx_idx2 in range(self.K):
                                    self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2] = <double **> malloc(self.K * sizeof(double *))
                                    for nxt_admx_idx3 in range(self.K):
                                        self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3] = <double *> malloc(self.K * sizeof(double))
                                        for nxt_admx_idx4 in range(self.K):
                                            if win_recomb[win_idx] > 0:
                                                self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4] = (self.g - 1) * win_recomb[win_idx]
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
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
                            for nxt_admx_idx1 in range(self.K):
                                for nxt_admx_idx2 in range(self.K):
                                    for nxt_admx_idx3 in range(self.K):
                                        for nxt_admx_idx4 in range(self.K):
                                            self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4] = \
                                            self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4] * scale
                                            if not c_isfinite(self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4]):
                                                self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4] = 0.0                        
    
    cpdef emission_prob_ibd_admx_mem_alloc(self, int win_idx, bool post):
        cdef int snp_idx 
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int snp_idx_win
        #cdef double sum
        
        snp_idx_win = 0
        if not post:
            self._emission_prob_ibd_admx = <double *********> malloc(self._win_size * sizeof(double ********))
            start_snp = self.start_snp(win_idx)
            end_snp = self.end_snp(win_idx)
        else:
            self._emission_prob_ibd_admx = <double *********> malloc(3 * self._win_size * sizeof(double ********))
            if self._ibd_segment_start[win_idx]:
                start_snp = self.start_snp(win_idx-2)
                end_snp = self.end_snp(win_idx)
            else:
                start_snp = self.start_snp(win_idx)
                end_snp = self.end_snp(win_idx+2)
        for snp_idx in range(start_snp, end_snp):
            self._emission_prob_ibd_admx[snp_idx_win] = <double ********> malloc(self.K * sizeof(double*******))
            for admx_idx1 in range(self.K):
                self._emission_prob_ibd_admx[snp_idx_win][admx_idx1] = <double *******> malloc(self.K * sizeof(double******))
                for admx_idx2 in range(self.K):
                    self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2] = <double ******> malloc(self.K * sizeof(double*****))
                    for admx_idx3 in range(self.K):
                        self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3] = <double *****> malloc(self.K * sizeof(double****))
                        for admx_idx4 in range(self.K):
                            self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4] = <double ****> malloc((self._layer_state_nums[admx_idx1][snp_idx]) * sizeof(double***))
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1] = <double ***> malloc((self._layer_state_nums[admx_idx2][snp_idx]) * sizeof(double**))
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2] = <double **> malloc((self._layer_state_nums[admx_idx3][snp_idx]) * sizeof(double*))
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3] = <double *> malloc((self._layer_state_nums[admx_idx4][snp_idx]) * sizeof(double))
            snp_idx_win += 1
    
    cpdef emission_prob_ibd_admx_mem_free(self, int win_idx, bool post):
        cdef int snp_idx 
        cdef int node_idx1
        cdef int node_idx2
        cdef int node_idx3
        cdef int node_idx4
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int snp_idx_win
        
        snp_idx_win = 0
        if not post:
            start_snp = self.start_snp(win_idx)
            end_snp = self.end_snp(win_idx)
        else:
            if self._ibd_segment_start[win_idx]:
                start_snp = self.start_snp(win_idx-2)
                end_snp = self.end_snp(win_idx)
            else:
                start_snp = self.start_snp(win_idx)
                end_snp = self.end_snp(win_idx+2)
        for snp_idx in range(start_snp, end_snp):
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        free(self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3])
                                    free(self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2])
                                free(self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1])
                            free(self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4])
                        free(self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3])
                    free(self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2])
                free(self._emission_prob_ibd_admx[snp_idx_win][admx_idx1])
            free(self._emission_prob_ibd_admx[snp_idx_win])
            snp_idx_win += 1
        free(self._emission_prob_ibd_admx)
    
    cpdef calc_emission_probs_ibd_admx(self, int win_idx, bool post):
        
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
        cdef int snp_idx_win
        
        snp_idx_win = 0
        if not post:
            start_snp = self.start_snp(win_idx)
            end_snp = self.end_snp(win_idx)
        else:
            if self._ibd_segment_start[win_idx]:
                start_snp = self.start_snp(win_idx-2)
                end_snp = self.end_snp(win_idx)
            else:
                start_snp = self.start_snp(win_idx)
                end_snp = self.end_snp(win_idx+2)
        for snp_idx in range(start_snp, end_snp):
            #print "calculating emission probs in layer: " + str(snp_idx)
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                            #print "admx_idx1: " + str(admx_idx1) + " admx_idx2: " + str(admx_idx2) + " admx_idx3: " + str(admx_idx3) + " admx_idx4: " + str(admx_idx4) + " snp_idx: " + str(snp_idx) + " node_idx1: " + str(node_idx1) + " node_idx2: " + str(node_idx2) + " node_idx3: " + str(node_idx3) + " node_idx4: " + str(node_idx4)
                                            if not self._phased:
                                                if self._chr1[snp_idx] == self._chr2[snp_idx]:
                                                    self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] = \
                                                    self._states[admx_idx1][snp_idx][node_idx1].prob_em[self._chr1[snp_idx]] * self._states[admx_idx2][snp_idx][node_idx2].prob_em[self._chr2[snp_idx]]
                                                else:
                                                    self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] = \
                                                    self._states[admx_idx1][snp_idx][node_idx1].prob_em[self._chr1[snp_idx]] * self._states[admx_idx2][snp_idx][node_idx2].prob_em[self._chr2[snp_idx]] + \
                                                    self._states[admx_idx1][snp_idx][node_idx1].prob_em[self._chr2[snp_idx]] * self._states[admx_idx2][snp_idx][node_idx2].prob_em[self._chr1[snp_idx]]
                                                
                                                if self._chr3[snp_idx] == self._chr4[snp_idx]:
                                                    self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] *= \
                                                    self._states[admx_idx3][snp_idx][node_idx3].prob_em[self._chr3[snp_idx]] * self._states[admx_idx4][snp_idx][node_idx4].prob_em[self._chr4[snp_idx]]
                                                else:
                                                    self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] *= \
                                                    self._states[admx_idx3][snp_idx][node_idx3].prob_em[self._chr3[snp_idx]] * self._states[admx_idx4][snp_idx][node_idx4].prob_em[self._chr4[snp_idx]] + \
                                                    self._states[admx_idx3][snp_idx][node_idx3].prob_em[self._chr4[snp_idx]] * self._states[admx_idx4][snp_idx][node_idx4].prob_em[self._chr3[snp_idx]]
                                            else:
                                                self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] = \
                                                self._states[admx_idx1][snp_idx][node_idx1].prob_em[self._chr1[snp_idx]] * \
                                                self._states[admx_idx2][snp_idx][node_idx2].prob_em[self._chr2[snp_idx]] * \
                                                self._states[admx_idx3][snp_idx][node_idx3].prob_em[self._chr3[snp_idx]] * \
                                                self._states[admx_idx4][snp_idx][node_idx4].prob_em[self._chr4[snp_idx]]
            snp_idx_win += 1                                                     
    
    cpdef scale_factors_mem_alloc(self, int win_idx, bool post):
        cdef int snp_idx 
        cdef int ibd
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int snp_idx_win
        
        snp_idx_win = 0
        if not post:
            self._scale_factor = <double ******> malloc(self._win_size * sizeof(double *****))
            self._backward_scale_factor = <double ******> malloc(self._win_size * sizeof(double *****))
            start_snp = self.start_snp(win_idx)
            end_snp = self.end_snp(win_idx)
        else:
            self._scale_factor = <double ******> malloc(3 * self._win_size * sizeof(double *****))
            self._backward_scale_factor = <double ******> malloc(3 * self._win_size * sizeof(double *****))
            if self._ibd_segment_start[win_idx]:
                start_snp = self.start_snp(win_idx-2)
                end_snp = self.end_snp(win_idx)
            else:
                start_snp = self.start_snp(win_idx)
                end_snp = self.end_snp(win_idx+2)
        for snp_idx in range(start_snp, end_snp):
            self._scale_factor[snp_idx_win] = <double *****> malloc(self.K * sizeof(double****))
            self._backward_scale_factor[snp_idx_win] = <double *****> malloc(self.K * sizeof(double****))
            for admx_idx1 in range(self.K):
                self._scale_factor[snp_idx_win][admx_idx1] = <double ****> malloc(self.K * sizeof(double***))
                self._backward_scale_factor[snp_idx_win][admx_idx1] = <double ****> malloc(self.K * sizeof(double***))
                for admx_idx2 in range(self.K):
                    self._scale_factor[snp_idx_win][admx_idx1][admx_idx2] = <double ***> malloc(self.K * sizeof(double**))
                    self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2] = <double ***> malloc(self.K * sizeof(double**))
                    for admx_idx3 in range(self.K):
                        self._scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3] = <double **> malloc(self.K * sizeof(double*))
                        self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3] = <double **> malloc(self.K * sizeof(double*))
                        for admx_idx4 in range(self.K):
                            self._scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4]= <double *> malloc((self._layer_state_nums[admx_idx1][snp_idx]) * sizeof(double))
                            self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4]= <double *> malloc((self._layer_state_nums[admx_idx1][snp_idx]) * sizeof(double))
                            for ibd in range(2):
                                self._scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
                                self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
            snp_idx_win += 1
            
    cpdef scale_factors_mem_free(self, int win_idx, bool post):
        cdef int snp_idx 
        cdef int ibd
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int snp_idx_win
        
        snp_idx_win = 0
        if not post:
            start_snp = self.start_snp(win_idx)
            end_snp = self.end_snp(win_idx)
        else:
            if self._ibd_segment_start[win_idx]:
                start_snp = self.start_snp(win_idx-2)
                end_snp = self.end_snp(win_idx)
            else:
                start_snp = self.start_snp(win_idx)
                end_snp = self.end_snp(win_idx+2)
        for snp_idx in range(start_snp, end_snp):
            #print "freeing in snp: " + str(snp_idx)
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
                            free(self._scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4])
                            free(self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4])
                        free(self._scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3])
                        free(self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3])
                    free(self._scale_factor[snp_idx_win][admx_idx1][admx_idx2])
                    free(self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2])
                free(self._scale_factor[snp_idx_win][admx_idx1])
                free(self._backward_scale_factor[snp_idx_win][admx_idx1])
            free(self._scale_factor[snp_idx_win])
            free(self._backward_scale_factor[snp_idx_win])
            snp_idx_win += 1
        free(self._scale_factor)
        free(self._backward_scale_factor)
        
    cdef scale_factors_init(self, int win_idx, bool post):
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
        cdef int snp_idx_win
    
        snp_idx_win = 0
        if not post:
            start_snp = self.start_snp(win_idx)
            end_snp = self.end_snp(win_idx)
        else:
            if self._ibd_segment_start[win_idx]:
                start_snp = self.start_snp(win_idx-2)
                end_snp = self.end_snp(win_idx)
            else:
                start_snp = self.start_snp(win_idx)
                end_snp = self.end_snp(win_idx+2)
        for snp_idx in range(start_snp, end_snp):
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
                            for ibd in range(2):
                                self._scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
                                self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
            snp_idx_win += 1

    
    cpdef forward_probs_mem_alloc(self, int win_idx, bool post):
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
        cdef int snp_idx_win
        
        snp_idx_win = 0
        if not post:
            self._forward_probs_ibd_admx = <double **********> malloc(self._win_size * sizeof(double *********))
            start_snp = self.start_snp(win_idx)
            end_snp = self.end_snp(win_idx)
        else:
            self._forward_probs_ibd_admx = <double **********> malloc(3 * self._win_size * sizeof(double *********))
            if self._ibd_segment_start[win_idx]:
                start_snp = self.start_snp(win_idx-2)
                end_snp = self.end_snp(win_idx)
            else:
                start_snp = self.start_snp(win_idx)
                end_snp = self.end_snp(win_idx+2)
        for snp_idx in range(start_snp, end_snp):
            self._forward_probs_ibd_admx[snp_idx_win] = <double *********> malloc(self.K * sizeof(double********))
            for admx_idx1 in range(self.K):
                self._forward_probs_ibd_admx[snp_idx_win][admx_idx1] = <double ********> malloc(self.K * sizeof(double*******))
                for admx_idx2 in range(self.K):
                    self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2] = <double *******> malloc(self.K * sizeof(double******))
                    for admx_idx3 in range(self.K):
                        self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3] = <double ******> malloc(self.K * sizeof(double*****))
                        for admx_idx4 in range(self.K):
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
            
    cpdef forward_probs_mem_free(self, int win_idx, bool post):
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
        cdef int snp_idx_win
        
        snp_idx_win = 0
        if not post:
            start_snp = self.start_snp(win_idx)
            end_snp = self.end_snp(win_idx)
        else:
            if self._ibd_segment_start[win_idx]:
                start_snp = self.start_snp(win_idx-2)
                end_snp = self.end_snp(win_idx)
            else:
                start_snp = self.start_snp(win_idx)
                end_snp = self.end_snp(win_idx+2)
        for snp_idx in range(start_snp, end_snp):
            #print "freeing in snp: " + str(snp_idx)
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
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
                                                
    cdef forward_probs_init(self, int win_idx, bool post):
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
        cdef int snp_idx_win
    
        snp_idx_win = 0
        if not post:
            start_snp = self.start_snp(win_idx)
            end_snp = self.end_snp(win_idx)
        else:
            if self._ibd_segment_start[win_idx]:
                start_snp = self.start_snp(win_idx-2)
                end_snp = self.end_snp(win_idx)
            else:
                start_snp = self.start_snp(win_idx)
                end_snp = self.end_snp(win_idx+2)
        for snp_idx in range(start_snp, end_snp):
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                            for ibd in range(2):
                                                self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = 0
            snp_idx_win += 1
    
    cpdef calc_forward_probs_ibd_admx(self, int win_idx, bool post):
        
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
        cdef double temp_mult
        cdef double sum
        cdef double eps_or_1_eps
        cdef int snp_idx_win
        cdef int start_snp
        cdef int end_snp

        if not post:
            first_win_idx = win_idx
            start_snp = self.start_snp(win_idx)
            end_snp = self.end_snp(win_idx)
        else:
            first_win_idx = max(0,win_idx-1)
            if self._ibd_segment_start[win_idx]:
                start_snp = self.start_snp(win_idx-2)
                end_snp = self.end_snp(win_idx)
            else:
                start_snp = self.start_snp(win_idx)
                end_snp = self.end_snp(win_idx+2)

        # first layer
        for admx_idx1 in range(self.K):
            for admx_idx2 in range(self.K):
                for admx_idx3 in range(self.K):
                    for admx_idx4 in range(self.K):
                        for node_idx1 in range(self._layer_state_nums[admx_idx1][start_snp]):
                            for node_idx2 in range(self._layer_state_nums[admx_idx2][start_snp]):
                                for node_idx3 in range(self._layer_state_nums[admx_idx3][start_snp]):
                                    for node_idx4 in range(self._layer_state_nums[admx_idx4][start_snp]):
                                        for ibd in range(2):
                                            if ibd == 0:
                                                self._forward_probs_ibd_admx[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = \
                                                self._pi[admx_idx1][first_win_idx][node_idx1] * self._pi[admx_idx2][first_win_idx][node_idx2] * self._pi[admx_idx3][first_win_idx][node_idx3] * self._pi[admx_idx4][first_win_idx][node_idx4] * \
                                                self._emission_prob_ibd_admx[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4]
                                            else:
                                                if admx_idx1 == admx_idx3:
                                                    self._forward_probs_ibd_admx[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = \
                                                    self._pi[admx_idx1][first_win_idx][node_idx1] * self._pi[admx_idx2][first_win_idx][node_idx2] * self._pi[admx_idx4][first_win_idx][node_idx4] * \
                                                    self._emission_prob_ibd_admx[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4]
        # rescaling to avoid underflow
        for admx_idx1 in range(self.K):
            for admx_idx2 in range(self.K):
                for admx_idx3 in range(self.K):
                    for admx_idx4 in range(self.K):
                        for ibd in range(2):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][start_snp]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][start_snp]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][start_snp]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][start_snp]):
                                            self._scale_factor[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] += self._forward_probs_ibd_admx[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd]
                            if self._scale_factor[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] > 0:                 
                                self._scale_factor[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 1.0 / self._scale_factor[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]
                            else: 
                                self._scale_factor[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = DBL_MAX
    
        # all other layers
        snp_idx_win = 0
        for snp_idx in range(start_snp, end_snp-1):
            # calculate forward probabilities
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
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
                                                                        self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] += \
                                                                        self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][prev_node1][prev_node2][prev_node3][prev_node4][prev_ibd] * \
                                                                        self._emission_prob_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] * \
                                                                        self._back_trans[admx_idx1][snp_idx+1][node_idx1][prev_node_idx1] * \
                                                                        self._back_trans[admx_idx2][snp_idx+1][node_idx2][prev_node_idx2] * \
                                                                        self._back_trans[admx_idx3][snp_idx+1][node_idx3][prev_node_idx3] * \
                                                                        self._back_trans[admx_idx4][snp_idx+1][node_idx4][prev_node_idx4]
                                                                    else:
                                                                        #if admx_idx1 == admx_idx3 and get_likely_allele(self._states[admx_idx1][snp_idx+1][node_idx1]) == get_likely_allele(self._states[admx_idx3][snp_idx+1][node_idx3]) and get_likely_allele(self._states[admx_idx1][snp_idx][prev_node_idx1]) == get_likely_allele(self._states[admx_idx3][snp_idx][prev_node_idx3]):
                                                                        if admx_idx1 == admx_idx3 and node_idx1 == node_idx3 and prev_node1 == prev_node3:
                                                                            self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] += \
                                                                            self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][prev_node1][prev_node2][prev_node3][prev_node4][prev_ibd] * \
                                                                            self._emission_prob_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4] * \
                                                                            self._back_trans[admx_idx1][snp_idx+1][node_idx1][prev_node_idx1] * \
                                                                            self._back_trans[admx_idx2][snp_idx+1][node_idx2][prev_node_idx2] * \
                                                                            self._back_trans[admx_idx4][snp_idx+1][node_idx4][prev_node_idx4]
            
            # rescaling to avoid underflow
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
                            for ibd in range(2):
                                for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx+1]):
                                    for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx+1]):
                                        for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx+1]):
                                            for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx+1]):
                                                self._scale_factor[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] += self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd]
                                if self._scale_factor[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] > 0:                 
                                    self._scale_factor[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 1.0 / self._scale_factor[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]
                                else: 
                                    self._scale_factor[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = DBL_MAX
#                                 if self._scale_factor[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] < 1:
#                                     print " snp_idx: " + str(snp_idx) + " admx_idx1: " + str(admx_idx1) + " admx_idx2: " + str(admx_idx2) + " admx_idx3: " + str(admx_idx3) + " admx_idx4: " + str(admx_idx4) + " ibd: " + str(ibd) + "scale factor: " + str(self._scale_factor[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                                    
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx+1]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx+1]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx+1]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx+1]):
                                            for ibd in range(2):
#                                                 if self._scale_factor[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] < 1: #and self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] > 1:
#                                                     print " snp_idx: " + str(snp_idx+1) + " admx_idx1: " + str(admx_idx1) + " admx_idx2: " + str(admx_idx2) + " admx_idx3: " + str(admx_idx3) + " admx_idx4: " + str(admx_idx4) + \
#                                                     " node_idx1: " + str(node_idx1) + " node_idx2: " + str(node_idx2) + " node_idx3: " + str(node_idx3) + " node_idx4: " + str(node_idx4) + " ibd: " + str(ibd) + " forward prob: " + str(self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd])
                                                self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = \
                                                self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] * self._scale_factor[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]
                                                
                                                #if snp_idx == (self.end_snp(win_idx) - 1):
                                                #    self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] += \ 
                                                #    self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd]
            snp_idx_win += 1
    
    cpdef backward_probs_mem_alloc(self, int win_idx, bool post):
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
        cdef int snp_idx_win
        cdef int start_snp
        cdef int end_snp

        snp_idx_win = 0
        if not post:
            self._backward_probs_ibd_admx = <double **********> malloc(self._win_size * sizeof(double *********))
            start_snp = self.start_snp(win_idx)
            end_snp = self.end_snp(win_idx)
        else:
            self._backward_probs_ibd_admx = <double **********> malloc(3 * self._win_size * sizeof(double *********))
            if self._ibd_segment_start[win_idx]:
                start_snp = self.start_snp(win_idx-2)
                end_snp = self.end_snp(win_idx)
            else:
                start_snp = self.start_snp(win_idx)
                end_snp = self.end_snp(win_idx+2)
        for snp_idx in range(start_snp, end_snp):
            self._backward_probs_ibd_admx[snp_idx_win] = <double *********> malloc(self.K * sizeof(double********))
            for admx_idx1 in range(self.K):
                self._backward_probs_ibd_admx[snp_idx_win][admx_idx1] = <double ********> malloc(self.K * sizeof(double*******))
                for admx_idx2 in range(self.K):
                    self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2] = <double *******> malloc(self.K * sizeof(double******))
                    for admx_idx3 in range(self.K):
                        self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3] = <double ******> malloc(self.K * sizeof(double*****))
                        for admx_idx4 in range(self.K):
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
    
    cdef backward_probs_init(self, int win_idx, bool post):
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
        cdef int snp_idx_win
        cdef int start_snp
        cdef int end_snp

        if not post:
            start_snp = self.start_snp(win_idx)
            end_snp = self.end_snp(win_idx)
        else:
            if self._ibd_segment_start[win_idx]:
                start_snp = self.start_snp(win_idx-2)
                end_snp = self.end_snp(win_idx)
            else:
                start_snp = self.start_snp(win_idx)
                end_snp = self.end_snp(win_idx+2)
        
        snp_idx_win = 0
        for snp_idx in range(start_snp, end_snp):
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                            for ibd in range(2):
                                                self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = 0
            snp_idx_win += 1
                             
    cpdef calc_backward_probs_ibd_admx(self, int win_idx, bool post):
        
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
        cdef int snp_idx_win
        cdef int start_snp
        cdef int end_snp
        cdef int win_size

        if not post:
            start_snp = self.start_snp(win_idx)
            end_snp = self.end_snp(win_idx)
        else:
            if self._ibd_segment_start[win_idx]:
                start_snp = self.start_snp(win_idx-2)
                end_snp = self.end_snp(win_idx)
            else:
                start_snp = self.start_snp(win_idx)
                end_snp = self.end_snp(win_idx+2)
        win_size=end_snp-start_snp
        
        # last layer
#         print "backward probs calc in snp num " + str(end_snp - 1)
        for admx_idx1 in range(self.K):
            for admx_idx2 in range(self.K):
                for admx_idx3 in range(self.K):
                    for admx_idx4 in range(self.K):
                        for node_idx1 in range(self._layer_state_nums[admx_idx1][end_snp - 1]):
                            for node_idx2 in range(self._layer_state_nums[admx_idx2][end_snp - 1]):
                                for node_idx3 in range(self._layer_state_nums[admx_idx3][end_snp - 1]):
                                    for node_idx4 in range(self._layer_state_nums[admx_idx4][end_snp - 1]):
                                        for ibd in range(2):
                                            if ibd == 0 or admx_idx1 == admx_idx3:
                                                self._backward_probs_ibd_admx[win_size - 1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = 1
                                                
        # rescaling to avoid underflow
        for admx_idx1 in range(self.K):
            for admx_idx2 in range(self.K):
                for admx_idx3 in range(self.K):
                    for admx_idx4 in range(self.K):
                        for ibd in range(2):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][end_snp - 1]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][end_snp - 1]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][end_snp - 1]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][end_snp - 1]):
                                            self._backward_scale_factor[win_size - 1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] += self._backward_probs_ibd_admx[win_size - 1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd]
                            if self._backward_scale_factor[win_size - 1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] > 0:                 
                                self._backward_scale_factor[win_size - 1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 1.0 / self._backward_scale_factor[win_size - 1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]
                            else: 
                                self._backward_scale_factor[win_size - 1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = DBL_MAX
#                             if ibd ==0 and self._backward_scale_factor[win_size - 1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] < 1: #and self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] > 1:
#                             print " snp_idx: " + str(end_snp - 1) + " admx_idx1: " + str(admx_idx1) + " admx_idx2: " + str(admx_idx2) + " admx_idx3: " + str(admx_idx3) + " admx_idx4: " + str(admx_idx4) + \
#                             " node_idx1: " + str(node_idx1) + " node_idx2: " + str(node_idx2) + " node_idx3: " + str(node_idx3) + " node_idx4: " + str(node_idx4) + " ibd: " + str(ibd) + " backward prob: " + str(self._backward_probs_ibd_admx[win_size - 1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd])                                    
        
        # all other layers
        snp_idx_win = win_size - 2 
        for snp_idx in reversed(range(start_snp, end_snp - 1)):
#             print "backward probs calc in snp num " + str(snp_idx)
            # calculate forward probabilities
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
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
                                                                        self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] += \
                                                                        self._backward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_node1][nxt_node2][nxt_node3][nxt_node4][nxt_ibd] * \
                                                                        self._trans[admx_idx1][snp_idx][node_idx1][nxt_node_idx1] * \
                                                                        self._trans[admx_idx2][snp_idx][node_idx2][nxt_node_idx2] * \
                                                                        self._trans[admx_idx3][snp_idx][node_idx3][nxt_node_idx3] * \
                                                                        self._trans[admx_idx4][snp_idx][node_idx4][nxt_node_idx4] * \
                                                                        self._emission_prob_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_node1][nxt_node2][nxt_node3][nxt_node4] #* \
                                                                        #self._s[snp_idx][ibd][nxt_ibd] 
                                                                    else:
#                                                                         if admx_idx1 == admx_idx3 and node_idx1 == node_idx3 and nxt_node1 == nxt_node3:
                                                                        self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] += \
                                                                        self._backward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_node1][nxt_node2][nxt_node3][nxt_node4][nxt_ibd] * \
                                                                        self._trans[admx_idx1][snp_idx][node_idx1][nxt_node_idx1] * \
                                                                        self._trans[admx_idx2][snp_idx][node_idx2][nxt_node_idx2] * \
                                                                        self._trans[admx_idx4][snp_idx][node_idx4][nxt_node_idx4] * \
                                                                        self._emission_prob_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_node1][nxt_node2][nxt_node3][nxt_node4]
                                                                        #self._s[snp_idx][ibd][nxt_ibd] * \
            # rescaling to avoid underflow
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
                            for ibd in range(2):
                                for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                    for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                        for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                            for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                                self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] += self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd]
                                if self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] > 0:                 
                                    self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 1.0 / self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]
                                else: 
                                    self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = DBL_MAX
#                                 if self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] < 1:
#                                     print " backward scale factor, snp_idx: " + str(snp_idx) + " admx_idx1: " + str(admx_idx1) + " admx_idx2: " + str(admx_idx2) + " admx_idx3: " + str(admx_idx3) + " admx_idx4: " + str(admx_idx4) + " ibd: " + str(ibd) + "scale factor: " + str(self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                                    
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                            for ibd in range(2):
#                                                 if ibd ==0 and self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] < 1: #and self._forward_probs_ibd_admx[snp_idx_win+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] > 1:
#                                                 if self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] > 1:
#                                                     print " snp_idx: " + str(snp_idx) + " admx_idx1: " + str(admx_idx1) + " admx_idx2: " + str(admx_idx2) + " admx_idx3: " + str(admx_idx3) + " admx_idx4: " + str(admx_idx4) + \
#                                                     " node_idx1: " + str(node_idx1) + " node_idx2: " + str(node_idx2) + " node_idx3: " + str(node_idx3) + " node_idx4: " + str(node_idx4) + " ibd: " + str(ibd) + " backward prob: " + str(self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd])
                                                self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] = \
                                                self._backward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd] * self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]                                                
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
        self._top_level_backtrack = <int ******> malloc(self.get_num_windows() * sizeof(int *****))
        self._top_level_scale_factor = <double*> malloc(self.get_num_windows() * sizeof(double))
        for win_idx in range(self.get_num_windows()):
            self._top_level_scale_factor[win_idx] = 0
            if self._ibs[win_idx]:
                self._top_level_ems_prob[win_idx] = <double *****> malloc(self.K * sizeof(double ****))
                self._top_level_forward_probs[win_idx] = <double *****> malloc(self.K * sizeof(double ****))
                self._top_level_backward_probs[win_idx] = <double *****> malloc(self.K * sizeof(double ****))
                self._top_level_backtrack[win_idx] = <int *****> malloc(self.K * sizeof(int ****))
                for admx_idx1 in range(self.K):
                    self._top_level_ems_prob[win_idx][admx_idx1] = <double ****> malloc(self.K * sizeof(double ***))
                    self._top_level_forward_probs[win_idx][admx_idx1] = <double ****> malloc(self.K * sizeof(double ***))
                    self._top_level_backward_probs[win_idx][admx_idx1] = <double ****> malloc(self.K * sizeof(double ***))
                    self._top_level_backtrack[win_idx][admx_idx1] = <int ****> malloc(self.K * sizeof(int ***))
                    for admx_idx2 in range(self.K):
                        self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2] = <double ***> malloc(self.K * sizeof(double **))
                        self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2] = <double ***> malloc(self.K * sizeof(double **))
                        self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2] = <double ***> malloc(self.K * sizeof(double **))
                        self._top_level_backtrack[win_idx][admx_idx1][admx_idx2] = <int ***> malloc(self.K * sizeof(int **))
                        for admx_idx3 in range(self.K):
                            self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3] = <double **> malloc(self.K * sizeof(double *))
                            self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3] = <double **> malloc(self.K * sizeof(double *))
                            self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3] = <double **> malloc(self.K * sizeof(double *))
                            self._top_level_backtrack[win_idx][admx_idx1][admx_idx2][admx_idx3] = <int **> malloc(self.K * sizeof(int *))
                            for admx_idx4 in range(self.K):
                                self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4] = <double *> malloc(2 * sizeof(double))
                                self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4] = <double *> malloc(2 * sizeof(double))
                                self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4] = <double *> malloc(2 * sizeof(double))
                                self._top_level_backtrack[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4] = <int *> malloc(2 * sizeof(int))
                                for ibd in range(2):
                                    self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
                                    self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = -DBL_MAX
                                    self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = -DBL_MAX
                                    self._top_level_backtrack[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
    
    cpdef top_level_init(self):
        cdef int win_idx
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int ibd 
         
        for win_idx in range(self.get_num_windows()):
            #print "initializing window %d" % win_idx
            if self._ibs[win_idx]:
                self._top_level_scale_factor[win_idx] = 0
                for admx_idx1 in range(self.K):
                    for admx_idx2 in range(self.K):
                        for admx_idx3 in range(self.K):
                            for admx_idx4 in range(self.K):
                                for ibd in range(2):
                                    self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
                                    self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = -DBL_MAX
                                    self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = -DBL_MAX   
                                    self._top_level_backtrack[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
    
    cpdef top_level_print(self):
        cdef int win_idx
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int ibd
        
        self._probs_file.write("ind1 ind2 win admx1 admx2 admx3 admx4 ibd forward backward gamma emission alpha1 alpha2 alpha3 alpha4 s1 s2 s3 s4\n")
        for win_idx in range(self.get_num_windows()):
            #print "initializing window %d" % win_idx
            if self._ibs[win_idx]:
                for admx_idx1 in range(self.K):
                    for admx_idx2 in range(self.K):
                        for admx_idx3 in range(self.K):
                            for admx_idx4 in range(self.K):
                                for ibd in range(2):
                                    curr_gamma = self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] + self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] #* (1.0/self._top_level_scale_factor[win_idx])
                                    self._probs_file.write(self._prefix_string + " " + str(win_idx) + " " + 
                                    str(admx_idx1) + " " +  
                                    str(admx_idx2) + " " +  
                                    str(admx_idx3) + " " +  
                                    str(admx_idx4) + " " +  
                                    str(ibd) + " " + 
                                    str(self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]) + " " + 
                                    str(self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]) + " " + 
                                    str(curr_gamma) + " " +
                                    str(self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]) + " " +
                                    str(self._alphas[admx_idx1]) + " " + 
                                    str(self._alphas[admx_idx2]) + " " + 
                                    str(self._alphas[admx_idx3]) + " " + 
                                    str(self._alphas[admx_idx4]) + " " +
                                    str(self._s[0][win_idx][0][0]) + " " + 
                                    str(self._s[0][win_idx][0][1]) + " " + 
                                    str(self._s[0][win_idx][1][0]) + " " + 
                                    str(self._s[0][win_idx][1][1]) + "\n")
        self._probs_file.flush()
    
    cpdef print_inner_probs(self, win_idx, bool post): 
        snp_idx_win = 0
        if not post:
            start_snp = self.start_snp(win_idx)
            end_snp = self.end_snp(win_idx)
        else:
            if self._ibd_segment_start[win_idx]:
                start_snp = self.start_snp(win_idx-2)
                end_snp = self.end_snp(win_idx)
            else:
                start_snp = self.start_snp(win_idx)
                end_snp = self.end_snp(win_idx+2)
        for snp_idx in range(start_snp, end_snp):
            for admx_idx1 in range(self.K):
                for admx_idx2 in range(self.K):
                    for admx_idx3 in range(self.K):
                        for admx_idx4 in range(self.K):
                            for node_idx1 in range(self._layer_state_nums[admx_idx1][snp_idx]):
                                for node_idx2 in range(self._layer_state_nums[admx_idx2][snp_idx]):
                                    for node_idx3 in range(self._layer_state_nums[admx_idx3][snp_idx]):
                                        for node_idx4 in range(self._layer_state_nums[admx_idx4][snp_idx]):
                                            for ibd in range(2):
                                                self._inner_probs_file.write(self._prefix_string + " " + str(snp_idx) + " " + 
                                                          str(admx_idx1) + " " +
                                                          str(admx_idx2) + " " +
                                                          str(admx_idx3) + " " +
                                                          str(admx_idx4) + " " +
                                                          str(node_idx1) + " " +
                                                          str(node_idx2) + " " +
                                                          str(node_idx3) + " " +
                                                          str(node_idx4) + " " +
                                                          str(ibd) + " " + 
                                                          str(self._forward_probs_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4][ibd]) + " " +                                            
                                                          str(self._scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]) + " " + 
                                                          str(self._emission_prob_ibd_admx[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][node_idx1][node_idx2][node_idx3][node_idx4]) + "\n")
            snp_idx_win+=1
        self._inner_probs_file.flush()
    
    cpdef calc_top_level_ems_probs_inner(self):
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
        cdef double tmp
        
        if self._debug:
            self._inner_probs_file.write("ind1 ind2 snp admx1 admx2 admx3 admx4 node1 node2 node3 node4 ibd forward scaling emission\n")
            
        for win_idx in range(self.get_num_windows()):
            if self._ibs[win_idx]:
                self.emission_prob_ibd_admx_mem_alloc(win_idx,False)
                self.calc_emission_probs_ibd_admx(win_idx,False)
                self.scale_factors_mem_alloc(win_idx,False)
                self.forward_probs_mem_alloc(win_idx,False)
                self.calc_forward_probs_ibd_admx(win_idx,False)
                self.backward_probs_mem_alloc(win_idx,False)
                self.calc_backward_probs_ibd_admx(win_idx,False)
                
                if self._debug:
                    self.print_inner_probs(win_idx,False)
                #self.backward_probs_mem_alloc(win_idx)
                #self.calc_backward_probs_ibd_admx(chr1,chr2,chr3,chr4,win_idx)
                for admx_idx1 in range(self.K):
                    for admx_idx2 in range(self.K):
                        for admx_idx3 in range(self.K):
                            for admx_idx4 in range(self.K):
                                for ibd in range(2):
                                    snp_idx_win = 0
                                    tmp = 0
                                    for snp_idx in range(self.start_snp(win_idx), self.end_snp(win_idx)):
                                        self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = \
                                        self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] - \
                                        log(self._scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                                        
                                        tmp = \
                                        tmp - \
                                        log(self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
#                                         print "forward scale factor " + str(snp_idx) + " " + str(admx_idx1) + " " + str(admx_idx2) + " " + str(admx_idx3) + " " + str(admx_idx4) + " " + str(ibd) + " " + str(self._scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
#                                         print "backward scale factor " + str(snp_idx) + " " + str(admx_idx1) + " " + str(admx_idx2) + " " + str(admx_idx3) + " " + str(admx_idx4) + " " + str(ibd) + " " + str(self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                                        snp_idx_win += 1
                                    print str(win_idx) + " " + str(ibd) + " " + str(self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]) + " " + str(tmp)
                                    #self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = exp(self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                self.emission_prob_ibd_admx_mem_free(win_idx,False)
                self.forward_probs_mem_free(win_idx,False)
                self.scale_factors_mem_free(win_idx,False)
           
    cpdef calc_top_level_ems_probs(self, int hap_idx1, int hap_idx2, int hap_idx3, int hap_idx4):
        self.set_chrs(self._haplos[hap_idx1],self._haplos[hap_idx2],self._haplos[hap_idx3],self._haplos[hap_idx4])
        self.calc_top_level_ems_probs_inner()
    
    cpdef calc_post_probs(self):
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
        cdef double before_prob
        cdef double after_prob
        cdef total_prob
        cdef max_prob
        cdef max_breakpoint
        cdef first_win
        cdef last_win
        cdef breakpoint_start
        cdef breakpoint_end
        
        if self._debug:
            self._inner_probs_file.write("ind1 ind2 snp admx1 admx2 admx3 admx4 node1 node2 node3 node4 ibd forward emission\n")
            self._inner_probs_file.flush()
        
        for win_idx in range(self.get_num_windows()):
            if self._ibs[win_idx]:
                self._exact_ibd_starts[win_idx] = self.start_snp(win_idx)
                self._exact_ibd_ends[win_idx] = self.end_snp(win_idx)
        
        for win_idx in range(self.get_num_windows()):
            if self._ibd_segment_start[win_idx] or self._ibd_segment_end[win_idx]:
                
                if self._ibd_segment_start[win_idx]:
                    print "calculating post probs in start window: " + str(win_idx)
                else:
                    print "calculating post probs in end window: " + str(win_idx)
                
                self.emission_prob_ibd_admx_mem_alloc(win_idx, True)
                self.calc_emission_probs_ibd_admx(win_idx, True)
                self.scale_factors_mem_alloc(win_idx, True)
                self.forward_probs_mem_alloc(win_idx, True)
                self.calc_forward_probs_ibd_admx(win_idx, True)
                self.backward_probs_mem_alloc(win_idx, True)
                self.calc_backward_probs_ibd_admx(win_idx, True)
#                 if self._debug:
#                     self.print_inner_probs(win_idx, True)
                #self.backward_probs_mem_alloc(win_idx)
                #self.calc_backward_probs_ibd_admx(chr1,chr2,chr3,chr4,win_idx)
                max_prob = -DBL_MAX
                max_breakpoint = -1
                
                if self._debug:
                    self._breakpoints_file.write(str(win_idx) + " ")
                
                if self._ibd_segment_start[win_idx]:
                    first_win = max(0,win_idx-2)
                    last_win = win_idx
                else:
                    first_win = win_idx
                    last_win = min(win_idx+2,self.get_num_windows()-1)
#                 last_win = min(win_idx+1,self.get_num_windows())
                
                breakpoint_start = self.start_snp(first_win)
                breakpoint_end = self.end_snp(last_win)
                for breakpoint in range(breakpoint_start,breakpoint_end):
#                     print "breakpoint: " + str(breakpoint)
                    for admx_idx1 in range(self.K):
                        for admx_idx2 in range(self.K):
                            for admx_idx3 in range(self.K):
                                for admx_idx4 in range(self.K):
                                    for ibd in range(2):
                                        self._top_level_ems_prob[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = -DBL_MAX
                                        self._top_level_ems_prob[1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = -DBL_MAX
                                        self._top_level_forward_probs[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = -DBL_MAX
                                        self._top_level_forward_probs[1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = -DBL_MAX
                    
#                     if breakpoint > self.start_snp(inner_win_idx):
#                                                 for snp_idx in range(self.start_snp(inner_win_idx), min(self.end_snp(inner_win_idx), breakpoint)):
                    
#                     for inner_win_idx in range(first_win,last_win+1):
                    for admx_idx1 in range(self.K):
                        for admx_idx2 in range(self.K):
                            for admx_idx3 in range(self.K):
                                for admx_idx4 in range(self.K):
                                    ibd = 0 if self._ibd_segment_start[win_idx] else 1
                                    snp_idx_win = 0
                                    self._top_level_ems_prob[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
                                    for snp_idx in range(self.start_snp(first_win), breakpoint):
#                                         if ibd == 1 and self._scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] < 1:
#                                             print "scale 0 : " + str(snp_idx) + " " + str(ibd) + " " + str(log(self._scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]))
                                        self._top_level_ems_prob[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = \
                                        self._top_level_ems_prob[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] - \
                                        log(self._scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                                        snp_idx_win += 1
#                                     print "ems 0 : " + str(self._top_level_ems_prob[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]) 
#                                     print " ems_prob 0: " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_ems_prob[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
#                                     self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = exp(self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
#                                     before_prob = logsumexp(before_prob,self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] + log(self._alphas[admx_idx1]) + log(self._alphas[admx_idx2]) + log(self._alphas[admx_idx3]) + log(self._alphas[admx_idx4]))

                    for admx_idx1 in range(self.K):
                        for admx_idx2 in range(self.K):
                            for admx_idx3 in range(self.K):
                                for admx_idx4 in range(self.K):
                                    ibd = 0 if self._ibd_segment_end[win_idx] else 1
                                    snp_idx_win = 0
                                    self._top_level_ems_prob[1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0
                                    for snp_idx in range(breakpoint, self.end_snp(last_win)):
#                                         if ibd == 1 and self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] < 1:
#                                             print "scale 1 : " + str(snp_idx) + " " + str(ibd) + " " + str(log(self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]))
                                        self._top_level_ems_prob[1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = \
                                        self._top_level_ems_prob[1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] - \
                                        log(self._backward_scale_factor[snp_idx_win][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                                        snp_idx_win += 1
#                                     print "ems 1 : " + str(self._top_level_ems_prob[1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
#                                     print " ems_prob 1: " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_ems_prob[1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                    
                    for admx_idx1 in range(self.K):
                        for admx_idx2 in range(self.K):
                            for admx_idx3 in range(self.K):
                                for admx_idx4 in range(self.K):
                                    for ibd in range(2):
                                        self._top_level_forward_probs[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = \
                                        logsumexp(self._top_level_forward_probs[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd],\
                                                              self._top_level_ems_prob[0][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] + \
                                                              log(self._alphas[admx_idx1]) + log(self._alphas[admx_idx2]) + log(self._alphas[admx_idx3]) + log(self._alphas[admx_idx4]))
                    
                    d = self._genetic_dist[breakpoint] - self._genetic_dist[self.start_snp(first_win)]
                    for anc in range(self.K):
                        self._s[anc][0][1][1] = exp(-self._t_1_0[anc] * d)  
                        self._s[anc][0][0][0] = exp(-self._t_0_1[anc] * d)
                        self._s[anc][0][1][0] = 1 - self._s[anc][0][1][1]
                        self._s[anc][0][0][1] = 1 - self._s[anc][0][0][0]
                    
                    for admx_idx1 in range(self.K):
                        for admx_idx2 in range(self.K):
                            for admx_idx3 in range(self.K):
                                for admx_idx4 in range(self.K):
                                    for ibd in range(2):
                                        for prev_admx_idx1 in range(self.K):
                                            for prev_admx_idx2 in range(self.K):
                                                for prev_admx_idx3 in range(self.K):
                                                    for prev_admx_idx4 in range(self.K):
                                                        for prev_ibd in range(2):
                                                            self._top_level_forward_probs[1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = \
                                                            logsumexp(self._top_level_forward_probs[1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd], \
                                                                      self._top_level_forward_probs[0][prev_admx_idx1][prev_admx_idx2][prev_admx_idx3][prev_admx_idx4][prev_ibd] + \
                                                                      self._top_level_ems_prob[1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] + \
                                                                      log(self._alphas[admx_idx1]) + log(self._alphas[admx_idx2]) + log(self._alphas[admx_idx3]) + log(self._alphas[admx_idx4]))
#                                         print str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_forward_probs[1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                    
                    #calc_top_level_backward_probs(first_win,last_win+1)
                    total_prob = -DBL_MAX
                    for admx_idx1 in range(self.K):
                        for admx_idx2 in range(self.K):
                            for admx_idx3 in range(self.K):
                                for admx_idx4 in range(self.K):
                                    for ibd in range(2):
                                        total_prob = logsumexp(total_prob,self._top_level_forward_probs[1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                    if self._debug:
                        self._breakpoints_file.write(str(total_prob) + "\n")
                    print "breakpoint: " + str(breakpoint) + " no ibd prob:" + str(self._top_level_forward_probs[1][0][0][0][0][0]) +  " ibd prob:" + str(self._top_level_forward_probs[1][0][0][0][0][1]) + " no ibd ems 0:" + str(self._top_level_ems_prob[0][0][0][0][0][0]) +  " ibd ems 0:" + str(self._top_level_ems_prob[0][0][0][0][0][1]) +" no ibd ems 1:" + str(self._top_level_ems_prob[1][0][0][0][0][0]) +  " ibd ems 1:" + str(self._top_level_ems_prob[1][0][0][0][0][1]) + " total prob: " + str(total_prob)
                    if total_prob > max_prob:
                        max_prob = total_prob
                        max_breakpoint = breakpoint
                
                if self._debug:
                    self._breakpoints_file.write("\n")
                    self._breakpoints_file.flush()
                
                if self._ibs[win_idx]:
                    if self._ibd_segment_start[win_idx]:
                        self._exact_ibd_starts[win_idx] = max_breakpoint
#                         if max_breakpoint <= self.end_snp(win_idx):
#                             
#                         else:
#                             self._exact_ibd_starts[win_idx+1] = max_breakpoint
                    else:
                        self._exact_ibd_ends[win_idx] = max_breakpoint
#                         if max_breakpoint >= self.start_snp(win_idx):
#                             
#                         elif max_breakpoint > self._exact_ibd_starts[win_idx-1]:
#                             self._exact_ibd_ends[win_idx-1] = max_breakpoint
                                        
                self.emission_prob_ibd_admx_mem_free(win_idx, True)
                self.forward_probs_mem_free(win_idx, True)
                self.scale_factors_mem_free(win_idx, True)
        
        pairIBD = cPairIBD()
        for win_idx in range(self.get_num_windows()):
            if self._ibs[win_idx]:
                if self._lod_scores[win_idx] > self._min_score:
                    pairIBD.add_interval(self._exact_ibd_starts[win_idx],self._exact_ibd_ends[win_idx],self._lod_scores[win_idx])
                    
        return (pairIBD)
    
    cpdef calc_top_level_viterbi(self):
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
        cdef bool first_ibs_window = False
       
        for win_idx in range(self.get_num_windows()-1):
            pass
#             if self._ibs[win_idx]:
#                 #print "calculating top level forward probs for window: " + str(win_idx)
#                 
#                 if win_idx == 0 or (win_idx > 0 and self._ibs[win_idx] and not self._ibs[win_idx-1]):
#                     for admx_idx1 in range(self.K):
#                         for admx_idx2 in range(self.K):
#                             for admx_idx3 in range(self.K):
#                                 for admx_idx4 in range(self.K):
#                                     for ibd in range(2):
#                                         self._top_level_backtrack[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] += \
#                                         self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] * \
#                                         self._alphas[admx_idx1] * self._alphas[admx_idx2] * self._alphas[admx_idx3] * self._alphas[admx_idx4] * self._ibd_prior[admx_idx1][ibd]
#                 
#                 sum = 0
#                 for admx_idx1 in range(self.K):
#                     for admx_idx2 in range(self.K):
#                         for admx_idx3 in range(self.K):
#                             for admx_idx4 in range(self.K):
#                                 for ibd in range(2):
#                                     for prev_admx_idx1 in range(self.K):
#                                         for prev_admx_idx2 in range(self.K):
#                                             for prev_admx_idx3 in range(self.K):
#                                                 for prev_admx_idx4 in range(self.K):
#                                                     for prev_ibd in range(2):
#                                                         self._top_level_backtrack[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] += \
#                                                         self._top_level_backtrack[win_idx][prev_admx_idx1][prev_admx_idx2][prev_admx_idx3][prev_admx_idx4][prev_ibd] * \
#                                                         self._top_level_ems_prob[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] * \
#                                                         self._s[admx_idx1][win_idx][prev_ibd][ibd] * \
#                                                         self._anc_trans[win_idx][prev_admx_idx1][prev_admx_idx2][prev_admx_idx3][prev_admx_idx4][admx_idx1][admx_idx2][admx_idx3][admx_idx4] * \
#                                                         self._alphas[admx_idx1] * self._alphas[admx_idx2] * self._alphas[admx_idx3] * self._alphas[admx_idx4]
#                                     sum += self._top_level_backtrack[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]
#                                     if False: #ibd == 0:
#                                         print str(win_idx) + " top level ems probs " + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_ems_prob[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
#                                         print str(win_idx) + " top level ibd prior and alphas: " + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._ibd_prior[0]) + " " + str(self._alphas[admx_idx1]) + " " + str(self._alphas[admx_idx2]) + " " + str(self._alphas[admx_idx3]) + " "  + str(self._alphas[admx_idx4])
#                                         print str(win_idx) + " top level ibd trans: " + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._s[admx_idx1][win_idx][0][0]) + " " + str(self._s[admx_idx1][win_idx][0][1]) + " " + str(self._s[admx_idx1][win_idx][1][0]) + " " + str(self._s[admx_idx1][win_idx][1][1])
#                                         for prev_admx_idx1 in range(self.K):
#                                             for prev_admx_idx2 in range(self.K):
#                                                 for prev_admx_idx3 in range(self.K):
#                                                     for prev_admx_idx4 in range(self.K):
#                                                         print str(win_idx) + " top level anc trans: " + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + ": " + str(self._anc_trans[win_idx][prev_admx_idx1][prev_admx_idx2][prev_admx_idx3][prev_admx_idx4][admx_idx1][admx_idx2][admx_idx3][admx_idx4])
#                                         print str(win_idx) + " top level viterbi probs " + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_forward_probs[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
    
    cpdef calc_top_level_forward_probs(self, int start_window, int end_window):
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
        cdef bool first_ibs_window = False
       
        for win_idx in range(start_window, end_window):
            if self._ibs[win_idx]:
                if win_idx == start_window or (win_idx > start_window and self._ibs[win_idx] and not self._ibs[win_idx-1]):
                    for admx_idx1 in range(self.K):
                        for admx_idx2 in range(self.K):
                            for admx_idx3 in range(self.K):
                                for admx_idx4 in range(self.K):
                                    for ibd in range(2):
                                        self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = \
                                        logsumexp(self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd],\
                                                              self._top_level_ems_prob[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] + \
                                                              log(self._alphas[admx_idx1]) + log(self._alphas[admx_idx2]) + log(self._alphas[admx_idx3]) + log(self._alphas[admx_idx4]))
                                                              
                if win_idx < end_window - 1:
                    for admx_idx1 in range(self.K):
                        for admx_idx2 in range(self.K):
                            for admx_idx3 in range(self.K):
                                for admx_idx4 in range(self.K):
                                    for ibd in range(2):
                                        for prev_admx_idx1 in range(self.K):
                                            for prev_admx_idx2 in range(self.K):
                                                for prev_admx_idx3 in range(self.K):
                                                    for prev_admx_idx4 in range(self.K):
                                                        for prev_ibd in range(2):
                                                            self._top_level_forward_probs[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = logsumexp(self._top_level_forward_probs[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd], self._top_level_forward_probs[win_idx][prev_admx_idx1][prev_admx_idx2][prev_admx_idx3][prev_admx_idx4][prev_ibd] + self._top_level_ems_prob[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] + log(self._s[admx_idx1][win_idx][prev_ibd][ibd]) + log(self._anc_trans[win_idx][prev_admx_idx1][prev_admx_idx2][prev_admx_idx3][prev_admx_idx4][admx_idx1][admx_idx2][admx_idx3][admx_idx4]) + log(self._alphas[admx_idx1]) + log(self._alphas[admx_idx2]) + log(self._alphas[admx_idx3]) + log(self._alphas[admx_idx4]))
    
    cpdef calc_top_level_backward_probs(self, int start_window, int end_window):
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
         
        #for win_idx in reversed(range(self.get_num_windows()-1)):
        for win_idx in reversed(range(start_window, end_window-1)):
            if self._ibs[win_idx+1]:
                #print "calculating top level forward probs for window: " + str(win_idx)
                
                if win_idx == end_window-2 or (self._ibs[win_idx+1] and not self._ibs[win_idx+2]):
                    for admx_idx1 in range(self.K):
                        for admx_idx2 in range(self.K):
                            for admx_idx3 in range(self.K):
                                for admx_idx4 in range(self.K):
                                    for ibd in range(2):
#                                         print str(win_idx+1) + " top level backward probs before scaling " + str(win_idx+1) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(1)
                                        self._top_level_backward_probs[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0 #self._top_level_scale_factor[win_idx+1]
#                                         print str(win_idx+1) + " top level backward probs after scaling " + str(win_idx+1) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_backward_probs[win_idx+1][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                
                for admx_idx1 in range(self.K):
                    for admx_idx2 in range(self.K):
                        for admx_idx3 in range(self.K):
                            for admx_idx4 in range(self.K):
                                for ibd in range(2):
                                    for nxt_admx_idx1 in range(self.K):
                                        for nxt_admx_idx2 in range(self.K):
                                            for nxt_admx_idx3 in range(self.K):
                                                for nxt_admx_idx4 in range(self.K):
                                                    for nxt_ibd in range(2):
                                                        self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] =\
                                                        logsumexp(self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd],\
                                                        self._top_level_backward_probs[win_idx+1][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4][nxt_ibd] + \
                                                        self._top_level_ems_prob[win_idx+1][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4][nxt_ibd] + \
                                                        log(self._s[admx_idx1][win_idx][ibd][nxt_ibd]) + \
                                                        log(self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4]) + \
                                                        log(self._alphas[nxt_admx_idx1]) + log(self._alphas[nxt_admx_idx2]) + log(self._alphas[nxt_admx_idx3]) + \
                                                        log(self._alphas[nxt_admx_idx4]))
#                                                         self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] += \
#                                                         self._top_level_backward_probs[win_idx+1][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4][nxt_ibd] * \
#                                                         self._top_level_ems_prob[win_idx+1][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4][nxt_ibd] * \
#                                                         self._s[admx_idx1][win_idx][ibd][nxt_ibd] * \
#                                                         self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4] * \
#                                                         self._alphas[nxt_admx_idx1] * self._alphas[nxt_admx_idx2] * self._alphas[nxt_admx_idx3] * self._alphas[nxt_admx_idx4]
                                    
#                                     print str(win_idx) + " top level backward probs before scaling " + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
#                                     self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = \
#                                     self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] * \
#                                     self._top_level_scale_factor[win_idx+1]
#                                     if not c_isfinite(self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd]):
#                                         self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] = 0.0
#                                     print str(win_idx) + " top level backward probs after scaling " + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                                                        #print ">backward " + str(win_idx) + " " + str(admx_idx1) + " " + str(admx_idx2) + " " + str(admx_idx3) + " " + str(admx_idx4) + " " + str(ibd) + " " + str(nxt_admx_idx1) + " " + str(nxt_admx_idx2) + " " + str(nxt_admx_idx3) + " " + str(nxt_admx_idx4) + " " + str(nxt_ibd) + " ems_prob: " + str(self._top_level_ems_prob[win_idx+1][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4][nxt_ibd])
                                                        #print ">backward " + str(win_idx) + " " + str(admx_idx1) + " " + str(admx_idx2) + " " + str(admx_idx3) + " " + str(admx_idx4) + " " + str(ibd) + " " + str(nxt_admx_idx1) + " " + str(nxt_admx_idx2) + " " + str(nxt_admx_idx3) + " " + str(nxt_admx_idx4) + " " + str(nxt_ibd) + " _s: " + str(self._s[win_idx][ibd][nxt_ibd])
                                                        #print ">backward " + str(win_idx) + " " + str(admx_idx1) + " " + str(admx_idx2) + " " + str(admx_idx3) + " " + str(admx_idx4) + " " + str(ibd) + " " + str(nxt_admx_idx1) + " " + str(nxt_admx_idx2) + " " + str(nxt_admx_idx3) + " " + str(nxt_admx_idx4) + " " + str(nxt_ibd) + " _anc_trans: " + str(self._anc_trans[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4])
                                                        #print ">backward " + str(win_idx) + " " + str(admx_idx1) + " " + str(admx_idx2) + " " + str(admx_idx3) + " " + str(admx_idx4) + " " + str(ibd) + " " + str(nxt_admx_idx1) + " " + str(nxt_admx_idx2) + " " + str(nxt_admx_idx3) + " " + str(nxt_admx_idx4) + " " + str(nxt_ibd) + " _backward prob i+1: " + str(self._top_level_backward_probs[win_idx+1][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4][nxt_ibd])
                                                        #print "yyyyyyyyyyyyyyyyyyyy " + str(win_idx) + " " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + str(self._top_level_backward_probs[win_idx+1][nxt_admx_idx1][nxt_admx_idx2][nxt_admx_idx3][nxt_admx_idx4][nxt_ibd])
                                    #if ibd == 0:
                                    #    print "back top level params: " + str(self._ibd_prior[ibd]) + " " + str(self._alphas[admx_idx1]) + " " + str(self._alphas[admx_idx2]) + " " + str(self._alphas[admx_idx3]) + " "  + str(self._alphas[admx_idx4])
                                    #    print "back top level forward probs " + str(admx_idx1) + " " +  str(admx_idx2) + " " +  str(admx_idx3) + " " +  str(admx_idx4) + " " +  str(ibd) + ": " + str(self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd])
                                                        
    cpdef posterior_top_level_decoding(self):
        cdef int win_idx
        cdef int admx_idx1
        cdef int admx_idx2
        cdef int admx_idx3
        cdef int admx_idx4
        cdef int ibd
        cdef double curr_gamma
        
#         anc_pairs = list(combinations_with_replacement(range(self.K),2))
#         curr_anc1 = [0]*len(anc_pairs)
#         curr_anc2 = [0]*len(anc_pairs)

        ibd_probs = []
        non_ibd_probs = []
        
        pairIBD = cPairIBD()
        for win_idx in range(self.get_num_windows()):
            if self._ibs[win_idx]:
                curr_ibd_prob = -DBL_MAX
                curr_non_ibd_prob = -DBL_MAX
                for admx_idx1 in range(self.K):
                    for admx_idx2 in range(self.K):
                        for admx_idx3 in range(self.K):
                            for admx_idx4 in range(self.K):
                                for ibd in range(2):
#                                     if self._top_level_scale_factor[win_idx] == 0:
#                                         print "top level scare factor is zero, in win_idx: " + str(win_idx) + " ibd: " + str(ibd) 
                                    curr_gamma = self._top_level_forward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] + \
                                    self._top_level_backward_probs[win_idx][admx_idx1][admx_idx2][admx_idx3][admx_idx4][ibd] #* (1.0/self._top_level_scale_factor[win_idx])
                                    if ibd == 0:
                                        #curr_non_ibd_prob += curr_gamma
                                        curr_non_ibd_prob = logsumexp(curr_non_ibd_prob,curr_gamma)
                                    else:
#                                         curr_ibd_prob += curr_gamma
                                        curr_ibd_prob = logsumexp(curr_ibd_prob,curr_gamma)
#                 curr_tot = curr_non_ibd_prob + curr_ibd_prob
#                 if curr_tot > 0:
#                     curr_non_ibd_prob = curr_non_ibd_prob / curr_tot
#                     curr_ibd_prob = curr_ibd_prob / curr_tot
#                                     curr_anc1[anc_pairs.index(min(admx_idx1,admx_idx2),max(admx_idx1,admx_idx2))] += curr_gamma
#                                     curr_anc2[anc_pairs.index(min(admx_idx3,admx_idx4),max(admx_idx3,admx_idx4))] += curr_gamma
#                 #if curr_ibd_prob > curr_non_ibd_prob:
#                 i += '1'
#                 #else:
#                 #    i += '0'
#                     
# #                 (curr_a1,curr_a2) = anc_pairs[curr_anc1.index(max(curr_anc1))]
# #                 (curr_a3,curr_a4) = anc_pairs[curr_anc2.index(max(curr_anc2))]
#                 a1 += '0' #str(curr_a1)
#                 a2 += '0' #str(curr_a2)
#                 a3 += '0' #str(curr_a3)
#                 a4 += '0' #str(curr_a4)
                    
                self._ibd_probs[win_idx] = curr_ibd_prob
                self._no_ibd_probs[win_idx] = curr_non_ibd_prob
                self._lod_scores[win_idx] = 2*(self._ibd_probs[win_idx] - self._no_ibd_probs[win_idx]) 
                if self._lod_scores[win_idx] > self._min_score:
                    pairIBD.add_interval(win_idx*self._win_size,(win_idx+1)*self._win_size,self._lod_scores[win_idx])
            else:
                self._ibd_probs[win_idx] = -1
                self._no_ibd_probs[win_idx] = -1
            
            ibd_probs.append(self._ibd_probs[win_idx])
            non_ibd_probs.append(self._no_ibd_probs[win_idx])
        
#         a1_filt = ''
#         a2_filt = ''
#         a3_filt = ''
#         a4_filt = ''
#         i_filt = ''
#         for ind in range(len(i)):
#             start = max(0,ind-3)
#             end = min(ind+4,len(i))
#             a1_filt += str(int(np.median([int(x) for x in a1[start:end]])))
#             a2_filt += str(int(np.median([int(x) for x in a2[start:end]])))
#             a3_filt += str(int(np.median([int(x) for x in a3[start:end]])))
#             a4_filt += str(int(np.median([int(x) for x in a4[start:end]])))
#             i_filt += str(int(np.median([int(x) for x in i[start:end]])))
                      
        #pairIBD.merge_intervals()
        if self._debug:
            self.top_level_print()
        
        return (pairIBD,ibd_probs,non_ibd_probs)
    
    cdef int start_snp(self, int win_idx):
        return max(0,win_idx * self._win_size)
    
    cdef int end_snp(self, int win_idx):
        return min((win_idx + 1) * self._win_size, self._snp_num)
    
    cpdef int start_position(self):
        return self._position[0] 
    
    cpdef int end_position(self):
        return self._position[self._snp_num-1]
    
    cpdef int get_position(self, int snp_num):
        return self._position[snp_num]
    
    cpdef int get_num_windows(self):
        cdef int num_win
        num_win = int(self._snp_num / self._win_size)
        if self._snp_num % self._win_size > 0:
            num_win += 1
        return num_win
    
    cpdef generate_composite_individuals(self, num_inds):
    
        new_haplos = <char **> malloc(2 * num_inds * sizeof(char *))
        for hap_idx in range(2 * num_inds):
            new_haplos[hap_idx] = <char *> malloc(self._snp_num * sizeof(char))
    
        for i in range(num_inds*2):
            print "composing individual " + str(i)
            print "first position: " + str(self._position[0])
            print "last position: " + str(self._position[self._snp_num-1])
            start = 0
            end = 0
            length = int(200 * (self._position[self._snp_num-1] - self._position[0]) / self._snp_num) 
            same_j = 0
            last_j = -1
            while end < self._snp_num - 1:
                stdout.write("%d," % start)
                stdout.flush()
                start = end 
                end = min(self._snp_num-1,start+length) 
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
        
        self._haplos = new_haplos
        
        
