#cython: profile=True

from libcpp cimport bool
from GeneticMap import GeneticMap
from GeneticMap cimport GeneticMap

cdef extern from "string.h":
    char *strncpy(char *dest, char *src, size_t n)

cdef extern from "structs.h":
    int c_isfinite(double x)

    cdef struct struct_state:
        double prob_em[2]
        int out_trans_num
        int in_trans_num
        bool likely_allele

    ctypedef struct_state state
    
    state create_state(double allele_0_prob, double allele_1_prob, int out_trans_num, int in_trans_num)
    
    cdef struct struct_gen_map_entry:
        int position
        double recomb_rate
        double genetic_dist

    ctypedef struct_gen_map_entry gen_map_entry
    
    gen_map_entry create_gen_map_entry(int position, double recomb_rate, double genetic_dist)

cdef class LDModel(object):
    '''
    Hidden Markov Model for a single ancestral population 
    '''
    ## constants
    
    cdef double eps

    # the proportion of ancestry
    cdef double _alpha    
    # a number representing the ancestry of haplotypes modeled (used by other classes that contain LDModel)
    cdef int _anc
    
    # rate of transition from IBD to No-IBD and vice-versa (notation from the Browning paper)
    cdef double _t_0_1
    cdef double _t_1_0
    
    # initial probabilities - size of _pi is _layer_state_nums[0]
    cdef double *_pi 
    
    # states of the HMM - size of _states is _snp_num (each state holds its emission probabilities and number of transitions to next layer)
    cdef state **_states
    
    # transition probability matrix - size of _trans is _snps_num - 1 
    # _trans[i][j][k] is the probability of transition from state j in layer (SNP) i to state _trans_idx[i][j][k] in layer i+1 
    cdef double ***_trans
    
    # index of states of edges - size of _trans_idx is _snps_num - 1  
    # _trans_idx[i][j][k] is the index of the state in layer i+1 that the edge represented by _trans[i][j][k] points to
    cdef int ***_trans_idx
    
    # backwards transition probability matrix - size of _back_trans is _snps_num
    # _back_trans[i][j][k] is the probability of transition from state _back_trans_idx[i][j][k] in layer (SNP) i-1 to state j in layer i
    cdef double ***_back_trans
    
    # index of states of back edges - size of _back_trans_idx is _snps_num
    # _back_trans_idx[i][j][k] is the index of the state in layer i-1 that the edge represented by _back_trans[i][j][k] originates from
    cdef int ***_back_trans_idx
    
    # snp IDs - size of _snp_ids is _snp_num
    cdef char **_snp_ids
    
    # number of nodes in each layer (SNP) - size of _layer_state_nums is _snp_num
    cdef int *_layer_state_nums
    
    # genetic map
    cdef GeneticMap _gm
    
    # number of snps (layers) in the model
    cdef public int _snp_num
    
    cdef char _allele_0
    cdef char _allele_1
    
    cpdef LDModel get_slice_model(self, int start_snp, int snp_num)
       
    cpdef int start_position(self)
    
    cpdef int end_position(self)
    
    cpdef int get_position(self, int snp_num)
    
    cdef bool* generate_random_hap(self, int max_snp_num)
        
        
