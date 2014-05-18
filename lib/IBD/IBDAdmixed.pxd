#cython: profile=True

from libcpp cimport bool

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

cdef class IBDAdmixed(object):
    '''
    The IBD Admixed Model 
    '''
    ## constants

    # admixture fraction
    cdef double *_alphas
    # number of generations since begining of admixture
    cdef int g
    # recombination rate
    cdef double r
    
    # number of populations
    cdef int K
    # transition rate of each ancestry from non-IBD to IBD (per cM)
    cdef double *_t_0_1 
    # transition rate of each ancestry from IBD to non-IBD (per cM)
    cdef double *_t_1_0
    
    cdef int _win_size
    
    cdef char** _haplos
    
    cdef char** _true_ancs
    
    # probability of IBD of each ancestry in the first position (_ibd_prior[anc][0] is the prob. of no IBS in ancestry anc, _ibd_prior[anc][1] is the probability of IBS in ancestry anc)
    cdef double **_ibd_prior
    
    # _s[anc][i][j][k] is the probability of transition from IBD state j (0/1 - no-IBD/IBD) to state k, from snp k-1 to snp k
    cdef double ****_s
    
    cdef double *********_anc_trans
    
    # initial probabilities - size of _pi is _layer_state_nums[0]
    cdef double ***_pi 
    
    # states of the HMM - size of _states is _snp_num (each state holds its emission probabilities and number of transitions to next layer)
    cdef state ***_states
    
    # transition probability matrix - size of _trans is _snps_num - 1 
    # _trans[i][j][k] is the probability of transition from state j in layer (SNP) i to state _trans_idx[i][j][k] in layer i+1 
    cdef double ****_trans
    
    # index of states of edges - size of _trans_idx is _snps_num - 1  
    # _trans_idx[i][j][k] is the index of the state in layer i+1 that the edge represented by _trans[i][j][k] points to
    cdef int ****_trans_idx
    
    # backwards transition probability matrix - size of _back_trans is _snps_num
    # _back_trans[i][j][k] is the probability of transition from state _back_trans_idx[i][j][k] in layer (SNP) i-1 to state j in layer i
    cdef double ****_back_trans
    
    # index of states of back edges - size of _back_trans_idx is _snps_num
    # _back_trans_idx[i][j][k] is the index of the state in layer i-1 that the edge represented by _back_trans[i][j][k] originates from
    cdef int ****_back_trans_idx
    
    # snp IDs - size of _snp_ids is _snp_num
    cdef char **_snp_ids
    
    # number of nodes in each layer (SNP) - size of _layer_state_nums is _snp_num
    cdef int **_layer_state_nums
    
    # physical positions
    cdef int *_position
    
    # genetic distances (cM from 5' end)
    cdef double *_genetic_dist
    
    # number of snps (layers) in the model
    cdef int _snp_num
    
    # number of haplotypes to analyze
    cdef int _nr_haplos
    
    cdef double **********_forward_probs_ibd_admx
    
    cdef double **********_backward_probs_ibd_admx
    
    cdef double *********_emission_prob_ibd_admx
    
    cdef double ******_scale_factor
    cdef double ******_backward_scale_factor
    
    cdef double ******_top_level_ems_prob
    
    cdef double ******_top_level_forward_probs
    cdef double ******_top_level_backward_probs
    
    cdef int ******_top_level_backtrack
    
    cdef double *_top_level_scale_factor
       
    # genetic map
    cdef gen_map_entry *_genetic_map 
    
    # log directory
    cdef char* _log_dir
    
    # log prefix
    cdef char* _log_prefix
    
    # log files
    cdef _inner_probs_file
    cdef _probs_file
    cdef _trans_file
    cdef _ems_file
    cdef _ibs_file
    cdef _breakpoints_file
    
    cdef char* _prefix_string
    
    cdef bool* _ibs
    cdef bool* _ibd_segment_start
    cdef bool* _ibd_segment_end
    cdef int _tot_ibs_windows
    
    cdef int* _exact_ibd_starts
    cdef int* _exact_ibd_ends
    
    cdef double* _ibd_probs
    cdef double* _no_ibd_probs
    cdef double* _lod_scores
    
    cdef bool* _chr1
    cdef bool* _chr2
    cdef bool* _chr3
    cdef bool* _chr4
    
    cdef char _allele_0
    cdef char _allele_1
    
    cdef bool _debug
    
    cdef bool _phased
    
    cdef int _offset
            
    ########
    # methods
    ###########
    
    cpdef set_ibd_trans_rate(self, anc, t_0_1, t_1_0)
    
    cpdef set_alphas(self, alphas)
    
    cpdef set_chrs(self, char* chr1, char* chr2, char* chr3, char* chr4)
    
    cpdef generate_random_hap(self, int anc)
    
    cpdef generate_admixed_random_hap(self)
    
    cpdef calc_ibd_prior(self)
    
    cpdef calc_anc_trans(self)
    
    cpdef emission_prob_ibd_admx_mem_alloc(self, int win_idx, bool post)
    
    cpdef emission_prob_ibd_admx_mem_free(self, int win_idx, bool post)
    
    cpdef calc_emission_probs_ibd_admx(self, int win_idx, bool post)
    
    cpdef forward_probs_mem_alloc(self, int win_idx, bool post)
    
    cpdef forward_probs_mem_free(self, int win_idx, bool post)
    
    cdef forward_probs_init(self, int win_idx, bool post)
    
    cpdef calc_forward_probs_ibd_admx(self, int win_idx, bool post)
    
    cpdef scale_factors_mem_alloc(self, int win_idx, bool post)
    
    cpdef scale_factors_mem_free(self, int win_idx, bool post)
    
    cdef scale_factors_init(self, int win_idx, bool post)
    
    cpdef backward_probs_mem_alloc(self, int win_idx, bool post)
    
    cdef backward_probs_init(self, int win_idx, bool post)
    
    cpdef calc_backward_probs_ibd_admx(self, int win_idx, bool post)
    
    #cpdef posterior_decoding_ibd_admx(self, int win_idx)
    
    #cpdef top_level_ems_prob_alloc_mem(self)
    
    cpdef top_level_alloc_mem(self)
    
    cpdef top_level_init(self)
    
    cpdef top_level_print(self)
    
    cpdef calc_top_level_forward_probs(self, int start_window, int end_window)
    
    cpdef calc_top_level_backward_probs(self, int start_window, int end_window)
    
    cpdef print_inner_probs(self, win_idx, bool post)
    
    cpdef calc_top_level_ems_probs_inner(self)
    
    cpdef calc_top_level_ems_probs(self, int hap_idx1, int hap_idx2, int hap_idx3, int hap_idx4)
    
    cpdef calc_post_probs(self)
    
    cpdef calc_top_level_viterbi(self)
    
    cpdef posterior_top_level_decoding(self)
    
    cdef int start_snp(self, int win_idx)
    
    cdef int end_snp(self, int win_idx)
    
    cpdef int start_position(self)
    
    cpdef int end_position(self)
    
    cpdef int get_position(self, int snp_num)
    
    cpdef int get_num_windows(self)
    
    cpdef generate_composite_individuals(self, num_inds)
        
        