from libcpp cimport bool
from LDModel import LDModel
from LDModel cimport LDModel

cdef class TestSet(object):

    # number of snps (layers) in the model
    cdef public int _snp_num
    
    # number of haplotypes to analyze
    cdef public int _nr_haplos    
    
    cdef char _allele_0
    cdef char _allele_1
    
    cdef bool** _haplos
    
    cpdef TestSet get_slice(self, int start_snp, int snp_num)
    
    cpdef GenotypePair get_genotype_pair(self, int ind1, int ind2)

    cpdef Genotype get_genotype(self, int ind)

    #cpdef generate_random_hap(self, int anc)
    
    cpdef generate_composite_individuals(self, LDModel m, num_inds)
    
cdef class GenotypePair(TestSet):
    
    cdef inline bool chr1(self, snp_idx):
        return self._haplos[0][snp_idx]
    
    cdef inline bool chr2(self, snp_idx):
        return self._haplos[1][snp_idx]
    
    cdef inline bool chr3(self, snp_idx):
        return self._haplos[2][snp_idx]
    
    cdef inline bool chr4(self, snp_idx):
        return self._haplos[3][snp_idx]
    
    cpdef generate_random_haps_inplace(self, LDModel m, int snp_num=*)
    
    cpdef set_ibd_segment(self, int start_snp, int snp_num)
    
cdef class Genotype(TestSet):
    
    cdef inline bool chr1(self, snp_idx):
        return self._haplos[0][snp_idx]
    
    cdef inline bool chr2(self, snp_idx):
        return self._haplos[1][snp_idx]
    
    cpdef generate_random_haps_inplace(self, LDModel m, int snp_num=*)