'''
Created on May 25, 2014

@author: Itamar
'''
import os
from IBD.LDModel import LDModel
from IBD.TestSet import GenotypePair
from IBD.IBDAdmixedModel import ibdadmixed

resource_path = "../resources"
ldm = LDModel(os.path.join(resource_path, "HapMap3_CEU_chr2.map"),
                  os.path.join(resource_path, "HapMap3_CEU_chr2.HapMap3_CEU_chr2.bgl.dag.gz"),
                  max_snp_num=1000)

p = GenotypePair()
p.generate_random_haps_inplace(ldm)
#p.read_haplos(os.path.join(resource_path, "ceu.genos.dat"),max_snp_num=1000)
p.set_ibd_segment(0,1000)
(pairIBD,lod_scores) = ibdadmixed(os.path.join(resource_path, "HapMap3_CEU_chr2.map"),
           [os.path.join(resource_path, "HapMap3_CEU_chr2.HapMap3_CEU_chr2.bgl.dag.gz")],
           p,
           g=8,
           alphas=[1],
           ibd_trans=[1e-5,1],
           win_size=250,
           phased=False)

x = 1