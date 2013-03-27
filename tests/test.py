'''
Created on Mar 23, 2013

@author: eskin3
'''

from IBD.LDModel import LDModel
#from Logic.IBDGenoHMM import IBDGenoHMM
#import Logic.LDModelUtils as ldu
import math
import numpy as np

h = LDModel(5000,2,8,25)
h.read_from_bgl_file("/home/eskin3/workspace/SimulateAdmixedPopulation/SimulateAdmixedPopulation/Testing/hapmap.chr1.ceu.hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_ceu.unr.phased.all.bgl.dag",0)
h.read_from_bgl_file("/home/eskin3/workspace/SimulateAdmixedPopulation/SimulateAdmixedPopulation/Testing/hapmap.chr1.yri.hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_yri.unr.phased.all.bgl.dag",1)
h.read_haplos("/home/eskin3/workspace/SimulateAdmixedPopulation/SimulateAdmixedPopulation/Testing/haplos_test.txt",4)
#h.read_from_bgl_file("example.data.bgl.dag",1)
h.read_genetic_map("/home/eskin3/workspace/SimulateAdmixedPopulation/SimulateAdmixedPopulation/Testing/genetic_map_chr1_b36.txt")
#ldu.draw_HMM(h,start_level=90,level_num=100)

#chr1 = h.get_haplo(0)
#chr2 = h.get_haplo(1)
#chr3 = h.get_haplo(2)
#chr4 = h.get_haplo(3)

chr1 = h.generate_random_hap(0)
chr2 = h.generate_random_hap(0)
chr3 = h.generate_random_hap(1)
chr4 = h.generate_random_hap(1)
chr5 = chr3[0:2000] + chr1[2000:3000] + chr3[3000:5000]
#chr5 = h.generate_random_hap(0)
#chr6 = h.generate_random_hap(0)

h.calc_ibd_prior()
h.calc_anc_trans()
h.top_level_alloc_mem()
h.calc_top_level_ems_probs(chr1,chr2,chr5,chr4)
h.calc_top_level_forward_probs()
h.calc_top_level_backward_probs()
(a1,a2,a3,a4,i) = h.posterior_top_level_decoding()
a1_filt = ''
a2_filt = ''
a3_filt = ''
a4_filt = ''
i_filt = ''
for ind in range(len(a1)):
    start = max(0,ind-3)
    end = min(ind+3,len(a1))
    a1_filt += str(int(np.median([int(x) for x in a1[start:end]])))
    a2_filt += str(int(np.median([int(x) for x in a2[start:end]])))
    a3_filt += str(int(np.median([int(x) for x in a3[start:end]])))
    a4_filt += str(int(np.median([int(x) for x in a4[start:end]])))
    i_filt += str(int(np.median([int(x) for x in i[start:end]])))
x=1