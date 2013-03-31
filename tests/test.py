'''
Created on Mar 23, 2013

@author: eskin3
'''

from IBD.LDModel import LDModel
from IBD.cIBD import cPairIBD,cPopulationIBD
import IBD.LDModelUtils as ldu
#from Logic.IBDGenoHMM import IBDGenoHMM
#import Logic.LDModelUtils as ldu
import math
import numpy as np
from itertools import combinations

dir='/home/eskin3/workspace/SimulateAdmixedPopulation/SimulateAdmixedPopulation/Testing/'
h = LDModel(100,2,8,25)
h.set_alphas([0.2,0.8])
h.read_from_bgl_file("../scripts/HapMap3_CEU_chr1.HapMap3_CEU_chr1.01.bgl.dag",0)
h.read_from_bgl_file("../scripts/HapMap3_YRI_chr1.HapMap3_YRI_chr1.01.bgl.dag",1)
#h.read_from_bgl_file(dir+"hapmap.chr1.ceu.hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_ceu.unr.phased.all.bgl.dag",0)
#h.read_from_bgl_file(dir+"hapmap.chr1.yri.hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_yri.unr.phased.all.bgl.dag",1)
h.read_haplos("../scripts/AfricanAmericans4.genos.dat",200)
#h.read_from_bgl_file("example.data.bgl.dag",1)
h.read_genetic_map(dir+"genetic_map_chr1_b36.txt")
#ldu.draw_HMM(h,anc=1,start_level=0,level_num=100)
#ldu.draw_HMM(h,start_level=90,level_num=100)
#chr1 = h.get_haplo(0)
#chr2 = h.get_haplo(1)
#chr3 = h.get_haplo(2)
#chr4 = h.get_haplo(3)


#chr1 = h.generate_random_hap(0)
#chr2 = h.generate_random_hap(0)
#chr3 = h.generate_random_hap(1)
#chr4 = h.generate_random_hap(1)
#chr5 = chr3[0:2000] + chr1[2000:3000] + chr3[3000:5000]
#out = open(dir+"halpos.test2.dat", 'w')
#out.writelines(chr1+"\n") 
#out.writelines(chr2+"\n")
#out.writelines(chr1+"\n")
#out.writelines(chr3+"\n")
#out.writelines(chr3+"\n")
#out.writelines(chr4+"\n")
#out.writelines(chr3+"\n")
#out.writelines(chr5+"\n")
#out.close()
#chr5 = h.generate_random_hap(0)
#chr6 = h.generate_random_hap(0)

h.calc_ibd_prior()
h.calc_anc_trans()
h.top_level_alloc_mem()

popIBD = cPopulationIBD()
pairs = [(16,47)]
#pairs = combinations(range(h.get_haplo_num()/2),2)
for (ind1,ind2) in pairs:
    print str(ind1) + "," + str(ind2)
    chr1 = h.get_haplo(ind1*2+1)
    chr2 = h.get_haplo(ind1*2)
    chr3 = h.get_haplo(ind2*2)
    chr4 = h.get_haplo(ind2*2+1)
    h.top_level_init()
    h.calc_top_level_ems_probs_inner(chr1,chr2,chr1,chr3)
    #h.calc_top_level_ems_probs(ind1*2,ind1*2+1,ind2*2,ind2*2+1)
    h.calc_top_level_forward_probs()
    h.calc_top_level_backward_probs()
    (a1,a2,a3,a4,i,ibd) = h.posterior_top_level_decoding()
    if len(ibd.to_list()) > 0:
        popIBD.add_human_pair((ind1,ind2),ibd)

out = open("test.IBDAdmixed.ibd.dat", 'w')
out.writelines(popIBD.to_string())
out.close()

x=1