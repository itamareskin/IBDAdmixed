'''
Created on Mar 23, 2013

@author: eskin3
'''

from IBD.LDModel import LDModel
from IBD.cIBD import cPairIBD,cPopulationIBD
#import IBD.LDModelUtils as ldu
#from Logic.IBDGenoHMM import IBDGenoHMM
#import Logic.LDModelUtils as ldu
import math
import numpy as np
from itertools import combinations
from IBD.intersection import Interval, IntervalTree
import string

h = LDModel(map_file_name = "tests/HapMap3_CEU_chr1.map",log_dir = "tests",log_prefix = "tests", k = 1,g = 8,win_size=50,max_snp_num = 500,eps = 1e-4,min_score = -100000000,phased = False,debug = True)
h.set_alphas([1])
h.set_ibd_trans_rate(0,1e-5,1)
#h.set_alphas([0.2,0.8])
#h.read_from_bgl_file("../scripts/HapMap3_CEU_chr1.HapMap3_CEU_chr1.01.bgl.dag",0)
#h.read_from_bgl_file("../scripts/HapMap3_YRI_chr1.HapMap3_YRI_chr1.01.bgl.dag",1)
h.read_from_bgl_file("tests/HapMap.HapMap3_CEU_chr1.bgl.dag",0)
#h.read_from_bgl_file(dir+"hapmap.chr1.yri.hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_yri.unr.phased.all.bgl.dag",1)
#h.read_from_bgl_file("example.data.bgl.dag",1)
#ldu.draw_HMM(h,anc=0,start_level=0,level_num=25)
#ldu.draw_HMM(h,start_level=90,level_num=100)
#chr1 = h.get_haplo(0)
#chr2 = h.get_haplo(1)
#chr3 = h.get_haplo(2)
#chr4 = h.get_haplo(3)

#h.print_transitions(dir+"trans.txt");


h.generate_random_haps_inplace(0,4)
ind1 = 0
ind2 = 1

hap = h.get_haplo(2)[0:200] + h.get_haplo(0)[200:400] + h.get_haplo(2)[400:500]
h.set_haplo(2,hap)

h.calc_ibd_prior()
h.calc_anc_trans()
h.top_level_alloc_mem()

h.top_level_init()
h.set_prefix_string(str(ind1) + " " + str(ind2))
tree = IntervalTree()
tree.add_interval(Interval(h.start_position(),h.end_position()))
h.set_ibs(tree)
#h.calc_top_level_ems_probs_inner(chr1,chr2,chr1,chr3)
#for chr_pair in chr_pairs:
chr_pair = (0,1,0,1)
h.calc_top_level_ems_probs(ind1*2+chr_pair[0],ind1*2+chr_pair[1],ind2*2+chr_pair[2],ind2*2+chr_pair[3])
h.calc_top_level_forward_probs(0,h.get_num_windows())
h.calc_top_level_backward_probs(0,h.get_num_windows())
(ibd,ibd_probs,no_ibd_probs) = h.posterior_top_level_decoding()
#ibd.merge_intervals_fast(merge_diff_vals=True)
#h.set_ibs(ibd,by_position=False,post=True)
tree = IntervalTree()
tree.add_interval(Interval(200,400))
h.set_ibd_segment_endpoints(tree,by_position=False)
ibd_new = h.calc_post_probs()
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
pairs = [(29,81)]
#pairs = combinations(range(h.get_haplo_num()/2),2)
for (ind1,ind2) in pairs:
    #print str(ind1) + "," + str(ind2)
    #chr1 = h.get_haplo(ind1*2)
    #chr2 = h.get_haplo(ind1*2+1)
    #chr3 = h.get_haplo(ind2*2)
    #chr4 = h.get_haplo(ind2*2+1)
    h.top_level_init()
    #h.calc_top_level_ems_probs_inner(chr1,chr2,chr1,chr3,dir+"inner.probs.txt")
    h.calc_top_level_ems_probs(ind1*2,ind1*2+1,ind2*2,ind2*2+1,dir+"inner.probs.txt")
    h.calc_top_level_forward_probs()
    h.calc_top_level_backward_probs()
    (a1,a2,a3,a4,i,ibd,ibd_probs,non_ibd_probs) = h.posterior_top_level_decoding()
    h.top_level_print(dir+"probs.txt")
    if len(ibd.to_list()) > 0:
        popIBD.add_human_pair((ind1,ind2),ibd)

out = open("test.IBDAdmixed.ibd.dat", 'w')
out.writelines(popIBD.to_string())
out.close()


x=1