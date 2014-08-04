'''
Created on Mar 23, 2013

@author: eskin3
'''

from IBD.LDModel import LDModel
from IBD.GenotypePairModel import GenotypePairModel
from IBD.IBDSegments import PairIBD,PopIBD
from IBD.GeneticMap import GeneticMap
#import IBD.LDModelUtils as ldu
#from Logic.IBDGenoHMM import IBDGenoHMM
#import Logic.LDModelUtils as ldu
import math
import numpy as np
from itertools import combinations
from IBD.intersection import Interval, IntervalTree
import string

gm = GeneticMap("K:\\Projects\\ibdadmixed\\tests\\HapMap3_CEU_chr2.map", max_snp_num=500)

max_snp_num = 10000
start_ibd = 2150
end_ibd = 2550
true_ibd = PairIBD()
true_ibd.add_interval(start_ibd,end_ibd)
h = LDModel(map_file_name = "HapMap3_CEU_TSI_chr2.train.map",log_dir = ".",log_prefix = ".", k = 2,g = 8,win_size=200,max_snp_num = max_snp_num,eps = 1e-4,min_score = -50,phased = False,debug = True)
h.set_alphas([0.5,0.5])
h.set_ibd_trans_rate(0,1e-4,1)
h.set_ibd_trans_rate(1,1e-4,1)
#h.set_alphas([0.2,0.8])
#h.read_from_bgl_file("../scripts/HapMap3_CEU_chr1.HapMap3_CEU_chr1.01.bgl.dag",0)
#h.read_from_bgl_file("../scripts/HapMap3_YRI_chr1.HapMap3_YRI_chr1.01.bgl.dag",1)
h.read_from_bgl_file("tests/HapMap3_CEU_chr2.low.HapMap3_CEU_chr2.bgl.dag",0)
h.read_from_bgl_file("tests/HapMap3_YRI_LWK_chr2.low.HapMap3_YRI_LWK_chr2.train.bgl.dag",1)
#h.read_from_bgl_file(dir+"hapmap.chr1.yri.hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_yri.unr.phased.all.bgl.dag",1)
#h.read_from_bgl_file("example.data.bgl.dag",1)
#ldu.draw_HMM(h,anc=0,start_level=0,level_num=25)
#ldu.draw_HMM(h,start_level=90,level_num=100)
#chr1 = h.get_haplo(0)
#chr2 = h.get_haplo(1)
#chr3 = h.get_haplo(2)
#chr4 = h.get_haplo(3)

#h.print_transitions(dir+"trans.txt");

ind1 = 0
ind2 = 1
# h.read_haplos("tests/test.genos.dat")
# h.read_true_ancs("tests/test.true.ancs.dat")
if True:
    h.generate_admixed_random_haps_inplace(4)
    
    hap = h.get_haplo(2)[0:start_ibd] + h.get_haplo(0)[start_ibd:end_ibd] + h.get_haplo(2)[end_ibd:max_snp_num]
    true_anc = h.get_true_anc(2)[0:start_ibd] + h.get_true_anc(0)[start_ibd:end_ibd] + h.get_true_anc(2)[end_ibd:max_snp_num]
    h.set_haplo(2,hap,true_anc)
    
    out = open("tests/test.trueibd.txt", 'w')
    out.write(str(ind1) + "," + str(ind2) + ":" + str(start_ibd) + "," + str(end_ibd) + "\n")
    out.close()
    
    out = open("tests/test.ped", 'w')
    out_s = ""
    out_s += "0 0 0 0 1 1 " + string.join([x for t in zip([str(int(x)+1) for x in list(h.get_haplo(0))], [str(int(x)+1) for x in list(h.get_haplo(1))]) for x in t], " ") + "\n"
    out_s += "1 1 0 0 1 1 " + string.join([x for t in zip([str(int(x)+1) for x in list(h.get_haplo(2))], [str(int(x)+1) for x in list(h.get_haplo(3))]) for x in t], " ") + "\n"
    out.write(out_s)
    out.close()
    
    out = open("tests/test.genos.dat", 'w')
    out.write(h.get_haplo(0) + "\n")
    out.write(h.get_haplo(1) + "\n")
    out.write(h.get_haplo(2) + "\n")
    out.write(h.get_haplo(3) + "\n")
    out.close()
    
    out = open("tests/test.true.ancs.dat", 'w')
    out.write(h.get_true_anc(0) + "\n")
    out.write(h.get_true_anc(1) + "\n")
    out.write(h.get_true_anc(2) + "\n")
    out.write(h.get_true_anc(3) + "\n")
    out.close()
    
    map_file = open("tests/HapMap3_CEU_chr1.map", 'r')
    out = open("tests/test.map", 'w')
    for i in range(max_snp_num):
        out.write(map_file.readline())
    out.close()


#h.calc_top_level_ems_probs_inner(chr1,chr2,chr1,chr3)
#for chr_pair in chr_pairs:
chr_pair = (0,1,0,1)
all_ibd = PairIBD()
for offset in range(0,200,25):
    h.set_offset(offset)
    
    h.calc_ibd_prior()
    h.calc_anc_trans()
    h.top_level_alloc_mem()
    
    h.top_level_init()
    h.set_prefix_string(str(ind1) + " " + str(ind2))
    tree = IntervalTree()
    tree.add_interval(Interval(h.start_position(),h.end_position()))
    h.set_ibs(tree)
    h.smooth_ibs()
    
    h.calc_top_level_ems_probs(ind1*2+chr_pair[0],ind1*2+chr_pair[1],ind2*2+chr_pair[2],ind2*2+chr_pair[3])
    h.calc_top_level_forward_probs(0,h.get_num_windows())
    h.calc_top_level_backward_probs(0,h.get_num_windows())
    (ibd,ibd_probs,no_ibd_probs) = h.posterior_top_level_decoding()
    print offset, ibd.to_string()
    all_ibd.update(ibd,max_val=True)

out = open("tests/test.ibdadmixed.txt", 'w')
out.write(str(ind1) + "," + str(ind2) + ":" + ibd.to_string() + "\n")
out.close()


#ibd.merge_intervals_fast(merge_diff_vals=True)
#h.set_ibs(ibd,by_position=False,post=True)
tree = IntervalTree()
tree.add_interval(Interval(2300,2700))
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

popIBD = PopIBD()
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