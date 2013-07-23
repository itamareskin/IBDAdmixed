'''
Created on Mar 23, 2013

@author: eskin3
'''

from IBD.LDModel import LDModel
from IBD.cIBD import cPairIBD,cPopulationIBD
#from Logic.IBDGenoHMM import IBDGenoHMM
#import Logic.LDModelUtils as ldu
import math
import numpy as np
import os
import string
from itertools import combinations
import sys
from bx.intervals import Interval, IntervalTree
import subprocess
import multiprocessing
import time

def runPair(args):
    
    (ind1,ind2) = args[0]
    print "running pair: " + str(ind1) + "," + str(ind2)
    queue = args[1]
    
    chr_pairs = [(0,1,0,1),(0,1,1,0),(1,0,0,1),(1,0,1,0)]
    
#     if len(sys.argv) >= 4:
#         chr_pairs = [chr_pairs[int(sys.argv[3])]]
        
    h.top_level_init()
    h.set_prefix_string(str(ind1) + " " + str(ind2))
    h.set_ibs(pairs[(ind1,ind2)])
    #h.calc_top_level_ems_probs_inner(chr1,chr2,chr1,chr3)
    #for chr_pair in chr_pairs:
    chr_pair = (0,1,0,1)
    h.calc_top_level_ems_probs(ind1*2+chr_pair[0],ind1*2+chr_pair[1],ind2*2+chr_pair[2],ind2*2+chr_pair[3])
    h.calc_top_level_forward_probs()
    h.calc_top_level_backward_probs()
    (a1,a2,a3,a4,i,ibd,ibd_probs,no_ibd_probs) = h.posterior_top_level_decoding()
    #h.top_level_print()
    
        #popIBD.add_human_pair((ind1,ind2),ibd)
    queue.put(((ind1,ind2),ibd.to_list(),ibd_probs,no_ibd_probs))

num_cpus = 1

if len(sys.argv) != 5:
    print "usage: " + sys.argv[0] + " <num_cpus> <num-snps> <dir> <filename>"
    exit(-1)

num_cpus = int(sys.argv[1])
num_snps = int(sys.argv[2])
dir=sys.argv[3]
file_name = dir + "/" + sys.argv[4]

h = LDModel(num_snps,2,8,25,dir)
h.set_alphas([0.2,0.8])
h.set_ibd_trans_rate(0,1e-5,1)
h.set_ibd_trans_rate(1,2e-4,1)
h.read_from_bgl_file(dir + "/HapMap.HapMap3_CEU_chr1.bgl.01.dag",0)
h.read_from_bgl_file(dir + "/HapMap.HapMap3_YRI_chr1.bgl.01.dag",1)
#h.read_from_bgl_file("hapmap.chr1.ceu.hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_ceu.unr.phased.all.bgl.dag",0)
#h.read_from_bgl_file("hapmap.chr1.yri.hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_yri.unr.phased.all.bgl.dag",1)
h.read_haplos(file_name + ".genos.dat",200)
#h.read_from_bgl_file("example.data.bgl.dag",1)
h.read_genetic_map(dir + "/genetic_map_chr1_b36.txt")
#ldu.draw_HMM(h,start_level=90,level_num=100)

h.calc_ibd_prior()
h.calc_anc_trans()
h.top_level_alloc_mem()

#popIBD = cPopulationIBD()

##if len(sys.argv) >= 3:
#    #pairs = [(int(sys.argv[1]),int(sys.argv[2]))]
#pairs_f = open(file_name + ".pairs.dat")
#pairs = pairs_f.readlines()
#pairs = [x.strip("\n") for x in pairs]
#pairs = [x.split(" ") for x in pairs]
#pairs = [(int(x[0]),int(x[1])) for x in pairs]
#pairs_f.close()
##else:
##    pairs = combinations(range(h.get_haplo_num()/2),2)
#pairs = combinations(range(5),2)

# germline_input = open(file_name + ".germline.run","w")
# germline_input.writelines(["1\n",file_name+".genos.map\n",file_name+".genos.ped\n",file_name+".generated\n"])
# germline_input.close()
# retcode = subprocess.call("germline -silent -bits 50 -min_m 2 -err_hom 0 -err_het 0 < " + file_name + ".germline.run > " + file_name + ".generated.out 2> " + ".generated.err", shell=True)
germline = open(file_name + ".generated.match")
pairs = {}
for counter in range(100):
    line = germline.readline()
    if not line:
        break
    line = line.replace('\t',' ')
    line = line.split(' ')
    pair = (int(line[0]),int(line[2]))
    if not pairs.has_key(pair):
        pairs[pair] = IntervalTree()
    pairs[pair].add_interval(Interval(long(line[5]),long(line[6])))
    
out = open(file_name + ".IBDAdmixed3.dat", 'w')
out_windows = open(file_name + ".IBDAdmixed3.windows.dat", 'w')
out_ibdprobs = open(file_name + ".IBDAdmixed3.ibdprobs.dat", 'w')
out_no_ibdprobs = open(file_name + ".IBDAdmixed3.noibdprobs.dat", 'w')
num_win = int(num_snps / 25)

manager = multiprocessing.Manager()
q = manager.Queue()
if num_cpus == 1:
    result = map(runPair, [(x, q) for x in pairs.keys()])
else:
    pool = multiprocessing.Pool(num_cpus)
    result = pool.map_async(runPair, [(x, q) for x in pairs.keys()])
processed=0
while processed < len(pairs.keys()):
    time.sleep(0.5)
    if not q.empty():
        ((ind1,ind2),ibd,ibd_probs,no_ibd_probs) = q.get()
        if len(ibd) > 0:
            out.write(str(ind1) + "," + str(ind2) + ":" + cPairIBD.from_list(ibd).to_string() + "\n")
            out_windows.write(str(ind1) + " " + str(ind2))
            out_ibdprobs.write(str(ind1) + " " + str(ind2) + " " + string.join([str(x) for x in ibd_probs]," ") + "\n")
            out_ibdprobs.flush()
            out_no_ibdprobs.write(str(ind1) + " " + str(ind2) + " " + string.join([str(x) for x in no_ibd_probs]," ") + "\n")
            out_no_ibdprobs.flush()
            #popIBD.add_human_pair((ind1,ind2),cPairIBD.from_list(ibd))
            
            for win_idx in range(num_win):
                start_snp = win_idx * 25
                end_snp = min((win_idx + 1) * 25, num_snps)
                intersect = cPairIBD.from_list(ibd).find(start_snp,end_snp)
                win_ibd = 0
                if len(intersect) > 0:
                    if (intersect[0].end - intersect[0].start) > 15:
                        win_ibd = 1
                out_windows.write(" " + str(win_ibd))
            out_windows.write("\n") 
            out_windows.flush()
            out.flush()
        processed+=1
    
out.close()
out_windows.close()
out_ibdprobs.close()
out_no_ibdprobs.close()



