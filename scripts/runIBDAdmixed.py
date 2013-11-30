#!/home/nasheran/itamares/Software/anaconda/bin/python

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
from IBD.intersection import Interval, IntervalTree
import subprocess
import multiprocessing
import time
import argparse
from htcondor import job, autorun

@job
def runPair(pair):
    
    (ind1,ind2) = pair
    #chr_pairs = [(0,1,0,1),(0,1,1,0),(1,0,0,1),(1,0,1,0)]
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
    (ibd,ibd_probs,no_ibd_probs) = h.posterior_top_level_decoding()
    
    return ((ind1,ind2),ibd.to_list(),ibd_probs,no_ibd_probs)

autorun()

def runPairMultiproc(args):
    (ind1,ind2) = args[0]
    print "running pair " + str(args[2]) + "/" + str(args[3]) + ":" + str(ind1) + "," + str(ind2)
    queue = args[1]
    queue.put(runPair(args[0]))

parser = argparse.ArgumentParser()
parser.add_argument("mapfile", type=str, help="PLINK-format map file name")
parser.add_argument("genofile", type=str, help="genomtypes file name (IBDAdmixed/LAMP format)")
parser.add_argument("out", type=str, help="output prefix")
parser.add_argument("hapmodelfile", nargs='*', type=str, help=".dag beagle model file name (one for each ancestry)")
parser.add_argument("-k", "--num-anc", action='store', dest='K', required=True, help='set number of ancestries')
parser.add_argument("-a", "--set-alphas", nargs='+', dest='alphas', required=True, default=[],help="set the alphas")
parser.add_argument("-p", "--num-cpus", action="store", dest='num_cpus',help="number of cpus to be used")
parser.add_argument("-n", "--max-snps", action="store", dest='num_snps',help="maximal number of snps to be used")
parser.add_argument("-g", "--debug", action='store_true', default=False, dest='debug', help='print debugging information')
parser.add_argument("-e", "--epsilon", action='store', dest='epsilon', help='epsilon for error')
parser.add_argument("-m", "--min-score", action='store', dest='min_score', help='minimal score to report as IBD')
parser.add_argument("-w", "--win-size", action='store', dest='win_size', help='window size (in number of SNPs)')
parser.add_argument("--pairs-file", action='store', dest='pairs_file', help='file containing pairs of individuals to process')
parser.add_argument("--germline-file", action='store', dest='germlinefile', help="germline results file")
parser.add_argument("--set-ibd-trans", nargs='+', dest='ibd_trans', help='set ibd to no IBD probabilities')
parser.add_argument("--phased", action='store_true', default=False, dest='phased', help='use phased mode')
parser.add_argument("--scramble", action='store_true', default=False, dest='scramble', help='scramble phase of genotypes')

args = parser.parse_args()

logf = open(args.out + ".log", 'w')
logf.write(str(args) + "\n")
logf.close()
    
num_cpus = 1
if args.num_cpus != None:
    num_cpus = int(args.num_cpus)
    
num_snps = 1000000000
if args.num_snps != None:
    num_snps = int(args.num_snps)

K = 1
if args.K != None:
    K = int(args.K)
    
epsilon=1e-4
if args.epsilon != None:
    epsilon = float(args.epsilon)
    
min_score=0
if args.min_score != None:
    min_score = float(args.min_score)
    
win_size=25
if args.min_score != None:
    win_size = int(args.win_size)
    
h = LDModel(map_file_name = args.mapfile,log_dir = ".",log_prefix = args.out,k = K,g = 8,win_size=win_size,max_snp_num = num_snps,eps = epsilon,min_score = min_score,phased = args.phased,debug = args.debug)

pair_list = []
if args.pairs_file != None:
    pairs_f = open(args.pairs_file)
    pair_list = pairs_f.readlines()
    pair_list = [x.strip("\n") for x in pair_list]
    pair_list = [x.split(",") for x in pair_list]
    pair_list = [(int(x[0]),int(x[1])) for x in pair_list]
    pairs_f.close()
else:
    pair_list = list(combinations(range(h.get_haplo_num()/2),2))

if args.alphas != None:
    if len(args.alphas) > 0:
        alphas = [float(x) for x in args.alphas]
        h.set_alphas(alphas)

       
if args.ibd_trans != None:
    if len(args.ibd_trans) > 0:
        for k in range(K):
            h.set_ibd_trans_rate(k,float(args.ibd_trans[k*2]),float(args.ibd_trans[k*2+1]))
            
#h.set_ibd_trans_rate(0,1e-5,1)
#h.set_ibd_trans_rate(1,2e-4,1)

for anc in range(K):
    h.read_from_bgl_file(args.hapmodelfile[anc],anc)
#h.print_emissions()
#h.read_from_bgl_file("hapmap.chr1.ceu.hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_ceu.unr.phased.all.bgl.dag",0)
#h.read_from_bgl_file("hapmap.chr1.yri.hapmap3_r2_b36_fwd.consensus.qc.poly.chr1_yri.unr.phased.all.bgl.dag",1)
nr_haplos = h.read_haplos(args.genofile,scramble=args.scramble)
nr_inds = nr_haplos/2
#h.read_from_bgl_file("example.data.bgl.dag",1)
#h.read_genetic_map(dir + "/genetic_map_chr1_b36.txt")

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

# germline_input = open(args.genofile + ".germline.run","w")
# germline_input.writelines(["1\n",aergs.mapfile+"\n",file_name+".genos.ped\n",file_name+".generated\n"])
# germline_input.close()
# retcode = subprocess.call("germline -silent -bits 50 -min_m 2 -err_hom 0 -err_het 0 < " + file_name + ".germline.run > " + file_name + ".generated.out 2> " + ".generated.err", shell=True)

pairs = {}
if args.germlinefile != None:
    with open(args.germlinefile) as germline:
        #for counter in range(nr_inds*(nr_inds-1)/2):
        while True:
            line = germline.readline()
            if not line:
                break
            line = line.replace('\t',' ')
            line = line.split(' ')
            pair = (int(line[0]),int(line[2]))
            if len(pair_list) > 0 and not pair in pair_list:
                continue
            if not pairs.has_key(pair):
                pairs[pair] = IntervalTree()
            #print "adding interval: " + line[5] + "," + line[6]
            pairs[pair].add_interval(Interval(long(line[5]),long(line[6])))
else:
    for pair in pair_list:
        if not pairs.has_key(pair):
            pairs[pair] = IntervalTree()
        pairs[pair].add_interval(Interval(h.start_position(),h.end_position()))
    

out = open(args.out + ".ibdadmixed.txt", 'w')
# out_windows = open(args.out + ".IBDAdmixed3.windows.dat", 'w')
out_ibdprobs = open(args.out + ".ibdprobs.txt", 'w')
out_no_ibdprobs = open(args.out + ".noibdprobs.txt", 'w')
num_win = h.get_num_windows()

manager = multiprocessing.Manager()
q = manager.Queue()
if num_cpus == 1:
    result = map(runPairMultiproc, [(x, q, ind, len(pairs)) for ind,x in enumerate(pairs.keys())])
else:
    pool = multiprocessing.Pool(num_cpus)
    result = pool.map_async(runPairMultiproc, [(x, q, ind, len(pairs)) for ind,x in enumerate(pairs.keys())])
processed=0

print "start processing results from individual pairs"
while processed < len(pairs.keys()):
    time.sleep(0.5)
    if not q.empty():
        ((ind1,ind2),ibd,ibd_probs,no_ibd_probs) = q.get()
        print "processing results from individual pair: " + str(ind1) + " " + str(ind2)
#         out_windows.write(str(ind1) + " " + str(ind2))
        out_ibdprobs.write(str(ind1) + " " + str(ind2) + " " + string.join([str(x) for x in ibd_probs]," ") + "\n")
        out_ibdprobs.flush()
        out_no_ibdprobs.write(str(ind1) + " " + str(ind2) + " " + string.join([str(x) for x in no_ibd_probs]," ") + "\n")
        out_no_ibdprobs.flush()
        #popIBD.add_human_pair((ind1,ind2),cPairIBD.from_list(ibd))
        if len(ibd) > 0:
            out.write(str(ind1) + "," + str(ind2) + ":" + cPairIBD.from_list(ibd).to_string() + "\n")
            out.flush()
        processed+=1
    
out.close()
# out_windows.close()
out_ibdprobs.close()
out_no_ibdprobs.close()

