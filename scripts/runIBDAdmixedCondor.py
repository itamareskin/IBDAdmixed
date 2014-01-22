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
import time
import argparse
from htcondor import job, autorun, write_dag_atexit, set_filename
import shutil
import multiprocessing
from functools import partial

parser = argparse.ArgumentParser()
parser.add_argument("mapfile", type=str, help="PLINK-format map file name")
parser.add_argument("genofile", type=str, help="genomtypes file name (IBDAdmixed/LAMP format)")
parser.add_argument("out", type=str, help="output prefix")
parser.add_argument("hapmodelfile", nargs='*', type=str, help=".dag beagle model file name (one for each ancestry)")
parser.add_argument("-k", "--num-anc", action='store', dest='K', required=True, help='set number of ancestries')
parser.add_argument("-a", "--set-alphas", nargs='+', dest='alphas', required=True, default=[],help="set the alphas")
parser.add_argument("-n", "--max-snps", action="store", dest='num_snps',help="maximal number of snps to be used")
parser.add_argument("-p", "--num-cpus", action="store", dest='num_cpus',help="number of cpus to be used")
parser.add_argument("-g", "--debug", action='store_true', default=False, dest='debug', help='print debugging information')
parser.add_argument("-e", "--epsilon", action='store', dest='epsilon', help='epsilon for error')
parser.add_argument("-m", "--min-score", action='store', dest='min_score', help='minimal score to report as IBD')
parser.add_argument("-w", "--win-size", action='store', dest='win_size', help='window size (in number of SNPs)')
parser.add_argument("-o", "--offset", action='store', dest='offset', help='offset for windows (in number of SNPs)')
parser.add_argument("--pairs-file", action='store', dest='pairs_file', help='file containing pairs of individuals to process')
parser.add_argument("--germline-file", action='store', dest='germlinefile', help="germline results file")
parser.add_argument("--set-ibd-trans", nargs='+', dest='ibd_trans', help='set ibd to no IBD probabilities')
parser.add_argument("--phased", action='store_true', default=False, dest='phased', help='use phased mode')
parser.add_argument("--scramble", action='store_true', default=False, dest='scramble', help='scramble phase of genotypes')
parser.add_argument("--condor", action='store_true', default=False, dest='condor', help='send jobs to htcondor')

if len(sys.argv) > 1:
    args = parser.parse_args()

@job(output="output.txt",error="error.txt")
def runPair(pair, map_file_name, log_dir, log_prefix, K, g, win_size, max_snp_num, eps, min_score, phased, debug, alphas, ibd_trans, hapmodelfile, outfile, genofile, scramble, pairs, offset):
    
    (ind1,ind2) = pair
    ibd_all = cPairIBD()
    
    h = LDModel(map_file_name = map_file_name,log_dir = log_dir,log_prefix = log_prefix,k = K,g = g,win_size=win_size,max_snp_num = max_snp_num,eps = eps,min_score = min_score,phased = phased,debug = debug,offset=offset)
    h.set_alphas(alphas)
    for k in range(K):
        h.set_ibd_trans_rate(k,float(ibd_trans[k*2]),float(ibd_trans[k*2+1]))
    for anc in range(K):
        h.read_from_bgl_file(hapmodelfile[anc],anc)
        
    h.read_haplos(genofile,scramble=scramble)
    
    h.calc_ibd_prior()
    h.calc_anc_trans()
    h.top_level_alloc_mem()
    h.top_level_init()
    h.set_prefix_string(str(ind1) + " " + str(ind2))
    h.set_ibs(pairs[(ind1,ind2)])
    h.smooth_ibs()
    
    for offset in range(0,win_size,25):

        #h.calc_top_level_ems_probs_inner(chr1,chr2,chr1,chr3)
        #for chr_pair in chr_pairs:
        chr_pair = (0,1,0,1)
        h.calc_top_level_ems_probs(ind1*2+chr_pair[0],ind1*2+chr_pair[1],ind2*2+chr_pair[2],ind2*2+chr_pair[3])
        h.calc_top_level_forward_probs(0,h.get_num_windows())
        h.calc_top_level_backward_probs(0,h.get_num_windows())
        (ibd,ibd_probs,no_ibd_probs) = h.posterior_top_level_decoding()
        ibd_all.update(ibd,max_val=True)
    
#     ibd.merge_intervals_fast(merge_diff_vals=True)
#     h.set_ibs(ibd,by_position=False,post=True)
#     ibd_new = h.calc_post_probs()
    
    (outdir,outfilename) = os.path.split(os.path.abspath(outfile))
    
    out = open(outdir + "/tmp_outputs."+outfilename+"/"+ outfilename + "." + str(ind1) + "." + str(ind2) + ".ibdadmixed.txt", 'w')
    out.write(str(ind1) + "," + str(ind2) + ":" + ibd_all.to_string() + "\n")
    out.close()
    
    out_ibdprobs = open(outdir + "/tmp_outputs."+outfilename+"/"+ outfilename + "." + str(ind1) + "." + str(ind2) + ".ibdprobs.txt", 'w')
    out_no_ibdprobs = open(outdir + "/tmp_outputs."+outfilename+"/"+ outfilename + "." + str(ind1) + "." + str(ind2) + ".noibdprobs.txt", 'w')
    out_ibdprobs.write(str(ind1) + " " + str(ind2) + " " + string.join([str(x) for x in ibd_probs]," ") + "\n")
    out_no_ibdprobs.write(str(ind1) + " " + str(ind2) + " " + string.join([str(x) for x in no_ibd_probs]," ") + "\n")
    out_ibdprobs.close()
    out_no_ibdprobs.close()
    
@job(output="output.txt",error="error.txt")
def combine_results(outfile):
    (outdir,outfilename) = os.path.split(os.path.abspath(outfile))
    
    final_out = open(outfilename + ".ibdadmixed.txt","w") 
    final_out_ibdprobs = open(outfilename + ".ibdprobs.txt","w")
    final_out_noibdprobs = open(outfilename + ".noibdprobs.txt","w")
    
    for f in os.listdir(outdir + "/tmp_outputs."+outfilename+"/"):
        if f.endswith(".ibdadmixed.txt"):
            parts = string.split(f, ".")
            ind1 = parts[-4]
            ind2 = parts[-3]
            
            out = open(outdir + "/tmp_outputs."+outfilename+"/"+ outfilename + "." + str(ind1) + "." + str(ind2) + ".ibdadmixed.txt")
            out_ibdprobs = open(outdir + "/tmp_outputs."+outfilename+"/"+ outfilename + "." + str(ind1) + "." + str(ind2) + ".ibdprobs.txt")
            out_no_ibdprobs = open(outdir + "/tmp_outputs."+outfilename+"/"+ outfilename + "." + str(ind1) + "." + str(ind2) + ".noibdprobs.txt")
            
            ibd = out.read()
            ibdprobs = out_ibdprobs.read()
            noibdprobs = out_no_ibdprobs.read()
            
            final_out.write(ibd)
            final_out_ibdprobs.write(ibdprobs)
            final_out_noibdprobs.write(noibdprobs)
            
            out.close()
            out_ibdprobs.close()
            out_no_ibdprobs.close()
    
    final_out.close()
    final_out_ibdprobs.close()
    final_out_noibdprobs.close()
    
    shutil.rmtree(outdir + "/tmp_outputs."+outfilename+"/")

#if __name__=='__main__':

autorun(write_dag=False)
set_filename(os.path.basename(args.out)+".dag")

logf = open(args.out + ".log", 'w')
logf.write(str(args) + "\n")
logf.close()

temp_path = os.path.join(os.path.dirname(args.out),"tmp_outputs."+os.path.basename(args.out))
if os.path.exists(temp_path):
    shutil.rmtree(temp_path)
os.makedirs(temp_path)

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
    
offset=0
if args.offset != None:
    offset = int(args.offset)
    
#h = LDModel(map_file_name = args.mapfile,log_dir = ".",log_prefix = args.out,k = K,g = 8,win_size=win_size,max_snp_num = num_snps,eps = epsilon,min_score = min_score,phased = args.phased,debug = args.debug)

pair_list = []
if args.pairs_file != None:
    pairs_f = open(args.pairs_file)
    pair_list = pairs_f.readlines()
    pair_list = [x.strip("\n") for x in pair_list]
    pair_list = [x.split(",") for x in pair_list]
    pair_list = [(int(x[0]),int(x[1])) for x in pair_list]
    pairs_f.close()

if args.alphas != None:
    if len(args.alphas) > 0:
        alphas = [float(x) for x in args.alphas]
            
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


if args.condor:
    jobs = []
    for ind,x in enumerate(pairs.keys()):
        curr_job = runPair.queue(x,args.mapfile,".",args.out,K,8,win_size,num_snps,epsilon,min_score,args.phased,args.debug,alphas,args.ibd_trans,args.hapmodelfile,args.out,args.genofile,args.scramble,pairs,offset)
        jobs.append(curr_job)
    
    last_job = combine_results.queue(args.out)
    for curr_job in jobs:
        last_job.parent(curr_job)
        
    write_dag_atexit()
    
    retcode = subprocess.call("condor_submit_dag -f " + os.path.basename(args.out)+".dag", shell=True)
else:    
    runPairPartial = partial(runPair, map_file_name=args.mapfile, log_dir=".", log_prefix=args.out, K=K, g=8, win_size=win_size, max_snp_num=num_snps, eps=epsilon, min_score=min_score, phased=args.phased, debug=args.debug, alphas=alphas, ibd_trans=args.ibd_trans, hapmodelfile=args.hapmodelfile, outfile=args.out, genofile=args.genofile, scramble=args.scramble, pairs=pairs, offset=offset)
    if num_cpus == 1:
        result = map(runPairPartial, [x for ind,x in enumerate(pairs.keys())])
        print result
    else:
        pool = multiprocessing.Pool(num_cpus)
        results = pool.map_async(runPairPartial, [x for ind,x in enumerate(pairs.keys())])
#             for result in results.get():
#                 print result
        pool.close()
        pool.join()
    combine_results(args.out)
    
