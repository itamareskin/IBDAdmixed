#!/a/home/cc/cs/itamares/Software/anaconda/bin/python

'''
Created on Mar 23, 2013

@author: eskin3
'''

from IBD.TestSet import TestSet  # @UnresolvedImport
# from Logic.IBDGenoHMM import IBDGenoHMM
# import Logic.LDModelUtils as ldu
import os
import string
import sys
import subprocess
import argparse
from htcondor import job, autorun, write_dag_atexit, set_filename
import shutil
import multiprocessing
import logging
import logging.handlers
from functools import partial
from itertools import islice
from IBD.TestSet import GenotypePair  # @UnresolvedImport
from IBD.IBDAdmixedModel import ibdadmixed
from IBD.NaiveModel import naivemodel
from IBD.GeneticMap import GeneticMap

parser = argparse.ArgumentParser()
parser.add_argument("mapfile", type=str, help="PLINK-format map file name")
parser.add_argument("genofile", type=str, help="genomtypes file name (IBDAdmixed/LAMP format)")
parser.add_argument("out", type=str, help="output prefix")
parser.add_argument("hapmodelfile", nargs='*', type=str, help=".dag beagle model file name (one for each ancestry)")
parser.add_argument("-k", "--num-anc", type=int, dest='K', default=1, help='set number of ancestries')
parser.add_argument("-a", "--set-alphas", type=float, nargs='+', dest='alphas', default=[1], help="set the alphas")
parser.add_argument("-g", "--generations", type=float, dest='generations', default=8, help="set the number of generations of admixture")
parser.add_argument("-n", "--max-snps", type=int, dest='num_snps', default=1000000000, help="maximal number of snps to be used")
parser.add_argument("-p", "--num-cpus", type=int, dest='num_cpus', default=1, help="number of cpus to be used")
parser.add_argument("-e", "--epsilon", type=float, dest='epsilon', default=1e-4, help='epsilon for error')
parser.add_argument("-m", "--min-score", type=float, dest='min_score', default=0, help='minimal score to report as IBD')
parser.add_argument("-w", "--win-size", type=int, dest='win_size', default=250, help='window size (in number of SNPs)')
parser.add_argument("-o", "--offset", type=int, dest='offset', default=0, help='offset for windows (in number of SNPs)')
parser.add_argument("--pair", type=int, nargs=2, dest='pair', help='single pair to process')
parser.add_argument("--pairs-file", dest='pairs_file', help='file containing pairs of individuals to process')
parser.add_argument("--germline-file", dest='germlinefile', help="germline results file")
parser.add_argument("--set-ibd-trans", type=float, nargs='+', dest='ibd_trans', help='set ibd to no IBD probabilities')
parser.add_argument("--debug", action='store_true', default=False, dest='debug', help='print debugging information')
parser.add_argument("--phased", action='store_true', default=False, dest='phased', help='use phased mode')
parser.add_argument("--naive-model", action='store_true', default=False, dest='naive', help='use naive model')
parser.add_argument("--scramble", action='store_true', default=False, dest='scramble', help='scramble phase of genotypes')
parser.add_argument("--condor", action='store_true', default=False, dest='condor', help='send jobs to htcondor')
parser.add_argument("--keep-temp", action='store_true', default=False, dest='keeptemp', help='keep temoporary files')
parser.add_argument("--recover", action='store_true', default=False, dest='recover', help='recover unfinished run')

#parser.add_argument("--log-filename-prefix", dest='log_filename', default="ibdadmixed", help="name of log filename prefix")

if len(sys.argv) > 1:
    args = parser.parse_args()
    args.mapfile = os.path.normpath(args.mapfile)
    args.genofile = os.path.normpath(args.genofile)
    args.out = os.path.normpath(args.out)
    if args.germlinefile is not None:
        args.germlinefile = os.path.normpath(args.germlinefile)
    if args.pairs_file is not None:
        args.pairs_file = os.path.normpath(args.pairs_file)
    for anc in range(len(args.hapmodelfile)):
        args.hapmodelfile[anc] = os.path.normpath(args.hapmodelfile[anc])
    if len(args.hapmodelfile) != args.K:
        raise ValueError
    if args.K > 1 and len(args.alphas) != args.K:
        args.alphas = [1.0/args.K] * args.K
    #args.log_filename = os.path.normpath(args.log_filename)
    #args.log_filename = args.out + ".log"

    logger = logging.getLogger('logger')
    logger.addHandler(logging.FileHandler(args.out + ".log"))
    if args.debug: logger.setLevel(logging.DEBUG)

def get_input_file_name(prefix, pair):
    (outdir, outfilename) = os.path.split(os.path.abspath(args.out))
    return os.path.join(outdir, "tmp_outputs." + outfilename, outfilename + "." + str(pair[0]) + "." + str(pair[1]) + ".genos.dat")

def get_output_file_name(outdir, outfilename, pair, suffix):
    if pair is not None and len(pair) > 0:
        return os.path.join(outdir, "tmp_outputs." + outfilename, outfilename + "." + str(pair[0]) + "." + str(pair[1]) + suffix)
    else:
        return os.path.join(outdir, "tmp_outputs." + outfilename)
    
def intervals_from_germline_file(germlinefile, pairs, pos_dict):
    ibs_intervals = {}
    buffer_size = 5000
    print "reading IBS intervals from GERMLINE file: " + germlinefile
    with open(germlinefile) as germline:
        done = False
        while True:
            if done:
                break
            lines = list(islice(germline, buffer_size))
            if len(lines) == 0:
                done = True
            for line in lines:
                if not line:
                    break
                line = line.replace('\t', ' ')
                line = line.split(' ')
                pair = (int(line[0]), int(line[2]))
                if pair not in pairs:
                    continue
                if not ibs_intervals.has_key(pair):
                    ibs_intervals[pair] = []
                if pos_dict.has_key(long(line[5])) and pos_dict.has_key(long(line[6])):
                    ibs_intervals[pair].append((pos_dict[long(line[5])], pos_dict[long(line[6])]))
    print "finished reading IBS intervals."
    return ibs_intervals

@job(output="output.txt", error="error.txt")
def runPair(pair, input_file_name, args):
    print "Running on pair: " + str(pair)
    logger = logging.getLogger('logger')
    (outdir, outfilename) = os.path.split(os.path.abspath(args.out))
    logger.addHandler(logging.FileHandler(get_output_file_name(outdir, outfilename, pair, ".log")))
    if args.debug: logger.setLevel(logging.DEBUG)
    logger.debug("Running on pair: " + str(pair))
    (ind1, ind2) = pair
    obs_data = GenotypePair()
    obs_data.read_haplos(input_file_name)

    if not args.naive:
        (ibd_p,lod_scores) = ibdadmixed(args.mapfile,
                                     args.hapmodelfile,
                                     obs_data,
                                     ibs_intervals=args.ibs_intervals[pair],
                                     max_snp_num=args.num_snps,
                                     phased=args.phased,
                                     g=args.generations,
                                     alphas=args.alphas,
                                     ibd_trans=args.ibd_trans,
                                     win_size=args.win_size,
                                     min_score=args.min_score)
    else:
        (ibd_p,lod_scores) = naivemodel(args.mapfile,
                                     args.hapmodelfile,
                                     obs_data,
                                     ibs_intervals=args.ibs_intervals[pair],
                                     max_snp_num=args.num_snps,
                                     phased=args.phased,
                                     g=args.generations,
                                     alphas=args.alphas,
                                     win_size=args.win_size,
                                     min_score=args.min_score)

    (outdir, outfilename) = os.path.split(os.path.abspath(args.out))
    with open(get_output_file_name(outdir, outfilename, pair, ".ibdadmixed.txt"), 'w') as out:
        out.write(str(ind1) + "," + str(ind2) + ":" + ibd_p.to_string() + "\n")
    with open(get_output_file_name(outdir, outfilename, pair, ".lodscores.txt"), 'w') as out_lodscores:
        for interval in lod_scores.keys():
            for win in lod_scores[interval].keys():
                out_lodscores.write(str(ind1) + "\t" + \
                                    str(ind2) + "\t" + \
                                    str(interval[0]) + "\t" + \
                                    str(interval[1] )+ "\t" + \
                                    str(win[0]) + "\t" + \
                                    str(win[1]) + "\t" + \
                                    str(lod_scores[interval][win]) + \
                                    "\n")
    
@job(output="output.txt", error="error.txt")
def combine_results(args):
    (outdir, outfilename) = os.path.split(os.path.abspath(args.out))
    
    with open(args.out + ".ibdadmixed.txt", "w") as final_out, open(args.out + ".lodscores.txt", "w") as final_out_lodscores:
        for f in os.listdir(get_output_file_name(outdir, outfilename, [], "")):
            if f.endswith(".ibdadmixed.txt"):
                parts = string.split(f, ".")
                pair = (parts[-4], parts[-3])

                with open(get_output_file_name(outdir, outfilename, pair, ".ibdadmixed.txt")) as out:
                    ibd = out.read()
                    final_out.write(ibd)
                with open(get_output_file_name(outdir, outfilename, pair, ".lodscores.txt")) as out_lodscores:
                    lodscores = out_lodscores.read()
                    final_out_lodscores.write(lodscores)
    
    if not args.keeptemp:
        shutil.rmtree(os.path.join(outdir,"tmp_outputs." + outfilename))

autorun(write_dag=False)
set_filename(os.path.basename(args.out) + ".dag")

if not args.recover:
    with  open(args.out + ".cmdlog", 'w') as logf:
        logf.write(str(args) + "\n")

    temp_path = os.path.join(os.path.dirname(args.out), "tmp_outputs." + os.path.basename(args.out))
    if os.path.exists(temp_path):
        shutil.rmtree(temp_path)
    os.makedirs(temp_path)

    pairs = []
    if args.pair is not None:
        pairs = [tuple(args.pair)]

    if args.pairs_file != None:
        with open(args.pairs_file) as pairs_f:
            pairs = pairs_f.readlines()
            pairs = [x.strip("\n") for x in pairs]
            pairs = [x.split(",") for x in pairs]
            pairs = [(int(x[0]), int(x[1])) for x in pairs]

    ibs_intervals = {}
    gm = GeneticMap(args.mapfile, args.num_snps)
    pos_dict = gm.get_position_dict()
    if args.germlinefile != None:
        ibs_intervals = intervals_from_germline_file(args.germlinefile, pairs, pos_dict)
    args.ibs_intervals = ibs_intervals

    if len(args.ibs_intervals) == 0:
        print "No pair of haplotypes to analyze"
        exit(-1)

    # write inputs for jobs
    print "Writing inputs for individuals jobs..."
    ts = TestSet()
    ts.read_haplos(args.genofile, max_snp_num=args.num_snps)
    for idx, pair in enumerate(args.ibs_intervals.keys()):
        gp = ts.get_genotype_pair(pair[0], pair[1])
        input_file_name = get_input_file_name(args.out,pair)
        gp.write_haplos(input_file_name)
    print "Finished writing inputs for individuals jobs."

if args.condor:
    jobs = []

    if not args.recover:
        for idx, pair in enumerate(args.ibs_intervals.keys()):
            curr_job = runPair.queue(pair, get_input_file_name(args.out,pair), args)
            jobs.append(curr_job)
    
    last_job = combine_results.queue(args)
    for curr_job in jobs:
        last_job.parent(curr_job)
        
    write_dag_atexit()
    
    retcode = subprocess.call("condor_submit_dag -f " + os.path.basename(args.out) + ".dag", shell=True)
else:
    if not args.recover:
        runPairPartial = partial(runPair, args=args)
        if args.num_cpus == 1:
            result = map(runPairPartial,
                         [pair for idx, pair in enumerate(args.ibs_intervals.keys())],
                         [get_input_file_name(args.out,pair) for idx, pair in enumerate(args.ibs_intervals.keys())])
        else:
            pool = multiprocessing.Pool(args.num_cpus)
            results = pool.map_async(runPairPartial,
                                     [pair for idx, pair in enumerate(args.ibs_intervals.keys())],
                                     [get_input_file_name(args.out,pair) for idx, pair in enumerate(args.ibs_intervals.keys())])
            pool.close()
            pool.join()
    combine_results(args)
    
