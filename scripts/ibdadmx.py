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
from numpy import linspace
from IBD.TestSet import GenotypePair  # @UnresolvedImport
from IBD.IBDAdmixedModel import ibdadmixed
from IBD.NaiveModel import naivemodel
from IBD.GeneticMap import GeneticMap
from IBD.IBDSegments import PopIBD
from utils.FormatConversions import convert_ped_to_bgl,convert_ped_to_lamp
from utils.GermlineWrapper import germline

parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(help='sub-command help', dest='command')

parser_a = subparsers.add_parser('ibd', help='detect IBD segments')
parser_a.add_argument("input", type=str, help="prefix of input files (.dat,.map in IBDAdmixed/LAMP format)")
parser_a.add_argument("out", type=str, help="output prefix")
parser_a.add_argument("hapmodelfile", nargs='*', type=str, help=".dag beagle model file name (one for each ancestry)")
parser_a.add_argument("-k", "--num-anc", type=int, dest='K', default=1, help='set number of ancestries')
parser_a.add_argument("-a", "--set-alphas", type=float, nargs='+', dest='alphas', default=[1], help="set the alphas")
parser_a.add_argument("-g", "--generations", type=float, dest='generations', default=8, help="set the number of generations of admixture")
parser_a.add_argument("-n", "--max-snps", type=int, dest='num_snps', default=1000000000, help="maximal number of snps to be used")
parser_a.add_argument("-p", "--num-cpus", type=int, dest='num_cpus', default=1, help="number of cpus to be used")
parser_a.add_argument("-e", "--epsilon", type=float, dest='epsilon', default=1e-4, help='epsilon for error')
parser_a.add_argument("-m", "--min-score", type=float, dest='min_score', default=0, help='minimal score to report as IBD')
parser_a.add_argument("-w", "--win-size", type=int, dest='win_size', default=250, help='window size (in number of SNPs)')
parser_a.add_argument("-o", "--offset", type=int, dest='offset', default=0, help='offset for windows (in number of SNPs)')
parser_a.add_argument("--pair", type=int, nargs=2, dest='pair', help='single pair to process')
parser_a.add_argument("--pairs-file", dest='pairs_file', help='file containing pairs of individuals to process')
parser_a.add_argument("--germline-file", dest='germlinefile', help="germline results file")
parser_a.add_argument("--set-ibd-trans", type=float, nargs='+', dest='ibd_trans', help='set ibd to no IBD probabilities')
parser_a.add_argument("--debug", action='store_true', default=False, dest='debug', help='print debugging information')
parser_a.add_argument("--phased", action='store_true', default=False, dest='phased', help='use phased mode')
parser_a.add_argument("--naive-model", action='store_true', default=False, dest='naive', help='use naive model')
parser_a.add_argument("--scramble", action='store_true', default=False, dest='scramble', help='scramble phase of genotypes')
parser_a.add_argument("--condor", action='store_true', default=False, dest='condor', help='send jobs to htcondor')
parser_a.add_argument("--keep-temp", action='store_true', default=False, dest='keeptemp', help='keep temoporary files')
parser_a.add_argument("--recover", action='store_true', default=False, dest='recover', help='recover unfinished run')

parser_b = subparsers.add_parser('ped2bgl', help='convert ped file to bgl file')
parser_b.add_argument('prefix', type=str, help='plink prefix (name of ped/map files)')

parser_c = subparsers.add_parser('bglmodel', help='run beagle to create the LD model')
parser_c.add_argument('prefix', type=str, help='plink prefix (name of bgl/markers files)')

parser_c = subparsers.add_parser('germline', help='run germline')
parser_c.add_argument('prefix', type=str, help='germline prefix (name of bgl/markers files)')

parser_d = subparsers.add_parser('stats', help='calc stats')
parser_d.add_argument("mapfile", type=str, help="PLINK-format map file name")
parser_d.add_argument("trueibdfile", type=str, help="true ibd file name")
parser_d.add_argument("estimatedibdfile", type=str, help="estimated ibd file name")
parser_d.add_argument("-l", "--min-length", type=float, dest='min_length', default=0, help="minimum length of IBD segments to consider as true IBD (in cM)")
parser_d.add_argument("--filter-est-length", type=float, dest='min_est_length', default=0.8, help="filter short segments from estimated IBD")
parser_d.add_argument("-s", "--min-score", type=float, dest='min_score', default=0, help="minimum score of IBD segments to consider")
parser_d.add_argument("-m", "--max-score", type=float, dest='max_score', default=500, help="maximum score of IBD segments to consider")
parser_d.add_argument("--num-score-points", type=int, dest='num_score_points', default=20, help="Number of score points to calculate stats on")
parser_d.add_argument("--lod-score", action='store_true', default=False, dest='lod_score', help='Score is LOD (the higher the better)')
parser_d.add_argument("--compare-same-inds", action='store_true', default=False, dest='compare_same_inds', help='Calculate stats only on pairs of individuals that appear in both true and estimated IBD results')
parser_d.add_argument("--ibd-admixed", action='store_true', default=False, dest='ibd_admixed', help='')
parser_d.add_argument("--beagle", action='store_true', default=False, dest='beagle', help='')
parser_d.add_argument("--parente", action='store_true', default=False, dest='parente', help='')
parser_d.add_argument("--germline", action='store_true', default=False, dest='germline', help='')

parser_e = subparsers.add_parser('ped2lamp', help='convert plink format to LAMP format')
parser_e.add_argument('prefix', type=str, help='plink prefix (name of ped/map files)')
parser_e.add_argument("--phased", action='store_true', default=False, dest='phased', help='output in phased format')

#parser.add_argument("--log-filename-prefix", dest='log_filename', default="ibdadmixed", help="name of log filename prefix")

if len(sys.argv) > 1:
    args = parser.parse_args()

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
                line = line.replace('\t\t', ' ')
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
        (ibd_p,lod_scores) = ibdadmixed(args.input + ".map",
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
        (ibd_p,lod_scores) = naivemodel(args.input + ".map",
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
                out_lodscores.write(str(ind1) + "\t\t" + \
                                    str(ind2) + "\t\t" + \
                                    str(interval[0]) + "\t\t" + \
                                    str(interval[1] )+ "\t\t" + \
                                    str(win[0]) + "\t\t" + \
                                    str(win[1]) + "\t\t" + \
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

if args.command == "ped2bgl":
    convert_ped_to_bgl(args.prefix + ".ped", args.prefix + ".map", args.prefix + ".bgl", args.prefix + ".markers")

elif args.command == "bglmodel":
    dir = os.path.dirname(__file__)
    filename = os.path.join(dir, '../external/beagle.jar')
    subprocess.call(['java', '-Xmx5000m', '-Djava.io.tmpdir=.', '-jar', filename,
                     'data='+args.prefix+'.bgl',
                     'out='+args.prefix,
                     'scale-2.0',
                     'shift=0'])

if args.command == "germline":
    germline(['germline','-prefix', args.prefix, '-silent', '-bits', '64', '-min_m', '0.1', '-err_hom', '4', '-err_het', '2', '-map', args.prefix+".map", '-w_extend'])

elif args.command == "ibd":

    #args.log_filename = os.path.normpath(args.log_filename)
    #args.log_filename = args.out + ".log"
    logger = logging.getLogger('logger')
    logger.addHandler(logging.FileHandler(args.out + ".log"))
    if args.debug: logger.setLevel(logging.DEBUG)

    if args.input is not None:
        args.input = os.path.normpath(args.input)
    if args.out is not None:
        args.out = os.path.normpath(args.out)
    if args.germlinefile is not None:
        args.germlinefile = os.path.normpath(args.germlinefile)
    if args.pairs_file is not None:
        args.pairs_file = os.path.normpath(args.pairs_file)
    if args.hapmodelfile is not None:
        for anc in range(len(args.hapmodelfile)):
            args.hapmodelfile[anc] = os.path.normpath(args.hapmodelfile[anc])
    if len(args.hapmodelfile) != args.K:
        raise ValueError
    if args.K > 1 and len(args.alphas) != args.K:
        args.alphas = [1.0/args.K] * args.K

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
        gm = GeneticMap(args.input + ".map", args.num_snps)
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
        ts.read_haplos(args.input + ".dat", max_snp_num=args.num_snps)
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

elif args.command == "stats":

    if (args.ibd_admixed + args.beagle + args.parente + args.germline) > 1:
        raise ValueError

    gm = GeneticMap(args.mapfile)

    print "\nFiltering true IBD segments with length < " + str(args.min_length) + "cM"

    true_ibd = PopIBD.fast_deserialize(args.trueibdfile)
    true_ibd.filter_by_length(args.min_length,1e4,gm)

    ibd_est = PopIBD.fast_deserialize(args.estimatedibdfile)

    if args.compare_same_inds:
        print "Considering only individuals that were found in the estimated IBD"
        true_ibd.filter_by_human_pairs(ibd_est.keys())

    scores = [x[2] for x in ibd_est.to_list()]
    print "\nmin score in estimated ibd: " + str(min(scores))
    print "max score in estimated ibd: " + str(max(scores))
    stats = ibd_est.stats_win(true_ibd,gm)
    print "\nStatistics for unfiltered estimated ibd"
    print "power: " + str(stats['power'])
    print "FDR: " + str(stats['FDR'])
    print "FPR: " + str(stats['FPR'])

    ibd_est.filter_by_score(args.min_score,args.max_score)
    ibd_est.merge_all(max_val = args.lod_score)

    scores = linspace(args.min_score,args.max_score,args.num_score_points)
    if args.lod_score:
        scores = reversed(scores)

    print "\nStatistics by min score"
    print "Score\t\tPower\t\tFDR\t\tFPR"
    output_file_name = args.estimatedibdfile + ".stats.txt"
    with open(output_file_name, "w") as output_file:
        for score in scores:
            if args.lod_score:
                ibd_est.filter_by_score(args.min_score,score)
            else:
                ibd_est.filter_by_score(score,args.max_score)
            stats = ibd_est.stats_win(true_ibd,gm)

            line = "{:02.4f}".format(score) + "\t\t" + \
                   "{:02.4f}".format(stats['power']) + "\t\t" + \
                   "{:02.4f}".format(stats['FDR']) + "\t\t" + \
                   "{:02.4f}".format(stats['FPR']) + "\t\t"
            print line
            output_file.write(line+"\n")

elif args.command == "ped2lamp":
    convert_ped_to_lamp(args.prefix + ".ped", args.prefix + ".map", args.prefix + ".lamp", args.phased)