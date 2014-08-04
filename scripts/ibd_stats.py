#!/a/home/cc/cs/itamares/Software/anaconda/bin/python

'''
@author: eskin
'''

import os
import string
import sys
import argparse
from IBD.IBDSegments import PairIBD, PopIBD
from IBD.GeneticMap import GeneticMap

parser = argparse.ArgumentParser()
parser.add_argument("mapfile", type=str, help="PLINK-format map file name")
parser.add_argument("trueibdfile", type=str, help="true ibd file name")
parser.add_argument("estimatedibdfile", type=str, help="estimated ibd file name")
parser.add_argument("-l", "--min-length", type=float, dest='min_length', default=0, help="minimum length of IBD segments to consider as true IBD (in cM)")
parser.add_argument("--filter-est-length", type=float, dest='min_est_length', default=0.8, help="filter short segments from estimated IBD")
parser.add_argument("-s", "--min-score", type=float, dest='min_score', default=0, help="minimum score of IBD segments to consider")
parser.add_argument("-m", "--max-score", type=float, dest='max_score', default=1e10, help="maximum score of IBD segments to consider")
parser.add_argument("--score-jump", type=float, dest='score_jump', default=0.1, help="score jumps")
parser.add_argument("--lod-score", action='store_true', default=False, dest='lod_score', help='Score is LOD (the higher the better)')
parser.add_argument("--compare-same-inds", action='store_true', default=False, dest='compare_same_inds', help='Calculate stats only on pairs of individuals that appear in both true and estimated IBD results')
parser.add_argument("--ibd-admixed", action='store_true', default=False, dest='ibd_admixed', help='')
parser.add_argument("--beagle", action='store_true', default=False, dest='beagle', help='')
parser.add_argument("--parente", action='store_true', default=False, dest='parente', help='')
parser.add_argument("--germline", action='store_true', default=False, dest='germline', help='')
args = parser.parse_args()

if (args.ibd_admixed + args.beagle + args.parente + args.germline) > 1:
    raise ValueError

gm = GeneticMap(args.mapfile)

true_ibd = PopIBD.fast_deserialize(args.trueibdfile)
true_ibd.filter_by_length(args.min_length,1e4,gm)

ibd_est = PopIBD.fast_deserialize(args.estimatedibdfile)

if args.compare_same_inds:
    true_ibd.filter_by_human_pairs(ibd_est.keys())

ibd_est.filter_by_score(args.min_score,args.max_score)
ibd_est.merge_all(max_val = args.lod_score)
stats = ibd_est.stats_win(true_ibd,gm)
print stats['power'], stats['FDR'], stats['FPR']

scores = range(args.min_score,args.max_score,args.score_jump)
if args.lod_score:
    scores = reversed(scores)

for score in scores:

    if args.lod_score:
        ibd_est.filter_by_score(args.min_score,score)
    else:
        ibd_est.filter_by_score(score,args.max_score)
    stats = ibd_est.stats_win(true_ibd,gm)
    print score, stats['power'], stats['FDR'], stats['FPR']

    (power_sigs,detected_dict) = ibd_admixed3.calc_power(true_ibd)
    (power, FDR, FPR) = ibd_admixed3.stats_win(true_ibd,116430,25)
    #(power,FDR,FPR) = ibd_admixed2.stats(true_ibd,116430)
    print "IBDAdmixed", score, power_sigs, power, FDR, FPR
