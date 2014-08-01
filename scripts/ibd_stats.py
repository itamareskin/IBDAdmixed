#!/a/home/cc/cs/itamares/Software/anaconda/bin/python

'''
@author: eskin
'''

import os
import string
import sys
import argparse
from IBD.cIBD import cPairIBD, cPopulationIBD
from IBD.GeneticMap import GeneticMap

parser = argparse.ArgumentParser()
parser.add_argument("mapfile", type=str, help="PLINK-format map file name")
parser.add_argument("trueibdfile", type=str, help="true ibd file name")
parser.add_argument("estimatedibdfile", type=str, help="estimated ibd file name")
parser.add_argument("-l", "--min-length", type=float, dest='min_length', default=0, help="minimum length of IBD segments to consider as true IBD (in cM)")

args = parser.parse_args()

gm = GeneticMap(args.mapfile)
true_ibd = cPopulationIBD.fast_deserialize(args.trueibdfile)
true_ibd.filter_by_length(min_length,1e4,gm)

ibd_admixed = cPopulationIBD.fast_deserialize("K:\Data\IBDAdmixed\New4\ceu.tsi.yri.lwk.half2.naive20.ibdadmixed.txt")
true_ibd.filter_by_human_pairs(ibd_admixed.keys())
ibd_admixed.filter_by_human_pairs(true_ibd.keys())
#ibd_admixed.merge_all(overlap = 1, max_val = True)
# true_ibd.filter_by_human_pairs(ibd_admixed.keys())

ibd_admixed2 = cPopulationIBD.from_string(ibd_admixed.to_string())
scores = [x[2] for x in ibd_admixed.to_list()]
print min(scores),max(scores)
ibd_admixed2.filter_by_score(-25,1e10)
ibd_admixed2.merge_all_fast(overlap = 1, max_val = True, merge_diff_vals=True)
ibd_admixed2.filter_by_score(-25,1e10)
ibd_admixed2.filter_by_length(0.5,1000,gm)
(power_sigs,detected_dict) = ibd_admixed2.calc_power(true_ibd)
(power, FDR, FPR) = ibd_admixed2.stats_win(true_ibd,116430,25)
#(power,FDR,FPR) = ibd_admixed2.stats(true_ibd,116430)
print power_sigs, power, FDR, FPR

# f = open("/home/eskin/Data/IBDAdmixed/fromLecs/New/ceu.yri.tsilwkref.ibdadmixed.filt.txt","w")
# f.write(ibd_admixed2.to_string())
# f.close()

scores = [int(x[2]) for x in ibd_admixed2.to_list() if (not np.isnan(x[2]) and not np.isinf(x[2]))]
print min(scores),max(scores)
#for score in range(0,205,5):
#for score in range(-50,max(scores)+1,10):
for score in range(min(scores),max(scores)+1,2):
    ibd_admixed2 = cPopulationIBD.from_string(ibd_admixed.to_string(),dists)
    ibd_admixed2.filter_by_score(-25,max(scores)+100)
    ibd_admixed2.merge_all_fast(overlap = 1, max_val = True, merge_diff_vals=True)
    ibd_admixed2.filter_by_length(0.8,1000,gm)

    ibd_admixed3 = cPopulationIBD.from_string(ibd_admixed.to_string(),dists)
    ibd_admixed3.filter_by_score(score,max(scores)+100)
    ibd_admixed3.filter_by_other_ibd(ibd_admixed2)

    (power_sigs,detected_dict) = ibd_admixed3.calc_power(true_ibd)
    (power, FDR, FPR) = ibd_admixed3.stats_win(true_ibd,116430,25)
    #(power,FDR,FPR) = ibd_admixed2.stats(true_ibd,116430)
    print "IBDAdmixed", score, power_sigs, power, FDR, FPR
