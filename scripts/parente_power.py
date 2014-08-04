'''
Created on Aug 3, 2013

@author: eskin
'''


from IBD.IBDSegments import PairIBD, PopIBD
import pylab as P
import numpy as np
import math
import sys

true_ibd = PopIBD.fast_deserialize(sys.argv[1], sys.argv[3])
parente = PopIBD.fast_deserialize(sys.argv[2], sys.argv[3])
parente.filter_by_human_pairs(true_ibd.keys())
true_ibd.filter_by_human_pairs(parente.keys())
gm_f = open(sys.argv[3])
data = gm_f.readlines()
dists = [float(x.split(" ")[2]) for x in data]

parente.merge_all(max_val = True, merge_diff_vals=False)
parente2 = PopIBD.from_string(parente.to_string(),dists)

(power, FDR, FPR) = parente2.stats_win(true_ibd,116430,25)
print "PARENTE", power, FPR
  
scores = [x[2] for x in parente.to_list()]
print min(scores),max(scores)
for score in range(int(min(scores)),int(max(scores))+1):
    parente2 = PopIBD.from_string(parente.to_string(),dists)
    #parente2.filter_by_human_pairs(true_ibd.keys())
    parente2.filter_by_score(score,1000)
    (power, FDR, FPR) = parente2.stats_win(true_ibd,116430,25)
    print "PARENTE", score, power, FPR