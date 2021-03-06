'''
Created on Aug 3, 2013

@author: eskin
'''


from IBD.IBDSegments import PairIBD, PopIBD
import pylab as P
import numpy as np
import math
from IBD.GeneticMap import GeneticMap

map_file = "/home/eskin/Data/IBDAdmixed/fromLecs/HapMap3_CEU_chr2.map"

true_ibd = PopIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/New/single.ceu.trueibd.txt")
#true_ibd = PopIBD.fast_deserialize("/home/eskin/Data/beagle/fromLecs/artificial.admixed.test.trueibd.txt", map_file)
#true_ibd.filter_by_length(2.5,10)

gm = GeneticMap(map_file)

beagle = PopIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/New/single.ceu.beagle4.ibd.txt")
beagle.filter_by_human_pairs(true_ibd.keys())
true_ibd.filter_by_human_pairs(beagle.keys())

beagle2 = PopIBD.from_string(beagle.to_string())
(power, FDR, FPR) = beagle2.stats_win(true_ibd,116430,25)
print power, FDR, FPR

scores = [int(x[2]) for x in beagle.to_list()]
print min(scores),max(scores)
beagle2.merge_all(overlap = 1, max_val = True, merge_diff_vals=True)
beagle2.filter_by_length(0.1,1000,gm)
(power, FDR, FPR) = beagle2.stats_win(true_ibd,116430,25)
print power, FDR, FPR
 
for score in range(min(scores),max(scores)+1):
    beagle2 = PopIBD.from_string(beagle.to_string())
    beagle2.filter_by_score(score,205)
    beagle2.merge_all(overlap = 1, max_val = True, merge_diff_vals=True)
    beagle2.filter_by_length(0.5,1000,gm)
    (power, FDR, FPR) = beagle2.stats_win(true_ibd,116430,25)
    print "beagle4", score, power, FPR
