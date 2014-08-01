'''
Created on Oct 8, 2013

@author: eskin
'''

from IBD.cIBD import cPairIBD, cPopulationIBD
import pylab as P
import numpy as np
import math
import sys
from IBD.GeneticMap import GeneticMap

true_ibd = cPopulationIBD.fast_deserialize(sys.argv[1])
#true_ibd.filter_by_length(4.5,10)
beagle = cPopulationIBD.fast_deserialize(sys.argv[2])
beagle.filter_by_human_pairs(true_ibd.keys())
true_ibd.filter_by_human_pairs(beagle.keys())
beagle.merge_all()

gm = GeneticMap(sys.argv[3])

beagle2 = cPopulationIBD.from_string(beagle.to_string())
beagle2.filter_by_score(0,1e-5)
beagle2.merge_all(merge_diff_vals=True)
beagle2.filter_by_length(0.5,1000,gm)
f = open(sys.argv[2] + ".filt.ibd.txt","w")
f.write(beagle2.to_string())
f.close()

for score in range(0,80):
    beagle2 = cPopulationIBD.from_string(beagle.to_string())
    beagle2.filter_by_score(0,math.pow(10, -score))
    beagle2.merge_all(merge_diff_vals=True)
    beagle2.filter_by_length(0.8,1000,gm)

    (power, FDR, FPR) = beagle2.stats_win(true_ibd,116430,25)
    print "fastIBD", math.pow(10, -score), power, FDR, FPR
