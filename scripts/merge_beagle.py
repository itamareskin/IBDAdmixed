'''
Created on Oct 8, 2013

@author: eskin
'''

from IBD.cIBD import cPairIBD, cPopulationIBD
import pylab as P
import numpy as np
import math
import sys


true_ibd = cPopulationIBD.fast_deserialize(sys.argv[1], sys.argv[3])
beagle = cPopulationIBD.fast_deserialize(sys.argv[2], sys.argv[3])
beagle.merge_all()
gm_f = open(sys.argv[3])
data = gm_f.readlines()
dists = [float(x.split(" ")[2]) for x in data]
for score in range(0,40):
    beagle2 = cPopulationIBD.from_string(beagle.to_string(),dists)
    #beagle2 = cPopulationIBD.fast_deserialize(sys.argv[2], sys.argv[3])
    #beagle.filter_by_length(2,1000)
    beagle2.filter_by_score(0,math.pow(10, -score))
    (power, FDR, FPR) = beagle2.stats_win(true_ibd,116430,25)
    print "fastIBD", math.pow(10, -score), power, FPR
