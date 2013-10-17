'''
Created on Oct 15, 2013

@author: eskin
'''
from IBD.cIBD import cPairIBD, cPopulationIBD
import pylab as P
import numpy as np
import math
import sys

map_file = sys.argv[1]

true_ibd = cPopulationIBD.fast_deserialize(sys.argv[2], map_file)

gm_f = open(map_file)
data = gm_f.readlines()
dists = [float(x.split(" ")[2]) for x in data]

germline = cPopulationIBD.fast_deserialize(sys.argv[3], map_file)
germline.filter_by_length(0.1,100)
germline.merge_all(25)
(power, FDR, FPR) = germline.stats_win(true_ibd,116430,25)
print power, FDR, FPR