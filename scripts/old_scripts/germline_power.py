'''
Created on Oct 15, 2013

@author: eskin
'''
from IBD.IBDSegments import PairIBD, PopIBD
import pylab as P
import numpy as np
import math
import sys
from IBD.GeneticMap import GeneticMap

map_file = sys.argv[1]

true_ibd = PopIBD.fast_deserialize(sys.argv[2])

gm = GeneticMap(map_file)

germline = PopIBD.fast_deserialize(sys.argv[3])
germline.filter_by_length(0.1,100,gm)
germline.merge_all(25)
(power, FDR, FPR) = germline.stats_win(true_ibd,116430,25)
print power, FDR, FPR