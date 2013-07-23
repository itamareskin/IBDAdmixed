'''
Created on May 27, 2013

@author: eskin
'''

import os, sys
import random
import string
import logging
from IBD.FoundersContainer import FoundersContainer
from IBD.cIBD import cPairIBD, cPopulationIBD 
from itertools import product

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

if len(sys.argv) != 4:
    print "usage: " + sys.argv[0] + " <num_snps> <inputfile> <outputfile>"
    exit(-1) 

num_snps = int(sys.argv[1])
ibd_f = open(sys.argv[2])
ibd_string = ibd_f.readlines()
ibd_f.close()
ibd = cPopulationIBD.from_string(ibd_string)

ibd_win_f = open(sys.argv[3], 'w')
num_win = int(num_snps / 25)
for pair in ibd.keys():
    ibd_win_f.write(str(pair[0]) + " " + str(pair[1]))
    for win_idx in range(num_win):
        start_snp = win_idx * 25
        end_snp = min((win_idx + 1) * 25, num_snps)
        intersect = ibd.get_value(pair).find(start_snp,end_snp)
        win_ibd = 0
        if len(intersect) > 0:
            if (intersect[0].end - intersect[0].start) >= 24:
                win_ibd = 1
        ibd_win_f.write(" " + str(win_ibd))
    ibd_win_f.write("\n")        
ibd_win_f.close()  