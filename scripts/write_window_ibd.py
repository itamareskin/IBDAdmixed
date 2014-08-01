'''
Created on May 27, 2013

@author: eskin
'''

import os, sys
import random
import string
import logging
from IBD.cIBD import cPairIBD, cPopulationIBD
from itertools import product

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

trueibd_out = open("../data/AfricanAmericans5.trueibd.dat")
ibd_string = trueibd_out.readlines()
trueibd_out.close()
ibd = cPopulationIBD.from_string(ibd_string)

data_f = open("../data/AfricanAmericans5.genos.dat")
data = data_f.readlines()
data_f.close()

ibs_out = open("../data/AfricanAmericans5.ibs.windows.dat", 'w')
num_win = int(1000 / 25)
for pair in ibd.keys():
    ibs_out.write(str(pair[0]) + " " + str(pair[1]))
    for win_idx in range(num_win):
        start_snp = win_idx * 25
        end_snp = min((win_idx + 1) * 25, 1000)
        win_ibs = 0
        if data[pair[0]*2][start_snp:end_snp] == data[pair[1]*2][start_snp:end_snp] or \
        data[pair[0]*2+1][start_snp:end_snp] == data[pair[1]*2][start_snp:end_snp] or \
        data[pair[0]*2][start_snp:end_snp] == data[pair[1]*2+1][start_snp:end_snp] or \
        data[pair[0]*2+1][start_snp:end_snp] == data[pair[1]*2+1][start_snp:end_snp]:
            win_ibs = 1
        ibs_out.write(" " + str(win_ibs))
    ibs_out.write("\n")        
ibs_out.close()  

trueibd_out = open("../data/AfricanAmericans5.trueibd.windows.dat", 'w')
num_win = int(1000 / 25)
for pair in ibd.keys():
    trueibd_out.write(str(pair[0]) + " " + str(pair[1]))
    for win_idx in range(num_win):
        start_snp = win_idx * 25
        end_snp = min((win_idx + 1) * 25, 1000)
        intersect = ibd.get_value(pair).find(start_snp,end_snp)
        win_ibd = 0
        if len(intersect) > 0:
            if (intersect[0].end - intersect[0].start) >= 24:
                win_ibd = 1
        trueibd_out.write(" " + str(win_ibd))
    trueibd_out.write("\n")        
trueibd_out.close()  