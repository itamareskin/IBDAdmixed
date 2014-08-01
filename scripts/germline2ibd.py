__author__ = 'Itamar'

import sys
import os
from IBD.GeneticMap import GeneticMap
from IBD.cIBD import cPopulationIBD
from itertools import islice

def intervals_from_germline_file(germlinefile, pairs, pos_dict):
    ibs_intervals = {}
    buffer_size = 5000
    print "reading IBS intervals from GERMLINE file: " + germlinefile
    with open(germlinefile) as germline:
        done = False
        while True:
            if done:
                break
            lines = list(islice(germline, buffer_size))
            if len(lines) == 0:
                done = True
            for line in lines:
                if not line:
                    break
                line = line.replace('\t', ' ')
                line = line.split(' ')
                pair = (int(line[0]), int(line[2]))
                if pair not in pairs:
                    continue
                if not ibs_intervals.has_key(pair):
                    ibs_intervals[pair] = []
                if pos_dict.has_key(long(line[5])) and pos_dict.has_key(long(line[6])):
                    ibs_intervals[pair].append((pos_dict[long(line[5])], pos_dict[long(line[6])]))
    print "finished reading IBS intervals."
    return ibs_intervals

gm = GeneticMap("K:\Data\IBDAdmixed\New4\ceu.tsi.yri.lwk.half2.genos.map", 116430)
pos_dict = gm.get_position_dict()

with open("K:\Data\IBDAdmixed\New4\ceu.tsi.yri.lwk.half2.trueibd.pairs.txt") as pairs_f:
    pairs = pairs_f.readlines()
    pairs = [x.strip("\n") for x in pairs]
    pairs = [x.split(",") for x in pairs]
    pairs = [(int(x[0]), int(x[1])) for x in pairs]

intervals = intervals_from_germline_file("K:\Data\IBDAdmixed\New4\ceu.tsi.yri.lwk.half2.match", pairs, pos_dict)
ibd = cPopulationIBD.from_dict(intervals)
with open("K:\Data\IBDAdmixed\New4\ceu.tsi.yri.lwk.half2.germline.ibd.txt","w") as outf:
    outf.write(ibd.to_string())
