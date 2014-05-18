'''
Created on Mar 16, 2013

@author: eskin3
'''

import simuOpt
simuOpt.setOptions(alleleType='binary', quiet=True, optimized=True)
import simuPOP as sim
from simuPOP import *
import os, sys
import random
import string
import logging
from utils import FoundersContainer.FoundersContainer
from IBD.cIBD import cPairIBD, cPopulationIBD 
from itertools import product, combinations 
from hapMapUtil import getHapMapMarkers
from sys import stdout
from IBD.LDModel import LDModel
import math
import numpy as np

def add_ibd(pop, ibd, pair, length, dists):
    start = random.randint(50,100000)
    dists_c = [abs(x-dists[start]-length) for x in dists]
    end = dists_c.index(min(dists_c))
    pop.individual(pair[1]).genotype()[start:end] = pop.individual(pair[0]).genotype()[start:end]
    pairIBD = cPairIBD()
    pairIBD.add_interval(start,end)
    ibd.add_human_pair(pair,pairIBD)
    
def create_artificial_pop(pop, num_composite, num_ibd_pairs, error_rate):
    haps_mat = np.around(np.random.rand(num_composite*2,len(breaks)-1)*pop.popSize()*2-0.5).astype(int)
    genos = [g for t in zip([list(x.genotype(0)) for x in pop.individuals()],[list(x.genotype(1)) for x in pop.individuals()]) for g in t]
    pop_composite = Population(size=num_composite, ploidy=2, loci=[len(pop.individual(0).lociPos())],
                     lociPos=pop.individual(0).lociPos(), lociNames=pop.individual(0).lociNames(), chromNames=pop.individual(0).chromNames(),
                     alleleNames=pop.individual(0).alleleNames(0), subPopNames=["composite"])
    pop_composite.dvars().geneticMap = pop.dvars().geneticMap
    for hap in range(num_composite*2):
        new_hap = []
        for b in range(len(breaks)-1):
            start = breaks[b]
            new_hap += genos[haps_mat[hap][b]][breaks[b]:breaks[b+1]]
        pop_composite.individual(int(hap/2)).setGenotype(new_hap, hap%2)
      
    pairs = list(combinations(range(pop_composite.popSize()),2))
    random.shuffle(pairs)
    pairs = set(pairs)
    ibd = cPopulationIBD()
    it = iter(pairs)
    for ibd_length in [1,2,3,4,5]:
        for i in range(num_ibd_pairs):
            pair = next(it)
            add_ibd(pop_composite, ibd, pair, ibd_length, dists)
    
    pop_composite.dvars().ibd = ibd
            
    error_num = int(len(pop_composite.lociNames())*error_rate)
    for i in range(pop_composite.popSize()):
        for c in range(2):
            k= max(0, int(random.gauss(error_num, math.sqrt(error_num))))
            error_positions = random.sample(range(len(pop_composite.lociNames())),k)
            for error_pos in error_positions:
                pop_composite.individual(i).genotype(c)[error_pos] = 1 - pop_composite.individual(i).genotype(c)[error_pos]
    
    return pop_composite

def save_pop(pop,prefix):
    pop.save(sys.argv[2] + ".pop")

    f = open(sys.argv[2] + "." + prefix + ".trueibd.txt","w")
    f.write(pop.dvars().ibd.to_string())
    f.close()
    
    map_out = open(sys.argv[2] + "." + prefix + ".genos.map", 'w')
    for locus in pop.lociNames():
        pos = '%d' % pop.locusPos(pop.locusByName(locus))
        map_out.writelines(pop.chromNames()[0] + " " + locus + " " + str(pop.dvars().geneticMap[locus]) + " " + pos + "\n")
    map_out.close()
    
    data_out = open(sys.argv[2] + "." + prefix + ".genos.dat", 'w')
    count=0
    for h1 in pop.individuals():
        count+=1
        #logger.info("Writing data of individual %d  ", int(h1.info('ind_id')))
        data_out.writelines(string.join([str(x) for x in h1.genotype(0)],'') + "\n")
        data_out.writelines(string.join([str(x) for x in h1.genotype(1)],'') + "\n")
    data_out.close()
    
    ped_out = open(sys.argv[2] + "." + prefix + ".genos.ped", 'w')
    count=0
    for h1 in pop.individuals():
        count+=1
        #logger.info("Writing data of individual %d  ", int(h1.info('ind_id')))
        ped_out.writelines(str(count-1) + " " + str(count-1) + " 0 0 " + str(h1.sex()) + " 1 " + string.join([str(x+1) for t in zip(h1.genotype(0), h1.genotype(1)) for x in t],' ') + "\n")    
    ped_out.close()

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

pop_files = string.split(sys.argv[1],",")
first = True
pop = None
for pop_file in pop_files:
    if first:
        pop = loadPopulation(pop_file)
        first = False;
    else:
        pop.addIndFrom(loadPopulation(pop_file))
        
train_size = 100
test_size = pop.popSize() - train_size
pop_train_base = pop.clone()
pop_test_base = pop.clone()
pop_train_base.removeIndividuals(range(train_size,pop.popSize()))
pop_test_base.removeIndividuals(range(train_size))

#m = LDModel(sys.argv[2] + ".genos.map",log_dir = ".",k = 1,g = 8,max_snp_num = 100000000,debug = False)
#m.generate_composite_individuals(5) 

dists = []
gm_f = open(sys.argv[3])
data = gm_f.readlines()
dists = [float(x.split(" ")[2]) for x in data]

breaks = [0]
last_break = 0
for i in range(1, len(pop.lociNames())):
    if dists[i] - dists[last_break] > 0.2:
        last_break = i;
        breaks.append(i)

if breaks[len(breaks)-1] < len(pop.lociNames()) - 1:
    breaks.append(len(pop.lociNames()) - 1)

pop_train = create_artificial_pop(pop_train_base,250,0,0)
pop_test = create_artificial_pop(pop_train_base,500,80,0.005)
        
save_pop(pop_train,"train")
save_pop(pop_test,"test")


x = 1