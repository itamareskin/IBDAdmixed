'''
Created on Mar 30, 2013

@author: eskin3
'''

from simuOpt import *
setOptions(optimized=True, alleleType='lineage', version='1.0.1')
from simuPOP import *
import os, sys, string, logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()


pop = loadPopulation("/home/eskin/Data/HapMap/HapMap3_CEU_chr1.pop")
pop.removeIndividuals(range(100,pop.popSize()))
pop1 = loadPopulation("/home/eskin/Data/HapMap/HapMap3_YRI_chr1.pop")
pop1.removeIndividuals(range(100,pop1.popSize()))
pop.addIndFrom(pop1)
map_out = open("/home/eskin/Data/HapMap/HapMap3_CEU_YRI_chr1.map", 'w')
for locus in pop.lociNames():
    pos = '%d' % pop.locusPos(pop.locusByName(locus))
    map_out.writelines(pop.chromNames()[0] + " " + locus + " " + str(pop.dvars().geneticMap[locus]) + " " + pos + "\n")
map_out.close()
ped_out = open("/home/eskin/Data/HapMap/HapMap3_CEU_YRI_chr1.ped", 'w')
count=0
for h1 in pop.individuals():
    count+=1
    if count > 200:
        break
    ped_out.writelines(str(count-1) + " " + str(count-1) + " 0 0 " + str(h1.sex()) + " 1 " + string.join([str(x+1) for t in zip(h1.genotype(0), h1.genotype(1)) for x in t],' ') + "\n")    
ped_out.close()



pop = loadPopulation("/home/eskin/Data/HapMap/HapMap3_CEU_chr1.pop")
map_out = open("/home/eskin/Data/HapMap/HapMap3_CEU_chr1.map", 'w')
for locus in pop.lociNames():
    pos = '%d' % pop.locusPos(pop.locusByName(locus))
    map_out.writelines(pop.chromNames()[0] + " " + locus + " " + str(pop.dvars().geneticMap[locus]) + " " + pos + "\n")
map_out.close()
ped_out = open("/home/eskin/Data/HapMap/HapMap3_CEU_chr1.ped", 'w')
count=0
for h1 in pop.individuals():
    count+=1
    if count > 100:
        break
    ped_out.writelines(str(count-1) + " " + str(count-1) + " 0 0 " + str(h1.sex()) + " 1 " + string.join([str(x+1) for t in zip(h1.genotype(0), h1.genotype(1)) for x in t],' ') + "\n")    
ped_out.close()

pop = loadPopulation("/home/eskin/Data/HapMap/HapMap3_YRI_chr1.pop")
map_out = open("/home/eskin/Data/HapMap/HapMap3_YRI_chr1.map", 'w')
for locus in pop.lociNames():
    pos = '%d' % pop.locusPos(pop.locusByName(locus))
    map_out.writelines(pop.chromNames()[0] + " " + locus + " " + str(pop.dvars().geneticMap[locus]) + " " + pos + "\n")
map_out.close()
ped_out = open("/home/eskin/Data/HapMap/HapMap3_YRI_chr1.ped", 'w')
count=0
for h1 in pop.individuals():
    count+=1
    if count > 100:
        break
    ped_out.writelines(str(count-1) + " " + str(count-1) + " 0 0 " + str(h1.sex()) + " 1 " + string.join([str(x+1) for t in zip(h1.genotype(0), h1.genotype(1)) for x in t],' ') + "\n")    
ped_out.close()



