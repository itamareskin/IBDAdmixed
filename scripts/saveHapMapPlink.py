'''
Created on Mar 30, 2013

@author: eskin3
'''

from simuOpt import *
setOptions(optimized=True, alleleType='lineage', version='1.0.1')
from simuPOP import *
import os, sys, string, logging
import random

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()


# dir = sys.argv[1]

# pop = loadPopulation(dir + "/HapMap3_CEU_chr1.pop")
# #pop.removeIndividuals(range(100,pop.popSize()))
# pop1 = loadPopulation(dir + "/HapMap3_TSI_chr2.pop")
# #pop1.removeIndividuals(range(100,pop1.popSize()))
# pop.addIndFrom(pop1)
# shuffeled_indexes = range(pop.popSize())
# random.shuffle(shuffeled_indexes)
# pop.addInfoFields("index")
# i=0
# for ind in pop.individuals():
#     ind.index = shuffeled_indexes[i]
#     i+=1
# pop.sortIndividuals("index") 
# 
# pop2 = pop.clone()
# pop.removeIndividuals(range(100))
# pop2.removeIndividuals(range(100,201))
# pop2.save(dir + "/HapMap3_CEU_TSI_chr2.train.pop")
# pop.save(dir + "/HapMap3_CEU_TSI_chr2.test.pop")

pop = loadPopulation(sys.argv[1])

path = sys.argv[1]
(dir,file) = os.path.split(path)
(file_prefix,ext) = os.path.splitext(file)

map_out = open(os.path.join([dir,string.join([file_prefix,"map"], ".")]), 'w')
for locus in pop.lociNames():
    pos = '%d' % pop.locusPos(pop.locusByName(locus))
    map_out.writelines(pop.chromNames()[0] + " " + locus + " " + str(pop.dvars().geneticMap[locus]) + " " + pos + "\n")
map_out.close()
ped_out = open(os.path.join([dir,string.join([file_prefix,"ped"], ".")]), 'w')
count=0
for h1 in pop.individuals():
    count+=1
    #if count > 200:
    #    break
    ped_out.writelines(str(count-1) + " " + str(count-1) + " 0 0 " + str(h1.sex()) + " 1 " + string.join([str(x+1) for t in zip(h1.genotype(0), h1.genotype(1)) for x in t],' ') + "\n")
    #ped_out.writelines(str(count-1) + " " + str(count-1) + " 0 0 " + str(h1.sex()) + " 1 " + string.join([h1.],' ') + "\n")
ped_out.close()


