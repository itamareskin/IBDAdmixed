'''
Created on Mar 16, 2013

@author: eskin3
'''

import simuOpt
simuOpt.setOptions(alleleType='lineage', quiet=True, optimized=True)
import simuPOP as sim
from simuPOP import *
import loadHapMap
import os, sys, urllib, gzip, exceptions, tempfile, shutil

#pop = sim.Population(1000, loci=[100]*1)
#
#pop.evolve(
#     initOps=[
#         sim.InitSex(),
#         sim.InitGenotype(freq=[0.5]*2),
#         sim.InitLineage(range(1000), mode=sim.PER_INDIVIDUAL),
#     ],
#     matingScheme=sim.RandomMating(ops=sim.Recombinator(rates=0.001)),
#     gen = 15
# )
#
#lineage = pop.lineage()
#x = 1

pop = createHapMapPop()
for chrom in [22]:
    for sample in ['CEU', 'YRI']:
        popFile = os.path.join("/home/eskin3/workspace/IBDAdmixed/scripts", "HapMap_%s_chr%d.pop" % (sample, chrom))
        try:
            if os.path.isfile(popFile):
                # test if this file is OK.
                loadHapMapPop(pop,chrom,popFile)
                pop.removeLoci(keep=range(100))
                continue
        except e:
            # continue to load file
            pass

pop.addInfoFields(['ancestry', 'migrate_to'])
# initialize ancestry
sim.initInfo(pop, [0]*pop.subPopSize(0) + [1]*pop.subPopSize(1),
    infoFields='ancestry')
# define two virtual subpopulations by ancestry value
pop.setVirtualSplitter(sim.InfoSplitter(field='ancestry', cutoff = [0.5]))
transmitters=[
    sim.MendelianGenoTransmitter(),
    sim.InheritTagger(mode=sim.MEAN, infoFields='ancestry')]
pop.evolve(
    initOps=sim.InitSex(),
    preOps=sim.Migrator(rate=[
        [0., 0], [0.05, 0]]), 
    matingScheme=sim.HeteroMating(
        matingSchemes=[
            sim.RandomMating(ops=transmitters),
            sim.RandomMating(subPops=[(0,0)], weight=-0.80, ops=transmitters),
            sim.RandomMating(subPops=[(0,1)], weight=-0.80, ops=transmitters)
        ],
    ),
    gen=10,
)
# remove the second subpop
pop.removeSubPops(1)