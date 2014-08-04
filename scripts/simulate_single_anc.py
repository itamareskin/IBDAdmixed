'''
Created on Mar 16, 2013

@author: eskin3
'''

import simuOpt
simuOpt.setOptions(alleleType='lineage', quiet=True, optimized=True)
import simuPOP as sim
from simuPOP import *
import os, sys
import random
import string
import logging
from utils import FoundersContainer.FoundersContainer
from IBD.IBDSegments import PairIBD, PopIBD
from itertools import product
from hapMapUtil import getHapMapMarkers

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

'''
custom population size function
'''


def demo(gen, pop):
    if gen<17 and pop.popSize() < 10000:
        return int(pop.popSize()*random.randint(130, 140)/100.0) 
    else:
        return int(pop.popSize())

if __name__ == '__main__':
    
    if len(sys.argv) != 3:
        print "usage: " + sys.argv[0] + " <filename> <num_snps>";
        exit(-1);
    
    filename = sys.argv[1]
    num_snps = int(sys.argv[2])
    
    if not os.path.isfile(filename + ".pop"):
    # Get the hapmap population samples
        pop = getHapMapMarkers('~/Data/IBDAdmixed', 
            names = [],
            chroms=[1], 
            HapMap_pops=['CEU'],
            startPos = [],
            endPos = [],
            numMarkers = [num_snps],
            minAF = 0,
            minDist = 0,
            logger=logger)
        
        #pop = sim.Population(size=(50,50), loci=[10]*4, )
            
        # define two virtual subpopulations by ancestry value
        # pop.setVirtualSplitter(sim.InfoSplitter(field='ancestry', cutoff = [0.5]))
        pop.addInfoFields(['ind_id'])
        sim.initInfo(pop,values=range(pop.popSize()),infoFields='ind_id')
        
        transmitters=[
            #sim.MendelianGenoTransmitter(),
            sim.Recombinator(intensity=1e-9),
            sim.InheritTagger(mode=sim.MEAN)]
        pop.evolve(
            initOps=[sim.InitSex(),
                     sim.InitLineage(mode=sim.PER_PLOIDY),
                     ],
            preOps=[sim.Stat(popSize=True),
                    #sim.Stat(effectiveSize=range(3), subPops=[0, 1], vars='Ne_demo_base_sp'),
                    #sim.Stat(inbreeding=sim.ALL_AVAIL, popSize=True, step=1),
                    #sim.PyEval(r'"gen %d: IBD freq %.4f, IBS freq %.4f, est: %.4f\n" % '
                    #'(gen, sum(IBD_freq.values()) /len(IBD_freq), '
                    #' sum(IBS_freq.values()) /len(IBS_freq), '
                    #' 1 - (1-1/(2.*popSize))**gen)', step=1),
                    sim.PyEval(r'"Gen %d:\t%s\n" % (gen, subPopSize)')],
            matingScheme=sim.RandomMating(ops=transmitters,subPopSize=demo),
            postOps=[sim.Stat(popSize=True),
                     #sim.Stat(effectiveSize=range(3), subPops=[0, 1], vars='Ne_demo_sp'),        
                     sim.PyEval(r'"%s\n" % subPopSize'),
                   #sim.PyEval(r'"Euro Ne: %.1f Afr Ne: %.1f \n" % (subPop[0]["Ne_demo"][0],subPop[1]["Ne_demo"][0])')
                   ],
            gen=25,
        )   
    
        #pop.addInfoFields(['ind_id'])
        sim.initInfo(pop,values=range(pop.popSize()),infoFields='ind_id')
        
        pop.removeIndividuals(range(100,pop.popSize()))
        
        pop.save(filename + ".pop")
    
    else:
        pop = loadPopulation(filename + ".pop")
    
    if True:
        map_out = open(filename + ".genos.map", 'w')
        for locus in pop.lociNames():
            pos = '%d' % pop.locusPos(pop.locusByName(locus))
            map_out.writelines(pop.chromNames()[0] + " " + locus + " " + str(pop.dvars().geneticMap[locus]) + " " + pos + "\n")
        map_out.close()
        
        data_out = open(filename + ".genos.dat", 'w')
        count=0
        for h1 in pop.individuals():
            count+=1
            if count > 100:
                break
            logger.info("Writing data of individual %d  ", int(h1.info('ind_id')))
            data_out.writelines(string.join([str(x) for x in h1.genotype(0)],'') + "\n")
            data_out.writelines(string.join([str(x) for x in h1.genotype(1)],'') + "\n")
        data_out.close()
        
        ped_out = open(filename + ".genos.ped", 'w')
        count=0
        for h1 in pop.individuals():
            count+=1
            if count > 100:
                break
            logger.info("Writing data of individual %d  ", int(h1.info('ind_id')))
            ped_out.writelines(str(count-1) + " " + str(count-1) + " 0 0 " + str(h1.sex()) + " 1 " + string.join([str(x+1) for t in zip(h1.genotype(0), h1.genotype(1)) for x in t],' ') + "\n")    
        ped_out.close()
        
        
        if True:
            lin_out = open(filename + ".lineages.dat", 'w')
            for h1 in pop.individuals():
                logger.info("Writing lineage of individual %d  ", int(h1.info('ind_id')))
                lin_out.writelines(string.join([str(x) for x in h1.lineage(0)]," ") + "\n")
                lin_out.writelines(string.join([str(x) for x in h1.lineage(1)]," ")+ "\n")
            lin_out.close()
    
    #pop.addInfoFields(['contributors','contributor_container'])
    # initialize ancestry
    #sim.initInfo(pop, [set(ind.lineage()) for ind in pop.individuals()],
    #infoFields='contributors')
    #sim.initInfo(pop, [FoundersContainer.from_founders_list(list(ind.lineage()),pop.locusPos(0),pop.locusPos(pop.totNumLoci)) for ind in pop.individuals()],
    #infoFields='contributor_container')
    
    cont_to_inds = {}
    ind_to_cont_list = {}
    ind_to_cont_struct = {}
    count=0
    for h1 in pop.individuals():
        count+=1
        if count > 100:
            break
        ind_to_cont_list[int(h1.info('ind_id'))] = set(h1.lineage())
        ind_to_cont_struct[int(h1.info('ind_id'))] = \
        [FoundersContainer.from_founders_list(list(h1.lineage(ploidy=0)),0,pop.totNumLoci()), 
         FoundersContainer.from_founders_list(list(h1.lineage(ploidy=1)),0,pop.totNumLoci())] 
        for cont in ind_to_cont_list[int(h1.info('ind_id'))]:
            if not cont_to_inds.has_key(cont):
                cont_to_inds[cont] = []
            if int(h1.info('ind_id')) not in cont_to_inds[cont]:  
                cont_to_inds[cont].append(int(h1.info('ind_id')))
            
    num_ibd_segments = 0
    total_ibd_length = 0
    count = 0
    ibd = PopIBD()
    ibd_in_admixed = PopIBD()
    for h1 in pop.individuals():
        count+=1
        if count > 100:
            break
            
        humans = []
        for cont in ind_to_cont_list[int(h1.info('ind_id'))]:
            for ind in cont_to_inds[cont]:
                if ind not in [i.info('ind_id') for i in humans]:
                    humans.append(pop.individual(ind))
                
        if logger is not None:
            if len(humans) > 0:
                logger.info("Calculating IBD for human with ID: %s (%d out of %d) with %d other humans ", int(h1.info('ind_id')), count, pop.popSize(), len(humans))
            
        for h2 in humans:
            if int(h1.info('ind_id')) >= int(h2.info('ind_id')):
                continue

            pairIBD = PairIBD()
            pairIBD_admixed = PairIBD()
            for c1,c2 in [(0,0)]: #product([0,1],[0,1]):
                pairIBD.add_ibd_from_containers(ind_to_cont_struct[int(h1.info('ind_id'))][c1],ind_to_cont_struct[int(h2.info('ind_id'))][c2],0)
                #pairIBD_admixed.add_ibd_from_containers_only_admixed(ind_to_cont_struct[int(h1.info('ind_id'))][c1],ind_to_cont_struct[int(h2.info('ind_id'))][c2],0,ind_to_cont_struct[int(h1.info('ind_id'))][1-c1],ind_to_cont_struct[int(h2.info('ind_id'))][1-c2], ceu_inds, yri_inds)
                #pairIBD.add_ibd_from_containers(ind_to_cont_struct[int(h1.info('ind_id'))][c1],ind_to_cont_struct[int(h2.info('ind_id'))][c2],0,ceu_inds,yri_inds)
                
            if pairIBD != None:
                if len(pairIBD.to_list()) > 0:
                    ibd.add_human_pair((int(h1.info('ind_id')),int(h2.info('ind_id'))), pairIBD)
            
            #if pairIBD_admixed != None:
            #    if len(pairIBD_admixed.to_list()) > 0:
            #        ibd_in_admixed.add_human_pair((int(h1.info('ind_id')),int(h2.info('ind_id'))), pairIBD_admixed)
                    
    
    trueibd_out = open(filename + ".trueibd.dat", 'w')
    trueibd_out.writelines(ibd.to_string())
    trueibd_out.close()
    
#    trueibd_out = open(filename + ".trueibd.admixed.dat", 'w')
#    trueibd_out.writelines(ibd_in_admixed.to_string())
#    trueibd_out.close()
    
#    pairs_f = open(filename + ".pairs.dat")
#    pairs = pairs_f.readlines()
#    pairs = [x.strip("\n") for x in pairs]
#    pairs = [x.split(" ") for x in pairs]
#    pairs = [(int(x[0]),int(x[1])) for x in pairs]
#    pairs_f.close()
    
    trueibd_out = open(filename + ".trueibd.windows.dat", 'w')
    num_win = int(num_snps / 25)
    for pair in ibd.keys():
        trueibd_out.write(str(pair[0]) + " " + str(pair[1]))
        for win_idx in range(num_win):
            start_snp = win_idx * 25
            end_snp = min((win_idx + 1) * 25, num_snps)
            if ibd.has_key(pair):
                intersect = ibd.get_value(pair).find(start_snp,end_snp)
            else:
                intersect = []
            win_ibd = 0
            if len(intersect) > 0:
                if (intersect[0].end - intersect[0].start) >= 24:
                    win_ibd = 1
            trueibd_out.write(" " + str(win_ibd))
        trueibd_out.write("\n")        
    trueibd_out.close()
    
    ibs_out = open(filename + ".ibs.windows.dat", 'w')
    for pair in ibd.keys():
        ibs_out.write(str(pair[0]) + " " + str(pair[1]))
        for win_idx in range(num_win):
            start_snp = win_idx * 25
            end_snp = min((win_idx + 1) * 25, num_snps)
            win_ibs = 0
            if pop.individual(pair[0]).genotype(0)[start_snp:end_snp] == pop.individual(pair[1]).genotype(0)[start_snp:end_snp] or \
            pop.individual(pair[0]).genotype(1)[start_snp:end_snp] == pop.individual(pair[1]).genotype(0)[start_snp:end_snp] or \
            pop.individual(pair[0]).genotype(0)[start_snp:end_snp] == pop.individual(pair[1]).genotype(1)[start_snp:end_snp] or \
            pop.individual(pair[0]).genotype(1)[start_snp:end_snp] == pop.individual(pair[1]).genotype(1)[start_snp:end_snp]:
                win_ibs = 1
            ibs_out.write(" " + str(win_ibs))
        ibs_out.write("\n")        
    ibs_out.close()      
    
    x = 1
    #export(pop, format='csv', output='AfricanAmericans.genos.dat', delimiter='', header=False, infoFields=[],sexFormatter=None, affectionFormatter=None, gui=None)
    #export(pop, format='PED', output='AfricanAmericans.genos.ped', gui=None)