'''
Created on Mar 16, 2013

@author: eskin3
'''

import simuOpt
simuOpt.setOptions(alleleType='binary', quiet=True, optimized=True)
from simuPOP import *
import random
import string
import logging
from IBD.IBDSegments import PairIBD, PopIBD
from itertools import combinations
import math
import argparse

def add_ibd(pop, ibd, pair, length, dists):
    start = random.randint(50,len(pop.individual(0).lociPos())-1000)
    dists_c = [abs(x-dists[start]-length) for x in dists]
    end = dists_c.index(min(dists_c))
    pop.individual(pair[1]).genotype()[start:end] = pop.individual(pair[0]).genotype()[start:end]
    pairIBD = PairIBD()
    pairIBD.add_interval(start,end)
    ibd.add_human_pair(pair,pairIBD)
    
def create_artificial_pop(pops, alphas, num_composite, num_ibd_pairs, error_rate):
    #haps_mat = np.around(np.random.rand(num_composite*2,len(breaks)-1)*pop.popSize()*2-0.5).astype(int)
    anc_inds = []
    ind = 0
    for alpha in alphas:
        anc_inds += [ind] * (int(round(100 * alpha)))
        ind += 1  
    #genos = [g for t in zip([list(x.genotype(0)) for x in pop.individuals()],[list(x.genotype(1)) for x in pop.individuals()]) for g in t]
    pop_composite = Population(size=num_composite, ploidy=2, loci=[len(pops[0].individual(0).lociPos())],
                     lociPos=pops[0].individual(0).lociPos(), lociNames=pops[0].individual(0).lociNames(), chromNames=pops[0].individual(0).chromNames(),
                     alleleNames=pops[0].individual(0).alleleNames(0), subPopNames=["composite"])
    pop_composite.dvars().geneticMap = pops[0].dvars().geneticMap
    for hap in range(num_composite*2):
        print "simulating haplotype "  + str(hap)
        
        breaks = [0]
        anc_switches = [0]
        last_break = 0
        last_anc_switch = 0
        for i in range(1, len(pops[0].lociNames())):
            if dists[i] - dists[last_break] > 0.2:
                last_break = i;
                breaks.append(i)
            if dists[i] - dists[last_anc_switch] > 4:
                last_anc_switch = i
                anc_switches.append(i)
        if breaks[len(breaks)-1] < len(pops[0].lociNames()) - 1:
            breaks.append(len(pops[0].lociNames()) - 1)
        
        new_hap = []
        anc_ind = random.choice(anc_inds)
        for b in range(len(breaks)-1):
            start = breaks[b]
            if breaks[b] in anc_switches:
                anc_ind = random.choice(anc_inds)
            ind_ind = random.randint(0,pops[anc_ind].popSize()-1)
            hap_ind = random.randint(0,1) 
            new_hap += pops[anc_ind].individual(ind_ind).genotype(hap_ind)[breaks[b]:breaks[b+1]]
            #new_hap += genos[haps_mat[hap][b]][breaks[b]:breaks[b+1]]
        pop_composite.individual(int(hap/2)).setGenotype(new_hap, hap%2)
      
    save_pop(pop_composite,scramble=args.scramble, suffix=".before")
    pairs = list(combinations(range(pop_composite.popSize()),2))
    random.shuffle(pairs)
    pairs = set(pairs)
    ibd = PopIBD()
    it = iter(pairs)
    for ibd_length in [1,2,3,4,5]:
        for i in range(num_ibd_pairs):
            print "simulating ibd pair " + str(i) + " with ibd length " + str(ibd_length) + " cM"
            pair = next(it)
            add_ibd(pop_composite, ibd, pair, ibd_length, dists)
    
    pop_composite.dvars().ibd = ibd
    
    if error_rate > 0:        
        error_num = int(len(pop_composite.lociNames())*error_rate)
        for i in range(pop_composite.popSize()):
            print "adding errors to individual " + str(i)
            for c in range(2):
                k= max(0, int(random.gauss(error_num, math.sqrt(error_num))))
                error_positions = random.sample(range(len(pop_composite.lociNames())),k)
                for error_pos in error_positions:
                    pop_composite.individual(i).genotype(c)[error_pos] = 1 - pop_composite.individual(i).genotype(c)[error_pos]
    
    return pop_composite

def save_pop(pop,scramble=False, suffix=""):
    pop.save(args.out + ".pop")

    if "ibd" in pop.vars():
        with open(args.out + suffix + ".trueibd.txt" ,"w") as f:
            f.write(pop.dvars().ibd.to_string())

        with open(args.out + suffix + ".trueibd.pairs.txt" ,"w") as f:
            f.write(string.join([str(x[0]) + "," + str(x[1]) for x in pop.dvars().ibd.keys()],"\n"))

    map_out = open(args.out + suffix + ".genos.map", 'w')
    for locus in pop.lociNames():
        pos = '%d' % pop.locusPos(pop.locusByName(locus))
        map_out.writelines(pop.chromNames()[0] + " " + locus + " " + str(pop.dvars().geneticMap[locus]) + " " + pos + "\n")
    map_out.close()
    
    if scramble:
        for h1 in pop.individuals():
            #logger.info("Writing data of individual %d  ", int(h1.info('ind_id')))
            scrambled = zip(*[(x[1],x[0]) if random.randint(0,1) == 1 else x for x in zip(h1.genotype(0),h1.genotype(1))])
            h1.setGenotype(list(scrambled[0]), 0)
            h1.setGenotype(list(scrambled[1]), 1)
    
    data_out = open(args.out + suffix + ".genos.dat", 'w')
    count=0
    for h1 in pop.individuals():
        count+=1
        #logger.info("Writing data of individual %d  ", int(h1.info('ind_id')))
        data_out.writelines(string.join([str(x) for x in h1.genotype(0)],'') + "\n")
        data_out.writelines(string.join([str(x) for x in h1.genotype(1)],'') + "\n")
    data_out.close()
    
    ped_out = open(args.out + suffix + ".genos.ped", 'w')
    count=0
    for h1 in pop.individuals():
        count+=1
        #logger.info("Writing data of individual %d  ", int(h1.info('ind_id')))
        ped_out.writelines(str(count-1) + " " + str(count-1) + " 0 0 " + str(h1.sex()) + " 1 " + string.join([str(x+1) for t in zip(h1.genotype(0), h1.genotype(1)) for x in t],' ') + "\n")    
    ped_out.close()

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

parser = argparse.ArgumentParser()
parser.add_argument("popfiles", nargs='*', type=str, help=".pop file name (one for each ancestry)")
parser.add_argument("map", type=str, help="map file")
parser.add_argument("out", type=str, help="output prefix")
parser.add_argument("-a", "--set-alphas", nargs='+', dest='alphas', required=True, default=[],help="set the alphas")
parser.add_argument("-n", "--num-inds", action="store", dest='num_inds', required=True, help="number of individuals to simulate")
parser.add_argument("-i", "--ibd-pairs", action="store", dest='ibd_pairs', required=True, help="number of ibd pairs to simulate")
parser.add_argument("-e", "--error-rate", action="store", dest='error_rate', required=True, help="genotyping error rate to add")
parser.add_argument("-s", "--scramble", action='store_true', default=False, dest='scramble', help='scramble phase of genotypes')


args = parser.parse_args()

logf = open(args.out + ".log", 'w')
logf.write(str(args) + "\n")
logf.close()
    
if args.alphas != None:
    if len(args.alphas) > 0:
        alphas = [float(x) for x in args.alphas]

pop_files = args.popfiles
pops = []
for pop_file in pop_files:
    print "loading population: " + pop_file
    pops.append(loadPopulation(pop_file))
    
num_inds = int(args.num_inds)
ibd_pairs = int(args.ibd_pairs)
error_rate = float(args.error_rate)
        
# train_size = 100
# test_size = pop.popSize() - train_size
# pop_train_base = pop.clone()
# pop_test_base = pop.clone()
# pop_train_base.removeIndividuals(range(15,113)+range(169,226))
# pop_test_base.removeIndividuals(range(0,15)+range(113,169))

#m = LDModel(sys.argv[2] + ".genos.map",log_dir = ".",k = 1,g = 8,max_snp_num = 100000000,debug = False)
#m.generate_composite_individuals(5) 

dists = []
gm_f = open(args.map)
data = gm_f.readlines()
dists = [float(x.split(" ")[2]) for x in data]

pop_out = create_artificial_pop(pops,alphas,num_inds,ibd_pairs,error_rate)
print "saving simulated population to file..."
save_pop(pop_out,scramble=args.scramble)

print "Finished!"

x = 1