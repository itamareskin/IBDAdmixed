#!/usr/bin/env python

from simuOpt import *
setOptions(optimized=True, alleleType='binary', version='1.0.1')
from simuPOP import *

import os, sys, urllib, gzip, tempfile, time


def _probeInfo(datafile, logger=None):
    '''Get population size of a sample
      This function also checks if each line has desired number of alleles.
    '''
    count = 0
    ll = 0
    data = open(datafile)
    
    chromNames = []
    name = []
    alleleNames = {}
    pos = []
    line_no = 0
    for line in data.readlines():
        data = line.split()
        chromNames.append(data[0])
        name.append(data[1])
        pos.append(int(data[3]))
        alleles = list(set(data[4:]))
        alleles.sort()
        alleleNames[data[1]] = alleles
        line_no+=1
    numInds = (len(data) - 4) / 2
    return numInds, name, alleleNames, pos, chromNames
    

def load_population(pop, diskfile, alleleNames, logger=None):
    '''Load population from file, with type (subpopulation type)'''
    # file format:
    #
    # rsID pos ind1_A ind2_A ....
    #
    data = open(diskfile)
    for line_no,line in enumerate(data.readlines()):
        fields = line.split()
        name = fields[1]
        alleleName = alleleNames[name]
        genotype = [alleleName.index(x) for x in fields[4:]]
        for col_no,geno in enumerate(genotype):
            ind = col_no / 2
            ploidy = col_no % 2
            if ploidy == 2:
                continue
            if ind == pop.popSize():
                print 'Warning: individual index %d greater than population size %d ' % (ind, pop.popSize())
            # always chromosome 0, because each population has only one chromosome
            assert pop.locusPos(line_no) == float(fields[3])
            assert pop.locusName(line_no) == name
            pop.individual(ind).setAllele(geno, line_no, ploidy)

if __name__ == '__main__':

    import logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('loadPlink')

    diskfile = sys.argv[1]
    (numInds, lociNames, alleleNames, lociPos, chromNames) = _probeInfo(diskfile)
    if logger:
        logger.info('Genotypes of %d individuals at %d loci are found' % (numInds, len(lociNames)))

    pop = Population(size=numInds, ploidy=2, loci=[len(lociPos)],
        lociPos=lociPos, lociNames=lociNames, chromNames=chromNames,
        alleleNames=[alleleNames[x] for x in lociNames], subPopNames=["plink"])
    load_population(pop, diskfile, alleleNames, logger)
    logger.info("Save population to %s." % sys.argv[2])
    
    
    
    pop.save(sys.argv[2])
