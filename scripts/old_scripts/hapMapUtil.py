
'''
This python module provides several utility functions that handles HapMap
populations. When used as a script, this module creates a population using
selected markers and populations.
'''

import simuOpt as opt
opt.setOptions(alleleType='lineage', quiet=True, optimized=True)
from simuPOP import *
from types import *

import os
import sys
import random
import string

def mergeHapMapPops(HapMap_dir, HapMap_pops, chrom, logger=None):
    '''
    Load HapMap dataset for multiple populations and merge them.
    The important step is to find a common set of markers and make sure
    alleles are recoded if necessary. (Alleles are sometimes coded differently
    in different HapMap populations.)

    HapMap_dir
        Directory where HapMap files (in simuPOP format) are saved.

    HapMap_pops
        HapMap populations to load.

    chrom
        Which chromosome to load and combine.

    logger
        A logger to record what is going on.
    '''
    pop = None
    for HapMap_pop in HapMap_pops:
        filename = os.path.join(os.path.expanduser(HapMap_dir), 'HapMap3_%s_chr%d.pop' % \
            (HapMap_pop, chrom))
        if logger is not None:
            logger.info('Loading HapMap Population %s' % filename)
        pop1 = loadPopulation(filename)
        if pop is None:
            pop = pop1
            continue
        # need to be merged.
        markers1 = set(pop.lociNames())
        markers2 = set(pop1.lociNames())
        common_markers = markers1 & markers2
        remove1 = markers1 - common_markers
        remove2 = markers2 - common_markers
        if len(remove1) > 0:
            if logger is not None:
                logger.info('Removing %d markers (%.2f percent) from population %s' % \
                    (len(remove1), len(remove1)*100./pop.totNumLoci(), pop.subPopNames()))
            pop.removeLoci([pop.locusByName(x) for x in remove1])
        if len(remove2) > 0:
            if logger is not None:
                logger.info('Removing %d markers (%.2f percent) from population %s' % \
                    (len(remove2), len(remove2)*100./pop1.totNumLoci(), pop1.subPopNames()))
            pop1.removeLoci([pop1.locusByName(x) for x in remove2])
        pop.addIndFrom(pop1)
    return pop



def getHapMapMarkers(HapMap_dir, names = [], chroms=[], HapMap_pops=['CEU'],
        startPos = [], endPos = [], numMarkers = [], minAF = 0, minDist = 0,
        logger=None):
    '''
    Return a population with specified HapMap markers.

    HapMap_dir
        Directory where the HapMap data has been saved, using script
        loadHapMap_r22.py from the simuPOP cookbook.

    names
        (Optional) A list of marker names. If given, only markers in this
        list will be selected.

    chroms
        A list of chromosomes to look in. If empty, all 22 autosomes
        will be tried. Chromosome index starts from 1. (1, ..., 22).

    HapMap_pops
        HapMap populations to load, can be one or both of 'CEU' and 'YRI'.

    startPos, endPos, numMarkers
        These list should be empty or match the length of ``chroms``.
        They specify the starting, ending position on each chromosome
        (in basepair), and number of markers to load.

    minAF
        Minimal minor allele frequency

    minDist
        Minimal distance between adjacent markers

    logger
        A logger to record what is going on.
    '''
    
    def paramExpandList(param, size, err=''):
        '''If parameter param is
        - a number: return a list of specified size.
        - a list of size 1: expand it to size and return.
        - a list of more than one element: raise an error if size mismatch.
        - an empty list: return param.
        '''
        if type(param) in [IntType, LongType, FloatType]:
            return [param]*size
        elif type(param) in [TupleType, ListType]:
            if len(param) == 1:
                return list(param)*size
            elif len(param) == 0:
                return param
            elif len(param) != size:
                raise exceptions.ValueError(err)
        return param
    
    if len(chroms) == 0:
        chs = range(1, 23)
    else:
        chs = chroms
    # read in HapMap data file
    pop = None
    geneticMap = {}
    #recombRates = {}
    sPos = paramExpandList(startPos, len(chs), 'Incorrect starting position')
    ePos = paramExpandList(endPos, len(chs), 'Incorrect ending position')
    nMarkers = paramExpandList(numMarkers, len(chs), 'Incorrect number of markers')
    #
    print sPos, ePos, chs, names
    for chIdx, ch in enumerate(chs):
        markers = []
        chPop = mergeHapMapPops(HapMap_dir, HapMap_pops, ch, logger)
        # Trim markers by marker names
        if len(names) != 0:
            if logger is not None:
                logger.info("Select markers using a list of %d markers from chromosome %s..." % (len(names), ch))
            # the markers may not be in order...
            indexes = []
            for name in names:
                try:
                    idx = chPop.locusByName(name)
                    if not idx in indexes:
                        indexes.append(idx)
                except:
                    pass
            chPop.removeLoci(keep = indexes)
        #
        if chPop.totNumLoci() == 0:
            continue
        if minAF > 0:
            stat(chPop, alleleFreq=range(chPop.totNumLoci()))
        # Trim by start, end position ...
        indexes = []
        lastPos = 0
        for loc in range(chPop.totNumLoci()):
            pos = chPop.locusPos(loc)
            if len(sPos) > 0 and pos < sPos[chIdx]:
                continue
            if len(ePos) > 0 and pos > ePos[chIdx]:
                continue
            if lastPos > 0 and pos - lastPos < minDist:
                continue
            if minAF > 0:
                maf = chPop.dvars().alleleFreq[loc][0]
                maf = min(maf, 1 - maf)
                if maf < minAF:
                    continue
            if len(nMarkers) > 0 and len(indexes) >= nMarkers[chIdx]:
                break
            indexes.append(loc)
            lastPos = pos
        if len(indexes) > 0:
            if logger is not None:
                logger.info('%s markers are found on chromosome %d ' % (len(indexes), ch))
            chPop.removeLoci(keep=indexes)
            geneticMap.update(chPop.dvars().geneticMap)
            #recombRates.update(chPop.dvars().recombRates)
            chPop.vars().clear()
            if pop is None:
                pop = chPop
            else:
                pop.addChromFrom(chPop)
        else:
            if logger is not None:
                logger.info('No qualified marker is found on chromosome %d ' % ch)
            del chPop
    pop.dvars().geneticMap = geneticMap
    #pop.dvars().recombRates = recombRates
    return pop
    
    # awk '{print substr($1,1,length($1)/2); print substr($1,length($1)/2+1,length($1)/2)}' test4.dat > test4.2line.dat