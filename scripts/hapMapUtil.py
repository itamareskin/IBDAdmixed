#!/usr/bin/env python
#
# Purpose:
#     This python module provides several utility functions that handles HapMap
#     populations.
#
# License:
#     This program is freely available in the simuPOP's online cookbook
#     (http://simupop.sourceforge.net/cookbook). You can redistribute it and/or
#     modify it freely, with or without this license notice. However, this
#     license notice should be present in the online version of this file.
#
#     This program is NOT part of simuPOP and is NOT protected by simuPOP's GPL
#     license. It is provided in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY. If you notice any bug or have some new ideas, you can
#     modify this file and, as a courtesy to other simuPOP users, incoporate
#     your changes to the online version of this file. If you are uncertain
#     about your changes, please feel free to discuss your changes in the
#     simuPOP mailinglist (simupop-list@lists.sourceforge.net, subscription
#     required).
#
# Change Log:
#     2009-06-25 Bo Peng <bpeng@mdanderson.org>
#
#         Add functions mergeHapMapPops and getHapMapMarkers
# 

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
        filename = os.path.join(HapMap_dir, 'HapMap_%s_chr%d.pop' % \
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
        #for i in range(pop.totNumLoci()):
        #    assert pop.alleleNames(i) == pop1.alleleNames(i)
        # merge two populations
        # check coding.
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
    genDist = {}
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
            genDist.update(chPop.dvars().genDist)
            chPop.vars().clear()
            if pop is None:
                pop = chPop
            else:
                pop.addChromFrom(chPop)
        else:
            if logger is not None:
                logger.info('No qualified marker is found on chromosome %d ' % ch)
            del chPop
    pop.dvars().genDist = genDist
    return pop

HapMap_pops = ['CEU', 'YRI', 'JPT+CHB']

options = [
    {
    'longarg': 'HapMap_dir=',
    'default': '../HapMap',
    'useDefault': True,
    'label': 'HapMap data directory',
    'description': '''Directory to store HapMap data in simuPOP format. The 
        HapMap data is expected to be downloaded and saved in simuPOP format
        using script loadHapMap_r22.py from the simuPOP online cookbook. The
        files have names such as HapMap_CEU_chr10.pop.''',
     'allowedTypes': [StringType],
     'validate': opt.valueValidDir(),
    },
    {
    'longarg': 'HapMap_pops=',
    'default' : HapMap_pops,
    'useDefault': False,
    'label' : 'Name of populations',
    'description': '''Which HapMap populations to use?''',
    'allowedTypes': [ListType, TupleType],
    'chooseFrom': HapMap_pops,
    'validate': opt.valueListOf(HapMap_pops),
    },
    {
    'longarg': 'markerList=',
    'default': '',
    'useDefault': True,
    'label': 'Marker list file',
    'description': '''A file with a list of marker names, in the form of
        "marker_name chrom_number marker_pos". The first line is ignored if it
        is a header line. Lines that start with '#' are ignored. An additional
        parameter 'markerListCols' can be used to specify the columns if the
        fields are not in order.''',
    'allowedTypes': [StringType],
    'validate': opt.valueOr(opt.valueEqual(''), opt.valueValidFile()),
    },
    {
    'longarg': 'markerListCols=',
    'default': [0, 1, 2],
    'useDefault': True,
    'label': 'Columns in markerList file',
    'description': '''Columns for marker name, chromosome number and position
        in the marker list file (start at 0). It should be [1, 9, 10] for an
        illumina annotation fiel.''',
    'allowedTypes': [TupleType, ListType],
    'validate': opt.valueListOf(opt.valueGE(0)),
    },
    {
    'longarg': 'chroms=',
    'default': [],
    'useDefault': True,
    'label': 'Chromosomes to use',
    'description': 'A list of chromosomes (1-22) to use.',
    'allowedTypes': [TupleType, ListType],
    'validate': opt.valueListOf(opt.valueBetween(1, 22)),
    },
    {
    'longarg': 'numMarkers=',
    'default': [],
    'useDefault': True,
    'label': 'Number of markers to use',
    'description': '''Number of markers to use for each chromosome. This
        parameter should be ignored if it is unspecified or is set to zero
        for some chromosomes.
        ''',
    'allowedTypes': [TupleType, ListType],
    'validate': opt.valueOr(opt.valueGT(0), opt.valueListOf(opt.valueGE(0)))
    },
    {
    'longarg': 'startPos=',
    'default': [],
    'useDefault': True,
    'label': 'staring position (bp)',
    'description': '''Starting position of the markers on each chromosome.
        The beginning of the chromosomes will be assumed if this parameter
        is unspecified or is set to zero.''',
    'allowedTypes': [TupleType, ListType],
    'validate': opt.valueOr(opt.valueGE(0), opt.valueListOf(opt.valueGE(0)))
    },
    {
    'longarg': 'endPos=',
    'default': [],
    'useDefault': True,
    'label': 'Ending position (bp)',
    'description': '''Ending position of the markers on each chromosome.
        The end of the chromosomes will be assumed if this parameter is
        unspecifed or is set to zero.''',
     'allowedTypes': [TupleType, ListType],
     'validate': opt.valueOr(opt.valueGE(0), opt.valueListOf(opt.valueGE(0)))
    },
    {
    'longarg': 'minAF=',
    'default': 0,
    'useDefault': True,
    'label': 'Minimal minor allele frequency',
    'description': '''Minimal allele frequency of selected markers.''',
    'allowedTypes': [IntType, LongType, FloatType],
    'validate': opt.valueBetween(0, 0.5)
    },
    {
    'longarg': 'minDist=',
    'default': 0,
    'useDefault': True,
    'label': 'Minimal distance between markers',
    'description': '''Minimal distance between adjacent markers''',
    'allowedTypes': [IntType, LongType, FloatType],
    'validate': opt.valueGE(0),
    },
    {
    'longarg': 'filename=',
    'label': 'Filename to save population',
    'default': '',
    'useDefault': True,
    'description': '''Name of the population or an absolute path to
        a file. This parameter will be ignored if an empty string or None
        is given.''',
    'allowedTypes': [StringType],
    }
]

def demo(pop):
    return [int(x*random.randint(130, 170)/100.0) for x in pop.subPopSizes()]

if __name__ == '__main__':
    import logging
    import simuPOP as sim
    from simuPOP.utils import export
    from IBD import cIBD
    from itertools import product
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger()
    pars = opt.Params(options,
        'This script chooses specified markers from one or more HapMap\n'
        'populations and saves them in simuPOP format.\n',
        __doc__)
    if not pars.getParam():
        sys.exit(1)
    if logger is not None:
        logger.info('Reading marker list %s' % pars.markerList)
    names = []
    if pars.markerList != '':
        mlist = open(pars.markerList)
        for line in mlist.readlines():
            if line.startswith('#') or line.strip() == '':
                continue
            fields = line.split(',')
            if len(fields) == 1:
                fields = line.split()
            try:
                name = fields[pars.markerListCols[0]]
                ch = int(fields[pars.markerListCols[1]].lstrip('chr'))
                pos = float(fields[pars.markerListCols[2]])
            except Exception, e:
                if logger is not None:
                    logger.debug("Ignoring line %s..." % line[:50])
                    logger.debug(e)
                continue
            if len(pars.chroms) > 0 and ch not in pars.chroms:
                continue
            print 'Adding locus %s (#%d)' % (name, len(names))
            names.append(name)
        if len(names) == 0:
            print 'Failed to locate any locus. Please check your parameters.'
            sys.exit(0)
    #
    pop = getHapMapMarkers(pars.HapMap_dir, 
        names = names,
        chroms=pars.chroms, 
        HapMap_pops=pars.HapMap_pops,
        startPos = pars.startPos,
        endPos = pars.endPos,
        numMarkers = pars.numMarkers,
        minAF = pars.minAF,
        minDist = pars.minDist,
        logger=logger)
    if pars.filename != '' and logger is not None:
        logger.info('Save population to file %s' % pars.filename)
        pop.save(pars.filename)
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
        initOps=[sim.InitSex(),
                 sim.InitLineage(range(pop.popSize()), mode=sim.PER_INDIVIDUAL),
                 ],
        matingScheme=sim.RandomMating(ops=transmitters,subPopSize=demo),
        postOps=[sim.Stat(popSize=True),
                 sim.PyEval(r'"%s\n" % subPopSize')],
        gen=10,
    )
    pop.evolve(
        preOps=[sim.Migrator(rate=[[0., 0], [0.05, 0]]),
                sim.Stat(popSize=True)], 
        matingScheme=sim.HeteroMating(
            matingSchemes=[
                #sim.RandomMating(sexMode=sim.RANDOM_SEX),
                sim.RandomMating(ops=transmitters),
                sim.RandomMating(subPopSize=list(pop.subPopSizes()), subPops=[(0,0)], weight=-0.80, ops=transmitters),
                sim.RandomMating(subPopSize=list(pop.subPopSizes()), subPops=[(0,1)], weight=-0.80, ops=transmitters)
            ],
        ),
        gen=10,
    )
    
    pop.addInfoFields(['ind_id'])
    sim.initInfo(pop,values=range(pop.popSize()),infoFields='ind_id')
    
    from IBD.FoundersContainer import FoundersContainer
    from IBD.cIBD import cPairIBD, cPopulationIBD 
    
    #pop.addInfoFields(['contributors','contributor_container'])
    # initialize ancestry
    #sim.initInfo(pop, [set(ind.lineage()) for ind in pop.individuals()],
    #infoFields='contributors')
    #sim.initInfo(pop, [FoundersContainer.from_founders_list(list(ind.lineage()),pop.locusPos(0),pop.locusPos(pop.totNumLoci)) for ind in pop.individuals()],
    #infoFields='contributor_container')
    
    pop.dvars().cont_to_inds = {}
    pop.dvars().ind_to_cont_list = {}
    pop.dvars().ind_to_cont_struct = {}
    for h1 in pop.individuals():
        pop.dvars().ind_to_cont_list[h1.info('ind_id')] = set(h1.lineage())
        pop.dvars().ind_to_cont_struct[h1.info('ind_id')] = \
        [FoundersContainer.from_founders_list(list(h1.lineage(ploidy=0)),0,pop.totNumLoci()), 
         FoundersContainer.from_founders_list(list(h1.lineage(ploidy=1)),0,pop.totNumLoci())] 
        for cont in pop.dvars().ind_to_cont_list[h1.info('ind_id')]:
            if not pop.dvars().cont_to_inds.has_key(cont):
                pop.dvars().cont_to_inds[cont] = []
            pop.dvars().cont_to_inds[cont].append(h1.info('ind_id'))
            
    num_ibd_segments = 0
    total_ibd_length = 0
    count = 0
    ibd = cPopulationIBD()
    for h1 in pop.individuals():
        count+=1
            
        humans = []
        for cont in pop.dvars().ind_to_cont_list[h1.info('ind_id')]:
            for ind in pop.dvars().cont_to_inds[cont]:
                humans.append(pop.individual(ind))
                
        if logger is not None:
            if len(humans) > 0:
                logger.info("Calculating IBD for human with ID: %s (%d out of %d) with %d other humans ", h1.info('ind_id'), count, pop.popSize(), len(humans))
            
        for h2 in humans:
            if h1.info('ind_id') >= h2.info('ind_id'):
                continue

            pairIBD = cPairIBD()
            for c1,c2 in product([0,1],[0,1]):
                pairIBD.add_ibd_from_containers(pop.dvars().ind_to_cont_struct[h1.info('ind_id')][c1],pop.dvars().ind_to_cont_struct[h2.info('ind_id')][c2],0)
            if pairIBD != None:
                if len(pairIBD.to_list()) > 0:
                    ibd.add_human_pair((h1.info('ind_id'),h2.info('ind_id')), pairIBD)
                    num_ibd_segments += len(pairIBD.to_list())
                    total_ibd_length += pairIBD.get_total_segments_length()
    
    x = 1
    export(pop, format='csv', output='test4.dat', delimiter='', header=False, infoFields=[],sexFormatter=None, affectionFormatter=None)
    # awk '{print substr($1,1,length($1)/2); print substr($1,length($1)/2+1,length($1)/2)}' test4.dat > test4.2line.dat