#load dependencies
import simuPOP as sim
import math as math
import random as random
import simuOpt, types
from simuPOP.utils import migrSteppingStoneRates

###constant variables here###
numSubPops = 60 #this should stay even
subPopSize = 1000
r1 = 0.5
r2 = 0.5
mutRate = 0
hybridPref = 1
reportGen = 10
gen=1000
dir = '/home/singhal/Desktop/hybridSim/results2/'

###setting variables from command line###
options = [
    {'arg':'m:', 
     'longarg':'m=', 
     'default':0., 
     'label':'migration rate', 
     'type': 'numbers',
     'validate':simuOpt.valueBetween(0.,1.),
     'description':'Migration Rate Between Adjoining Demes (proportional)'
    },
    {'arg':'a:', 
     'longarg':'am=', 
     'default':0.,
     'label':'Assortative Mating', 
     'type': 'numbers',
     'description':'Assortative Mating -- 0 means no assortative mating',
     'validate':simuOpt.valueBetween(0.,1.)     
    },
    {'arg':'l:', 
     'longarg':'nloci=', 
     'default':5, 
     'label':'Number of Loci',
      'type': 'numbers',
     'description':'Number of Loci',
     'validate':simuOpt.valueGT(0)
    },    
    {'arg':'s:', 
     'longarg':'selection=', 
     'default':0., 
     'label':'Selection against hybrids -- 0 means no selection',
      'type': 'numbers',
     'description':'Selection against hybrids',
     'validate':simuOpt.valueBetween(0.,1.)
    },
    {'arg':'r:', 
     'longarg':'rep=', 
     'default':0, 
     'label':'Replicate number',
      'type': 'numbers',
     'description':'Replicate number if running same conditions multiple times',
     'validate':simuOpt.valueBetween(0,100)
    },
]

###
allParam = simuOpt.Params(options)
# read in the parameters
allParam.getParam(gui=False)
allOptions = allParam.asDict()
selection = allOptions['selection'][0]
am = allOptions['am'][0]
m = allOptions['m'][0]
nloci = allOptions['nloci'][0]
rep = allOptions['rep'][0]

###define some functions that I will need###
def getAllele(pop):
	
	value = "%s\t" % (pop.dvars().gen)
	f.write(str(value))

	#to print clines
	dist = range(numSubPops)	
	for a in range(nloci):
		alleleFreqs = []
		for b in range(numSubPops):
			alleleFreqs.append(pop.dvars(b).alleleFreq[a][0])
		(a,b) = cline(alleleFreqs,dist)
		width = 4/a		
		center = -b/a
		value = "%.3f\t%.3f\t" % (width,center)
		f.write(str(value))

	cp1 = int(center)
	cp2 = cp1+1
	
	#to print LD
	for i in haploList:
		(a,b) = i[0], i[1]
		##to get ld 
		pA = (pop.dvars(cp1).alleleFreq[a][1] + pop.dvars(cp2).alleleFreq[a][1])/2
		pa = (pop.dvars(cp1).alleleFreq[a][0] + pop.dvars(cp2).alleleFreq[a][0])/2
		pB = (pop.dvars(cp1).alleleFreq[b][1] + pop.dvars(cp2).alleleFreq[b][1])/2
		pb = (pop.dvars(cp1).alleleFreq[b][0] + pop.dvars(cp2).alleleFreq[b][0])/2	
		pAB = (pop.dvars(cp1).haploFreq[i][(1,1)] + pop.dvars(cp2).haploFreq[i][(1,1)])/2
		if pA*pa*pB*pb > 0:
			r = (pAB - pA * pB)/math.sqrt(pA*pa*pB*pb)
		else:
			r = 0
		value = "%.3f\t" % (r)
		f.write(str(value))

	#to print HWE measures
	for a in range(nloci):
		exp = 2 * ( (pop.dvars(cp1).alleleFreq[a][1] + pop.dvars(cp2).alleleFreq[a][1])/2 ) * ( (pop.dvars(cp1).alleleFreq[a][0] + pop.dvars(cp2).alleleFreq[a][0])/2 )
		obs = (pop.dvars(cp1).genoFreq[a][(1,0)] + pop.dvars(cp1).genoFreq[a][(0,1)] + pop.dvars(cp2).genoFreq[a][(1,0)] + pop.dvars(cp2).genoFreq[a][(0,1)])/2
		if exp > 0:
			fis = 1 - (obs/exp)
		else:
			fis = 0
		value = "%.3f\t" % (fis)
		f.write(str(value))
	
	#to print hybrid indices	
	hiList1 = list(pop.indInfo('ancestry', subPop=[cp1]))
	hiList2 = list(pop.indInfo('ancestry', subPop=[cp2]))
	hiList1.extend(hiList2)

	for i in range(50):
		value = "%.3f\t" % (random.choice(hiList1))
		f.write(str(value))
	value = "\n"
	f.write(str(value)),
			
	return True

def logfunc (a):
	b = math.log(a/(1-a))
	return b

def cline(allele, dist):
	alleleNew = []
	distNew = []
	for i in range(len(allele)):
		if allele[i] > 0 and allele[i] < 1:
			alleleNew.append(allele[i])
			distNew.append(dist[i])
	
	y = map(logfunc,alleleNew)

	X = []
	Y = []
	for i in range(len(y)):
		if y[i] > -2 and y[i] < 2:
			X.append(distNew[i])
			Y.append(y[i])
		
	if len(X) < 2:
		a = 4
		b = -(numSubPops/2 - 1) * 4 
	else:
		N = len(X)
		Sx = Sy = Sxx = Syy = Sxy = 0.0
		for x, y in map(None, X, Y):
			Sx = Sx + x
			Sy = Sy + y
			Sxx = Sxx + x*x
			Syy = Syy + y*y
			Sxy = Sxy + x*y
		det = Sxx * N - Sx * Sx
		a, b = (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det
	return a, b
  
###define the component parts of my simulation, this should not be affected by what the user sets###
#need to add some sort of tracker to the name...
name = "sim_m%s_s%s_n%s_am%s_rep%s" % (m,selection,nloci,am,rep)
simuOpt.setOptions(quiet=True,alleleType='short')

s = 1 - selection
migrRate = migrSteppingStoneRates(m,numSubPops)
subPopList = [subPopSize]*numSubPops
rRates = [r1]*1 + [r2] * (nloci - 1)

lnames = []
for i in range(nloci): 
	lnames.append("l%s" % (i)) 

haploList = []
if nloci == 2:
	haploList.append((0,1))
elif nloci == 3:
	haploList.extend([(0,1),(0,2),(1,2)])
elif nloci > 3:
	haploList.extend([(0,1),(2,3),(3,4),(4,5)])

#getting my output file ready	
file = dir + name + '.out'
f = open(file, 'w')

header = ['gen']
for i in lnames:
	header.append(str('w_'+ i))
	header.append(str('c_'+ i))
for i in haploList:
	header.append(str('LD_' + str(i)))		
for i in lnames:
	header.append(str('HWE_' + i))
for i in range(50):
	header.append(str('hybInd_' + str(i)))
h = str("\t".join(header))
f.write(h)	 	
value = "\n"
f.write(str(value))

pop = sim.Population(
	size=subPopList,
	ploidy=2,
	loci=nloci,
	lociNames=lnames,
	infoFields=['migrate_to','fitness','ancestry','index'],
	)
pop.dvars().name = name

#need to initialize genotypes; assume that all loci start at one
sim.initGenotype(pop,freq=[0,1],subPops=range(len(subPopList)/2))
#need to make it 50/50 male and female. note that this is not probabilistic but forced.
sim.initSex(pop,maleProp=0.5)
sim.initInfo(pop, [0]*(subPopSize*numSubPops/2) + [1]*(subPopSize*numSubPops/2),infoFields='ancestry')
#define mating scheme
if hybridPref == 1 and am > 0:
	pop.setVirtualSplitter(sim.InfoSplitter(field='ancestry', cutoff = [0.05,0.95]))
	vsp = [0,1,2]
	vspList = []
	for i in range(numSubPops):
		for j in vsp: 
			vspList.append((i,j))
elif hybridPref == 0 and am > 0:
	pop.setVirtualSplitter(sim.InfoSplitter(field='ancestry', cutoff = [0.5]))
	vsp = [0,1]
	vspList = []
	for i in range(numSubPops):
		for j in vsp: 
			vspList.append((i,j))
	
"""
genotype transmitters; all loci are typical mendelian with recombination rates as set 
by sim.Recombinator
ancestry is a 100% accurate hybrid index value, where an individual's value is determined 
by the mean between its parents hybrid index
"""
transmitters=[
   sim.MendelianGenoTransmitter(),
   sim.InheritTagger(mode=sim.MEAN, infoFields='ancestry'),
   sim.Recombinator(rates=rRates,loci=lnames)]
		  
"""
defines all the operations that this population will undergo
1. migration
2. mutation
3. check sizes
4. kill any VSPs with too few individuals
5. selection
6. do mating
7. measure stats
"""
preOpList = [
			sim.Migrator(rate=migrRate),
			sim.SNPMutator(u=mutRate),	
			]
if am > 0:
	preOpList.append(sim.InitSex(maleProp=0.5,subPops=vspList))
if nloci < 3: 
	selOp=sim.MapSelector(loci=1,fitness={(0,0):1,(0,1):s,(1,1):1})
	preOpList.append(selOp)
else:
	#need a more complicated selection operator based on individual genotypes, heterozygosity
	locusSel = round(math.exp(math.log(s)/nloci),4)
	selList=[]
	for i in range(nloci)[2:]:
		op = sim.MapSelector(loci=i,fitness={(0,0):1,(0,1):locusSel,(1,1):1})
		selList.append(op)
	selOp=sim.MlSelector(selList,mode=sim.MULTIPLICATIVE)
	preOpList.append(selOp)

"""	
define mating scheme
"""
matingSchemeList=[]
#complete random mating
if am == 0:	
	matingSchemeList.append(sim.RandomMating(ops=transmitters))	
#some incomplete mating
else:
	matingSchemeList.append(sim.RandomMating(ops=transmitters))
	for i in range(numSubPops):
		for j in vsp:
			op = sim.RandomMating(subPops=[(i,j)],weight=-am,ops=transmitters)
			matingSchemeList.append(op)

"""
time to start the actual evolution
"""
simu = sim.Simulator(pop)			 
simu.evolve (		
	preOps=preOpList,
		
	matingScheme=sim.HeteroMating(
       matingSchemes=matingSchemeList,
   ),
	
	postOps = [	
		sim.Stat(genoFreq=lnames,alleleFreq=lnames,HWE=lnames, haploFreq=haploList, vars=['haploFreq_sp','HWE_sp','alleleFreq_sp','genoFreq_sp'],step=reportGen),	
		sim.PyOperator(func=getAllele,step=reportGen),				
	],		
	
	gen=gen
	)
f.close()	
