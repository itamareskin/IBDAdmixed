'''
Created on Aug 3, 2013

@author: eskin
'''


from IBD.cIBD import cPairIBD, cPopulationIBD
import pylab as P

map_file = "/home/eskin/Data/IBDAdmixed/fromLecs/HapMap3_CEU_chr2.map"

true_ibd = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/AfricanAmericans8.trueibd.dat", map_file)
#true_ibd.filter_by_length(0,100)
with open(map_file) as gm_f:
    dists = [float(x.split(" ")[2]) for x in gm_f.readlines()]
true_ibd.calc_dists(dists)

beagle  = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/beagle.ibd.txt", map_file)
# beagle.filter_by_score(-1,0)
beagle.filter_by_human_pairs(true_ibd.keys())
(power, FDR, FPR) = beagle.stats_win(true_ibd,116430,25)
print "beagle3", power, FDR, FPR

beagle  = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/AfricanAmericans.beagle.out3.ibd.txt", map_file)
# beagle.filter_by_score(600,1000000)
beagle.filter_by_human_pairs(true_ibd.keys())
(power, FDR, FPR) = beagle.stats_win(true_ibd,116430,25)
print "beagle4", power, FDR, FPR

ibd_admixed = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/AfricanAmericans8.IBDAdmixed3.dat", map_file)
ibd_admixed.filter_by_human_pairs(true_ibd.keys())
(power, FDR, FPR) = ibd_admixed.stats_win(true_ibd,116430,25)
print "IBDAdmixed", power, FDR, FPR

parente = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/parente.ibd.txt", map_file)
parente.filter_by_human_pairs(true_ibd.keys())
(power, FDR, FPR) = parente.stats_win(true_ibd,116430,25)
print "parente", power, FDR, FPR