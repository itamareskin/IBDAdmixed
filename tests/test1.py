'''
Created on Aug 3, 2013

@author: eskin
'''


from IBD.cIBD import cPairIBD, cPopulationIBD
import pylab as P
import numpy as np
import math

map_file = "/home/eskin/Data/IBDAdmixed/fromLecs/HapMap3_CEU_chr2.map"

true_ibd = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/New/single.ceu.randb.trueibd.txt", map_file)
#true_ibd = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/artificial.admixed.test.trueibd.txt", map_file)
#true_ibd.filter_by_length(1.5,10)

#true_ibd.filter_by_length(5,100)
#with open(map_file) as gm_f:
#    dists = [float(x.split(" ")[2]) for x in gm_f.readlines()]
#true_ibd.calc_dists(dists)
#P.figure()
#P.hist([x[3] for x in true_ibd.to_list()])
#P.show()

# germline = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/artificial.admixed.test.germline.ibd3.txt", map_file)
# germline.filter_by_length(1,100)
# germline.merge_all(50)
# (power, FDR, FPR) = germline.stats_win(true_ibd,116430,25)
# print power, FDR, FPR

gm_f = open(map_file)
data = gm_f.readlines()
dists = [float(x.split(" ")[2]) for x in data]

ibd_admixed = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/New/single.ceu.ceurandb.switch10.ibdadmixed.txt", map_file)
#true_ibd.filter_by_human_pairs(ibd_admixed.keys())
ibd_admixed.filter_by_human_pairs(true_ibd.keys())
true_ibd.filter_by_human_pairs(ibd_admixed.keys())

ibd_admixed2 = cPopulationIBD.from_string(ibd_admixed.to_string(),dists)
scores = [x[2] for x in ibd_admixed.to_list()]
print min(scores),max(scores)
#ibd_admixed2.filter_by_score(-55,1e10)
ibd_admixed2.merge_all_fast(overlap = 1, max_val = True, merge_diff_vals=True)
ibd_admixed2.calc_dists(dists)
ibd_admixed2.filter_by_length(0.5,1000)
(power, FDR, FPR) = ibd_admixed2.stats_win(true_ibd,116430,25)
print power, FDR, FPR

f = open("/home/eskin/Data/IBDAdmixed/fromLecs/New/ceu.yri.tsilwkref.ibdadmixed.filt.txt","w")
f.write(ibd_admixed2.to_string())
f.close()
 
scores = [int(x[2]) for x in ibd_admixed2.to_list()]
print min(scores),max(scores)
#for score in range(0,205,5):
for score in range(-30,max(scores)+1,10):
    ibd_admixed2 = cPopulationIBD.from_string(ibd_admixed.to_string(),dists)
    ibd_admixed2.filter_by_score(score,205)
    ibd_admixed2.merge_all_fast(overlap = 1, max_val = True, merge_diff_vals=True)
    ibd_admixed2.calc_dists(dists)
    ibd_admixed2.filter_by_length(0.5,1000)
    (power, FDR, FPR) = ibd_admixed2.stats_win(true_ibd,116430,25)
    print "IBDAdmixedUnphased", score, power, FPR

parente = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/single.test.genos.parente.ibd.txt", map_file)
(power, FDR, FPR) = parente.stats_win(true_ibd,116430,25)
print "PARENTE", power, FPR
  
scores = [x[2] for x in parente.to_list()]
print min(scores),max(scores)
#for score in range(0,205,5):
for score in range(int(min(scores)),int(max(scores))+1):
#for score in range(20,38,1):
    parente2 = cPopulationIBD.from_string(parente.to_string(),dists)
    #parente2.filter_by_human_pairs(true_ibd.keys())
    parente2.filter_by_score(score,1000)
    (power, FDR, FPR) = parente2.stats_win(true_ibd,116430,25)
    print "PARENTE", score, power, FPR

# 
# x=1
#beagle  = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/artificial.admixed.test.beagle.ibd.txt", map_file)
#print beagle.stats_win(true_ibd,116430,25)
#scores = [x[2] for x in beagle.to_list()]
for score in range(0,40):
    beagle = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/artificial.admixed.test.genos.parente.ibd.txt", map_file)
    #beagle.filter_by_length(2,1000)
    beagle.filter_by_score(0,math.pow(10, -score))
    (power, FDR, FPR) = beagle.stats_win(true_ibd,116430,25)
    print "fastIBD", math.pow(10, -score), power, FPR


#true_ibd.filter_by_human_pairs(ibd_admixed.keys())
#ibd_admixed.filter_by_human_pairs(ibd_admixed.keys())
(power, FDR, FPR) = ibd_admixed.stats_win(true_ibd,116430,25)
print power, FDR, FPR

beagle  = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/artificial.admixed.test.beagle.ibd.txt", map_file)
beagle.filter_by_score(0,1e-11)
print beagle.stats_win(true_ibd,116430,25)
beagle4  = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/aritficial2.beagle4.out7.ibd.txt", map_file)
beagle4.filter_by_score(0,1000000)
print beagle4.stats_win(true_ibd,116430,25)

#beagle_old = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/beagle.ibd.txt")

#beagle_new = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/beagle4.ibd5.txt")
#beagle_new2 = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/beagle4.ibd6.txt")
#beagle_new3 = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/beagle4.ibd3.txt")
#beagle_new4 = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/beagle4.ibd4.txt")

#ibd_admixed = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/Europeans.IBDAdmixed3.dat", map_file)
    
for score in range(0,20,1):
    parente = cPopulationIBD.fast_deserialize("/home/eskin/Data/IBDAdmixed/fromLecs/aritficial3.comp.parente.ibd.txt", map_file)
    parente.filter_by_human_pairs(true_ibd.keys())
    parente.filter_by_score(score)
    (power, FDR, FPR) = parente.stats_win(true_ibd,116430,25)
    print score, power, FPR
#parente.merge_all()
#f = open("/home/eskin/Data/IBDAdmixed/fromLecs/parente.artifical2.ibd.txt","w")
#f.write(parente.to_string())
#f.close()

#true_ibd.filter_by_human_pairs(ibd_admixed.keys())
#beagle.filter_by_human_pairs(ibd_admixed.keys())
# parente.filter_by_human_pairs(ibd_admixed.keys())
# beagle_old.filter_by_human_pairs(ibd_admixed.keys())
# beagle_new.filter_by_human_pairs(ibd_admixed.keys())
# beagle_new2.filter_by_human_pairs(ibd_admixed.keys())
# beagle_new3.filter_by_human_pairs(ibd_admixed.keys())
# beagle_new4.filter_by_human_pairs(ibd_admixed.keys())

#print beagle.stats_win(true_ibd,100000,25)
# print beagle_old.stats_win(true_ibd,116415,25)
#print ibd_admixed.stats_win(true_ibd,100000,25)
print parente.stats_win(true_ibd,116430,25)
# print beagle_new.stats_win(true_ibd,116415,25)
# print beagle_new2.stats_win(true_ibd,116415,25)
# print beagle_new3.stats_win(true_ibd,116415,25)
# print beagle_new3.stats_win(true_ibd,116415,25)