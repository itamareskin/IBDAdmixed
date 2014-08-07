'''
Created on Apr 1, 2013

@author: eskin3
'''
import scipy.stats as sc
#import numpy as nmp
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
import string
import math
import random
#from tables import *
#import tables
from IBD.IBDSegments import PairIBD,PopIBD
#from bx.intervals import IntervalTree, Interval
import time

ibdadmixed_f = open("ibdadmixed.ibd.txt", 'r')
ibdAdmixed = PopIBD.from_string(ibdadmixed_f.readlines())
ibdadmixed_f.close()
beagle_f = open("beagle4.ibd.txt", 'r')
beagle = PopIBD.from_string(beagle_f.readlines())
beagle_f.close()
true_f = open("AfricanAmericans4.trueibd.dat", 'r')
true_ibd = PopIBD.from_string(true_f.readlines())
true_f.close()

beagle_power = beagle.calc_power_intervals(true_ibd)
ibdadmixed_power = ibdAdmixed.calc_power_intervals(true_ibd)
x = np.arange(5)
w = 0.35
fig = plt.figure()
ax = fig.add_subplot(111)
rects1 = ax.bar(x, beagle_power, color='b', width=w)
rects2 = ax.bar(x+w, ibdadmixed_power, color='g', width=w)
ax.set_xlabel('Segment Length')
ax.set_ylabel('Power')
ax.legend( (rects1[0], rects2[0]), ('BEAGLE', 'IBD Admixed') )
ax.set_xticks(x+w)
ax.set_xticklabels( ('1cM', '2cM', '3cM', '4cM', '5cM') )
ax.grid(True)
def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.02*height, '%.2f'%(height*100),
                ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)
plt.show()