'''
Created on May 25, 2014

@author: Itamar
'''
import unittest
import os
import logging
import logging.handlers

from IBD.LDModel import LDModel
from IBD.TestSet import GenotypePair
from IBD.NaiveModel import naivemodel

posterior_probs_logger = logging.getLogger('posteriorprobs')
posterior_probs_logger.addHandler(logging.FileHandler("posteriorprobs.log"))
posterior_probs_logger.setLevel(logging.DEBUG)
inner_forward_probs_logger = logging.getLogger('innerforwardprobs')
inner_forward_probs_logger.addHandler(logging.FileHandler("innerforwardprobs.log"))
inner_forward_probs_logger.setLevel(logging.ERROR)

class Test(unittest.TestCase):


    @classmethod
    def setUpClass(cls):
        cls.resource_path = os.path.join(os.path.split(os.path.split(__file__)[0])[0], "resources")
        cls.ldm = LDModel(os.path.join(cls.resource_path, "HapMap3_CEU_chr2.map"), 
                          os.path.join(cls.resource_path, "HapMap3_CEU_chr2.HapMap3_CEU_chr2.bgl.dag.gz"),
                          max_snp_num=100)
        #cls.m = GenotypePairModel(cls.ldm,cls.ldm,cls.ldm,cls.ldm,False,0)

    def tearDown(self):
        pass

    def testnaivemodel(self):
        p = GenotypePair()
        p.generate_random_haps_inplace(self.ldm)
        p.set_ibd_segment(0,100)
        
        (pairIBD, lod_scores) = naivemodel(os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"),
                   [os.path.join(self.resource_path, "HapMap3_CEU_chr2.HapMap3_CEU_chr2.bgl.dag.gz")],
                    p, max_snp_num=100, win_size=25)

        x=1

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()