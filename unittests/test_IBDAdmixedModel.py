'''
Created on May 25, 2014

@author: Itamar
'''
import unittest
import os
from IBD.LDModel import LDModel
from IBD.GenotypePairModel import GenotypePairModel
from IBD.TestSet import GenotypePair
from IBD.IBDAdmixedModel import ibdadmixed
import logging
import logging.handlers
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
                          max_snp_num=1000)
        #cls.m = GenotypePairModel(cls.ldm,cls.ldm,cls.ldm,cls.ldm,False,0)

    def tearDown(self):
        pass

    def testibd_admixed(self):
        p = GenotypePair()
        p.generate_random_haps_inplace(self.ldm)
        p.set_ibd_segment(0,1000)
        (pairIBD,lod_scores) = ibdadmixed(os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"),
                   [os.path.join(self.resource_path, "HapMap3_CEU_chr2.HapMap3_CEU_chr2.bgl.dag.gz")],
                   p, 
                   g=8, 
                   alphas=[1], 
                   ibd_trans=[1e-5,1],
                   win_size=250,
                   phased=False,
                   max_snp_num=self.ldm._snp_num)
        self.assertGreater(pairIBD.get_IBD_percent(self.ldm._gm), 80, "IBD was not detected")

        p = GenotypePair()
        p.generate_random_haps_inplace(self.ldm)
        (pairIBD,lod_scores) = ibdadmixed(os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"),
                   [os.path.join(self.resource_path, "HapMap3_CEU_chr2.HapMap3_CEU_chr2.bgl.dag.gz")],
                   p,
                   g=8,
                   alphas=[1],
                   ibd_trans=[1e-5,1],
                   win_size=250,
                   phased=False,
                   max_snp_num=self.ldm._snp_num)
        self.assertLess(pairIBD.get_IBD_percent(self.ldm._gm), 20, "IBD was not detected")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()