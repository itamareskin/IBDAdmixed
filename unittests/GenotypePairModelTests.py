'''
Created on May 20, 2014

@author: Itamar
'''
import unittest
import os
from IBD.LDModel import LDModel
from IBD.GenotypePairModel import GenotypePairModel
from IBD.TestSet import GenotypePair

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.resource_path = os.path.join(os.path.split(os.path.split(__file__)[0])[0], "resources")
        cls.ldm = LDModel(os.path.join(cls.resource_path, "HapMap3_CEU_chr2.map"), 
                          os.path.join(cls.resource_path, "HapMap3_CEU_chr2.low.bgl.dag"),
                          max_snp_num=500)
        cls.m = GenotypePairModel(cls.ldm,cls.ldm,cls.ldm,cls.ldm,False,0)

    @classmethod
    def tearDownClass(cls):
        pass

    def testName(self):
        p = GenotypePair()
        p.generate_random_haps_inplace(self.ldm)
        self.m.calc_emission_probs(p)
        self.m.calc_forward_probs()
        likelihood = self.m.calc_likelihood()
        x = 1


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()