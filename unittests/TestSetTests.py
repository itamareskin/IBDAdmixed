'''
Created on May 19, 2014

@author: Itamar
'''
import unittest
import os
from IBD.LDModel import LDModel
from IBD.TestSet import TestSet, GenotypePair

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.resource_path = os.path.join(os.path.split(os.path.split(__file__)[0])[0], "resources")

#     def testread_haplos(self):
#         t = TestSet()
#         t.read_haplos("")
#     
#     def testgenerate_composite_individuals(self):
#         t = TestSet()
#         t.generate_composite_individuals()
        
    def testgenerate_random_haps_inplace(self):
        m = LDModel(os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"), 
                    os.path.join(self.resource_path, "HapMap3_CEU_chr2.low.bgl.dag"),
                    max_snp_num=500)
        p = GenotypePair()
        p.generate_random_haps_inplace(m)
        self.assertEqual(p._nr_haplos, 4, "wrong number of haplotpes created")
        self.assertEqual(p._snp_num, 500, "wrong number of snps created")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()