'''
Created on May 20, 2014

@author: Itamar
'''
import unittest
import os
from IBD.GeneticMap import GeneticMap

class Test(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.resource_path = os.path.join(os.path.split(os.path.split(__file__)[0])[0], "resources")
        
    def testread_map(self):
        gm = GeneticMap(os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"), max_snp_num=500)
        self.assertEqual(gm._snp_num, 500, "Not all snps were read correctly from map file")
        
    def testget_slice(self):
        gm = GeneticMap(os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"), max_snp_num=500)
        sliced = gm.get_slice(0,250)
        self.assertEqual(sliced._snp_num, 250, "Not all snps were read correctly from map file")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testread_map']
    unittest.main()