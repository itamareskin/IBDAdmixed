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

    def test_read_haplos(self):
        t = TestSet()
        t.read_haplos(os.path.join(self.resource_path, "ceu.genos.dat"), max_snp_num=500)
        t.write_haplos(os.path.join(self.resource_path, "ceu.genos.out.dat"))
        x = 1
#     
#     def testgenerate_composite_individuals(self):
#         t = TestSet()
#         t.generate_composite_individuals()
        
    def test_generate_random_haps_inplace(self):
        m = LDModel(os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"), 
                    os.path.join(self.resource_path, "HapMap3_CEU_chr2.low.bgl.dag"),
                    max_snp_num=500)
        p = GenotypePair()
        p.generate_random_haps_inplace(m)
        self.assertEqual(p._nr_haplos, 4, "wrong number of haplotpes created")
        self.assertEqual(p._snp_num, 500, "wrong number of snps created")

    def test_write_haplos(self):
        m = LDModel(os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"),
                    os.path.join(self.resource_path, "HapMap3_CEU_chr2.low.bgl.dag"),
                    max_snp_num=500)
        p = GenotypePair()
        p.generate_random_haps_inplace(m)
        p.write_haplos(os.path.join(self.resource_path, "test.write.haplos.dat"))

    def test_chr1(self):
        m = LDModel(os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"),
                    os.path.join(self.resource_path, "HapMap3_CEU_chr2.low.bgl.dag"),
                    max_snp_num=500)
        p = GenotypePair()
        p.generate_random_haps_inplace(m)
        p.write_haplos(os.path.join(self.resource_path, "test.write.haplos.dat"))

    def test_get_slice(self):
        t = TestSet()
        t.read_haplos(os.path.join(self.resource_path, "ceu.genos.dat"), max_snp_num=500)
        sliced = t.get_slice(0,50)
        self.assertEqual(sliced._snp_num, 50, "wrong number of snps created")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()