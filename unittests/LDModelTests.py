'''
Created on May 20, 2014

@author: Itamar
'''
import unittest
import pkgutil
import os
from IBD.LDModel import LDModel

class Test(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.resource_path = os.path.join(os.path.split(os.path.split(__file__)[0])[0], "resources")
        cls.m = LDModel(os.path.join(cls.resource_path, "HapMap3_CEU_chr2.map"), 
                        os.path.join(cls.resource_path, "HapMap3_CEU_chr2.low.bgl.dag"),
                        max_snp_num=500)

    def testread_from_bgl_file(self):
        self.m.read_from_bgl_file(os.path.join(self.resource_path, "HapMap3_CEU_chr2.low.bgl.dag"))
        self.assertEqual(self.m._snp_num, 500, "Not all snps were read correctly from map file") 
        
    def testget_slice_model(self):
        self.m.read_from_bgl_file(os.path.join(self.resource_path, "HapMap3_CEU_chr2.low.bgl.dag"))
        sliced = self.m.get_slice_model(0, 50)
        self.assertEqual(sliced._snp_num, 50, "wrong number of snps in slice")
    
#     def testread_from_bgl_file(self):
#         self.m.read_from_bgl_file("K:\\Projects\\ibdadmixed\\tests\\HapMap3_CEU_chr2.low.bgl.dag")
#         self.assertEqual(self.m._snp_num, 500, "Not all snps were read correctly from map file")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()