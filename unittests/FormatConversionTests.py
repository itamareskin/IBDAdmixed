__author__ = 'Itamar'

import unittest
import os
from utils.FormatConversions import convert_ped_to_bgl,convert_ped_to_lamp

class MyTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.resource_path = os.path.join(os.path.split(os.path.split(__file__)[0])[0], "resources")

    def test_ped_tp_bgl(self):
        convert_ped_to_bgl(os.path.join(self.resource_path, "HapMap3_CEU_chr2.ped"),
                           os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"),
                           os.path.join(self.resource_path, "HapMap3_CEU_chr2.bgl"),
                           os.path.join(self.resource_path, "HapMap3_CEU_chr2.markers"))

    def test_ped_tp_bgl(self):
        convert_ped_to_lamp(os.path.join(self.resource_path, "HapMap3_CEU_chr2.ped"),
                           os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"),
                           os.path.join(self.resource_path, "HapMap3_CEU_chr2.lamp.unphased.dat"),0)

        convert_ped_to_lamp(os.path.join(self.resource_path, "HapMap3_CEU_chr2.ped"),
                           os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"),
                           os.path.join(self.resource_path, "HapMap3_CEU_chr2.lamp.phased.dat"),1)
        x = 1

if __name__ == '__main__':
    unittest.main()
