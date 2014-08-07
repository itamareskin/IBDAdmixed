__author__ = 'Itamar'

import unittest
import os
from utils.FormatConversions import convert_ped_to_bgl

class MyTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.resource_path = os.path.join(os.path.split(os.path.split(__file__)[0])[0], "resources")

    def test_something(self):
        convert_ped_to_bgl(os.path.join(self.resource_path, "HapMap3_CEU_chr2.ped"),
                           os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"),
                           os.path.join(self.resource_path, "HapMap3_CEU_chr2.bgl"),
                           os.path.join(self.resource_path, "HapMap3_CEU_chr2.markers"))
        x = 1

if __name__ == '__main__':
    unittest.main()
