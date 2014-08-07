__author__ = 'Itamar'

import unittest
from utils.GermlineWrapper import germline


class MyTestCase(unittest.TestCase):
    def test_germline(self):
        germline(['germline','-prefix', '../resources/ceu.genos', '-bits', '64', '-min_m', '0.1', '-err_hom', '4', '-err_het', '2', '-map', '../resources/ceu.genos.map', '-w_extend'])


if __name__ == '__main__':
    unittest.main()
