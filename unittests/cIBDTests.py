__author__ = 'Itamar'

import unittest
import os
from IBD.cIBD import cPairIBD, cPopulationIBD
from IBD.GeneticMap import GeneticMap

class MyTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.resource_path = os.path.join(os.path.split(os.path.split(__file__)[0])[0], "resources")

    def test_add_interval(self):
        p = cPairIBD()
        p.add_interval(0,100,3)
        self.assertEqual(p.to_list(), [(0,100,3)], "adding interval does not work")

    def test_to_from_list(self):
        p = cPairIBD()
        p.add_interval(0,100,3)
        p.add_interval(50,150,4)
        self.assertEqual([(0,100,3),(50,150,4)], (cPairIBD.from_list(p.to_list())).to_list())

    def test_to_from_string(self):
        p = cPairIBD()
        p.add_interval(0,100,3)
        p.add_interval(50,150,4)
        pop = cPopulationIBD()
        pop.add_human_pair((1,2),p)
        p = cPairIBD()
        p.add_interval(200,300,5)
        p.add_interval(400,500,6)
        pop.add_human_pair((1,3),p)
        self.assertEqual(pop.to_list(), (cPopulationIBD.from_string(pop.to_string())).to_list())

    def test_merge_intervals(self):
        p = cPairIBD()
        p.add_interval(0,10,1)
        p.add_interval(20,30,2)
        p.add_interval(5,25,2)
        p.merge_intervals_fast(merge_diff_vals=True)
        self.assertEqual(p.to_list(), [(0,30,1)])

        p = cPairIBD()
        p.add_interval(0,10,1)
        p.add_interval(20,30,2)
        p.add_interval(5,25,2)
        p.merge_intervals_fast(merge_diff_vals=False)
        self.assertEqual(p.to_list(), [(0,10,1),(5,30,2)])

        p = cPairIBD()
        p.add_interval(20,30,2)
        p.add_interval(5,25,2)
        p.add_interval(0,10,1)
        p.merge_intervals_fast(merge_diff_vals=True)
        self.assertEqual(p.to_list(), [(0,30,1)])

        p2 = cPairIBD()
        p2.add_interval(0,30,1)
        p2.add_interval(5,15,2)
        p2.add_interval(20,25,3)
        p2.merge_intervals_fast(merge_diff_vals=True)
        self.assertEqual(p2.to_list(), [(0,30,1)])

    def test_filter_by_length(self):
        p = cPairIBD()
        p.add_interval(0,10,1)
        p.add_interval(20,30,2)
        p.add_interval(5,25,2)
        gm = GeneticMap(os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"), max_snp_num=100000)

        p.filter_by_length(1,2,gm)
        self.assertEqual(p.to_list(), [(0,30,1)])

if __name__ == '__main__':
    unittest.main()
