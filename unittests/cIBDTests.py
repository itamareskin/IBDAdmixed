__author__ = 'Itamar'

import unittest
import os
from IBD.cIBD import cPairIBD, cPopulationIBD
from IBD.GeneticMap import GeneticMap
from IBD.IntervalTree import Interval

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
        p.merge_intervals(merge_diff_vals=True)
        self.assertEqual(p.to_list(), [(0,30,1)])

        p = cPairIBD()
        p.add_interval(0,10,1)
        p.add_interval(20,30,2)
        p.add_interval(5,25,2)
        p.merge_intervals(merge_diff_vals=False)
        self.assertEqual(p.to_list(), [(0,10,1),(10,30,2)])

        p = cPairIBD()
        p.add_interval(0,30,1)
        p.add_interval(5,15,2)
        p.add_interval(20,25,3)
        p.merge_intervals(merge_diff_vals=False, max_val=True)
        self.assertEqual(p.to_list(), [(0,5,1),(5,15,2),(15,20,1),(20,25,3),(25,30,1)])

        p = cPairIBD()
        p.add_interval(20,30,2)
        p.add_interval(5,25,2)
        p.add_interval(0,10,1)
        p.merge_intervals(merge_diff_vals=True)
        self.assertEqual(p.to_list(), [(0,30,1)])

        p2 = cPairIBD()
        p2.add_interval(0,30,1)
        p2.add_interval(5,15,2)
        p2.add_interval(20,25,3)
        p2.merge_intervals(merge_diff_vals=True)
        self.assertEqual(p2.to_list(), [(0,30,1)])

    def test_merge_intervals_fast(self):
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
        p.add_interval(0,1000,1)
        p.add_interval(2000,2100,2)
        p.add_interval(3000,3300,2)
        gm = GeneticMap(os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"), max_snp_num=100000)

        p.filter_by_length(1,2,gm)
        self.assertEqual(p.to_list(), [(3000,3300,2)])

    def test_filter_by_score(self):
        p = cPairIBD()
        p.add_interval(0,1000,1)
        p.add_interval(2000,2100,2)
        p.add_interval(3000,3300,3)

        p.filter_by_score(1.5,2.5)
        self.assertEqual(p.to_list(), [(2000,2100,2)])

    def test_find(self):
        p = cPairIBD()
        p.add_interval(0,1000,1)
        p.add_interval(2000,2100,2)
        p.add_interval(3000,3300,3)

        intervals = p.find(1800,2300)
        self.assertEqual(intervals, [Interval(2000,2100,2)])

    def test_update(self):
        p = cPairIBD()
        p.add_interval(0,1000,1)
        p.add_interval(2000,2100,2)
        p.add_interval(3000,3300,3)

        other = cPairIBD()
        other.add_interval(800,1200,2)

        p.update(other, max_val=True)
        self.assertEqual(p.to_list(), [(0,800,1),(800,1200,2),(2000,2100,2),(3000,3300,3)])

    def test_stats_win(self):
        gm = GeneticMap(os.path.join(self.resource_path, "HapMap3_CEU_chr2.map"), max_snp_num=100000)

        p = cPairIBD()
        p.add_interval(0,1000,1)
        p.add_interval(2000,2100,2)
        p.add_interval(3000,3300,3)
        true_ibd = cPairIBD()
        true_ibd.add_interval(0,1000,1)
        true_ibd.add_interval(2000,2100,2)
        true_ibd.add_interval(3000,3300,3)
        stats = p.stats_win(true_ibd,gm)
        self.assertEqual(stats, {'FP': 0.0, 'power': 1.0, 'TP': 1400.0, 'TN': 98600.0, 'FDR': 0.0, 'FPR': 0.0})

        true_ibd = cPairIBD()
        true_ibd.add_interval(500,1000,1)
        true_ibd.add_interval(3200,3500,3)
        stats = p.stats_win(true_ibd,gm)
        self.assertEqual(stats, {'FP': 800.0, 'power': 0.75, 'TP': 600.0, 'TN': 98400.0, 'FDR': 0.5714285714285714, 'FPR': 0.008064516129032251})

        p.add_interval(0,1000,2)
        p.add_interval(2000,2100,3)
        p.add_interval(3000,3300,4)
        p.merge_intervals(merge_diff_vals=True)
        stats = p.stats_win(true_ibd,gm)
        self.assertEqual(stats, {'FP': 800.0, 'power': 0.75, 'TP': 600.0, 'TN': 98400.0, 'FDR': 0.5714285714285714, 'FPR': 0.008064516129032251})

        p = cPairIBD()
        p.add_interval(0,1000,2)
        p.add_interval(900,2000,3)
        p.add_interval(1900,3300,4)
        p.merge_intervals(merge_diff_vals=False)
        stats = p.stats_win(true_ibd,gm)
        self.assertEqual(stats, {'FP': 2700.0, 'power': 0.75, 'TP': 600.0, 'TN': 96500.0, 'FDR': 0.8181818181818182, 'FPR': 0.027217741935483875})



if __name__ == '__main__':
    unittest.main()
