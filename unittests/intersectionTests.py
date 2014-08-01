__author__ = 'Itamar'

import unittest
from IBD.IntervalTree import Interval, IntervalTree

class MyTestCase(unittest.TestCase):
    def test_find(self):
        t = IntervalTree()
        inter = Interval(0,100,3)
        t.insert_interval(inter)
        intersections = t.find(50,150)
        self.assertEqual(intersections, [inter])

    def test_tolist(self):
        t = IntervalTree()
        inter = Interval(0,100,3)
        t.insert_interval(inter)
        l = t.to_list()
        self.assertEqual(l, [inter])

if __name__ == '__main__':
    unittest.main()
