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

    def test_intersect(self):
        t = IntervalTree()
        t.insert_interval(Interval(0,100,1))
        t.insert_interval(Interval(200,300,2))
        t.insert_interval(Interval(350,400,3))
        other = IntervalTree()
        other.insert_interval(Interval(150,250,1))
        intersections = t.intersect(other)
        self.assertEqual(intersections, [(200,250)])

if __name__ == '__main__':
    unittest.main()
