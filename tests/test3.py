__author__ = 'Itamar'
from IBD.IBDSegments import PairIBD
from IBD.IntervalTree import Interval,IntervalTree

def merge_intervals_fast(p, overlap = 1, max_val=False, merge_diff_vals=False):
    intervals = p._tree.to_list()
    if len(intervals) == 0:
        return
    decorated = [(tup.start, tup) for tup in intervals]
    decorated.sort()
    max_end = max([x.end for x in intervals])
    intervals = [tup for second, tup in decorated]
    p._tree = IntervalTree()
    prev_score = -1e100
    prev_start = intervals[0].start-overlap-1
    prev_end = intervals[0].start-overlap-1

    for interval in intervals+[Interval(max_end+overlap+1,max_end+overlap+1,0)]:
        if interval.start - prev_end <= overlap and (interval.value == prev_score or merge_diff_vals):
            prev_end = max(prev_end,interval.end)
            if max_val:
                prev_score = max(prev_score,interval.value)
            else:
                prev_score = min(prev_score,interval.value)
        else:
            if prev_start != intervals[0].start-overlap-1:
                p._tree.insert_interval(Interval(prev_start,prev_end,prev_score))
            prev_start = interval.start
            prev_end = interval.end
            prev_score = interval.value

p2 = PairIBD()
p2.add_interval(0,30,1)
p2.add_interval(5,15,2)
p2.add_interval(20,25,3)
merge_intervals_fast(p2,merge_diff_vals=True)
x=1