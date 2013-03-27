
import copy_reg
from bx.intervals import IntervalTree, Interval
#from Infra.IntervalTree import IntervalTree, Interval

def IntersectIntervalTrees(tree,interval_list):
    intersections = []
    for interval in interval_list:
        intersections += [(max(interval[0], i.start), min(interval[1], i.end)) for i in tree.find(interval[0],interval[1])]
    return intersections
    
    
def get_interval_list(tree):
    l = []; 
    tree.traverse(lambda x: l.append((x.start,x.end)))
    return l

def exclude_interval_from_tree(tree,start,end):
    new_tree = IntervalTree()
    intersecting = tree.find(start,end)
    if len(intersecting) > 0:
        before_intervals = tree.before(start)
        for b in before_intervals:
            new_tree.add_interval(b)
        after_intervals = tree.after(end)
        for a in after_intervals:
            new_tree.add_interval(a)
        for inter in intersecting:
            l = substract_intervals(inter, Interval(start,end))
            for new_int in l:
                new_tree.add_interval(new_int)
    else:
        new_tree = tree
    return new_tree

def substract_intervals(interval1,interval2):
    result = []
    if interval1.start < interval2.start:
        result.append(Interval(interval1.start,interval2.start))
        if interval1.end > interval2.end:
            result.append(Interval(interval2.end,interval1.end))
    else:
        if interval2.end < interval1.end:
            result.append(Interval(interval2.end,interval1.end))
    return result

def merge_intervals(tree):
    new_tree = IntervalTree()
    intervals = get_interval_list(tree)
    for interval in intervals:
        intersections = [(x.start,x.end) for x in new_tree.find(interval[0]-1,interval[1]+1)]
        if len(intersections) > 0:
            non_intersecting = []; 
            new_tree.traverse(lambda x: non_intersecting.append((x.start,x.end)))
            for intersection in intersections:
                non_intersecting.remove(intersection)
                
            new_tree = IntervalTree()
            for non in non_intersecting:
                new_tree.add_interval(Interval(non[0],non[1]))
            
            intersections.append(interval)
            new_tree.add_interval(Interval(min([x[0] for x in intersections]),max([x[1] for x in intersections])))
        else:
            new_tree.add_interval(Interval(interval[0],interval[1]))
    return new_tree;

def main():
    
    t1 = IntervalTree()
    t2 = IntervalTree()
    t1.add_interval(Interval(0,20))
    t1.add_interval(Interval(30,32))
    t1.add_interval(Interval(33,40))
    
    t2.add_interval(Interval(10,13))
    t2.add_interval(Interval(14,17))
    t2.add_interval(Interval(26,35))
    
    r = IntersectIntervalTrees(t1,get_interval_list(t2))
    
    print r


if __name__ == "__main__":
     main()