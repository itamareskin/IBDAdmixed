'''
Created on Jan 13, 2012

@author: itamares
'''

from bx.intervals import IntervalTree,Interval
import Infra.IntervalUtils as iu 

class FoundersContainer(object):
    
    def __init__(self, start_snp, snp_num, founder=None):
        self._start_snp = start_snp
        self._snp_num = snp_num
        self._founder_to_tree = {}
        if founder != None:
            self._founder_to_tree[founder] = IntervalTree()
            self._founder_to_tree[founder].add_interval(Interval(start_snp,start_snp+snp_num))
    
    @property 
    def start_snp(self):
        return self._start_snp
    
    @property 
    def snp_num(self):
        return self._snp_num
    
    def set_founder_in_interval(self, founder, start, end):
        if not self._founder_to_tree.has_key(founder):
            self._founder_to_tree[founder] = IntervalTree()
        founders_to_remove = []
        for other_founder in self._founder_to_tree.keys():
            if other_founder == founder:
                continue
            new_tree = iu.exclude_interval_from_tree(self._founder_to_tree[other_founder], start, end)
            if len(iu.get_interval_list(new_tree)) > 0:
                self._founder_to_tree[other_founder] = new_tree
            else:
                if other_founder != founder:
                    founders_to_remove.append(other_founder)
        for f in founders_to_remove:
            self._founder_to_tree.pop(f)
        self._founder_to_tree[founder].add_interval(Interval(start,end))
        
    def get_founders(self):
        return self._founder_to_tree.keys()
    
    def get_founders_in_interval(self, start, end):
        l = []
        for founder in self._founder_to_tree.keys():
            intersections = self._founder_to_tree[founder].find(start,end)
            if len(intersections) > 0:
                l.append(founder)
        return l
    
    def get_intervals_of_founder_in_interval(self, founder, start, end):
        if not self._founder_to_tree.has_key(founder):
            return []
        l = []
        intersections = self._founder_to_tree[founder].find(start,end)
        for inter in intersections:
            l.append(Interval(max(inter.start,start),min(inter.end,end)))
        return l
    
    def get_total_length_of_founder(self, founder, start=None, end=None):
        if start == None:
            start = self.start_snp
        if end == None:
            end = self.start_snp + self.snp_num
        l = self.get_intervals_of_founder_in_interval(founder, start, end)
        total_length = 0
        for interval in l:
            total_length += interval.end-interval.start
        return total_length
    
    def copy_interval(self, founder_container2, start, end):
        founders = founder_container2.get_founders_in_interval(start,end)
        for founder in founders:
            intervals = founder_container2.get_intervals_of_founder_in_interval(founder,start,end)
            for inter in intervals:
                self.set_founder_in_interval(founder, inter.start, inter.end)
    
    def get_window(self, start, end):
        new_founder_container = FoundersContainer(self._start_snp, self._snp_num)
        new_founder_container.copy_interval(self, start, end)
        return new_founder_container
    
    def to_dict(self):
        d = {}
        for founder in self._founder_to_tree.keys():
            d[founder] = [] 
            self._founder_to_tree[founder].traverse(lambda x: d[founder].append((x.start,x.end)))
        return d
    
    def merge_all(self):
        founders = self.get_founders() 
        for f in founders:
            self._founder_to_tree[f] = iu.merge_intervals(self._founder_to_tree[f])
    
    @staticmethod
    def from_dict(d,start_snp,snp_num):
        f = FoundersContainer(start_snp,snp_num)
        for founder in d.keys():
            f._founder_to_tree[founder] = IntervalTree()
            for inter in d[founder]:
                f._founder_to_tree[founder].add_interval(Interval(inter[0],inter[1]))
        return f
    
    def to_serialize_string(self):
        s = []
        for founder in self._founder_to_tree.keys():
            l = [] 
            self._founder_to_tree[founder].traverse(lambda x: l.append((x.start,x.end)))
            if len(l)>0:
                s += [str(founder),":"]
                for inter in l:
                    s += [str(inter[0]),",",str(inter[1]),";"]
                s.pop()
                s.append("|")
        if len(s) > 0:
            s.pop() 
        return ''.join(s)
    
    @staticmethod
    def from_serialize_string(s,start_snp,snp_num):
        f = FoundersContainer(start_snp,snp_num)
        s = s.split("|")
        for founder_string in s:
            temp = founder_string.split(":")
            founder = int(temp[0])
            f._founder_to_tree[founder] = IntervalTree()
            intervals = temp[1].split(";")
            for inter in intervals:
                points = inter.split(",")
                f._founder_to_tree[founder].add_interval(Interval(int(points[0]),int(points[1])))
        return f
    
    @staticmethod
    def from_flat_string(s,start_snp,snp_num):
        f = FoundersContainer(start_snp,snp_num)
        s += '$'
        last_founder = s[0]
        curr_start = start_snp
        for i in range(len(s)):
            curr_founder = s[i]
            if curr_founder != last_founder:
                if not f._founder_to_tree.has_key(int(last_founder)):
                    f._founder_to_tree[int(last_founder)] = IntervalTree()
                f._founder_to_tree[int(last_founder)].add_interval(Interval(curr_start, start_snp+i))
                curr_start = start_snp+i
                last_founder = curr_founder
        return f
                
def main():
    
    f = FoundersContainer.from_flat_string("000111100010001111", 18)
    d = f.to_dict()
    t1 = IntervalTree()
    t1.add_interval(Interval(0,20))
    t1.add_interval(Interval(30,40))
    t1.add_interval(Interval(50,60))
    
    r = iu.exclude_interval_from_tree(t1,35,37)
    
    print r.find(-100,100)


if __name__ == "__main__":
     main()