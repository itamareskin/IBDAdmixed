#cython: profile=True
#cython: boundscheck=False
#cython: cdivision=True
#cython.wraparound=False
#cython.nonecheck=False

'''
Created on Jan 11, 2012

@author: itamares
'''
from __future__ import division
from IBD.IntervalTree cimport Interval, IntervalTree
from IBD.IntervalTree import Interval, IntervalTree
import numpy
import string
from libcpp.vector cimport vector as cpp_vector
from GeneticMap cimport GeneticMap

cdef class PairIBD:
    '''
    Represents segments of IBD shared by a pair of individuals
    Implemented using a IntervalTree to allow fast intersection queries
    '''
    
    cdef public IntervalTree _tree
    
    def __cinit__(self):
        self._tree = IntervalTree()

    property _list:
        def __get__(PairIBD self):
            return self.to_list()
        
    cpdef add_interval(PairIBD self, int start, int end, float value=0):
        self._tree.insert_interval(Interval(start,end,value))
    
    cpdef find(PairIBD self,int start,int end):
        return self._tree.find(start,end)
    
    cpdef update(PairIBD self, PairIBD other_ibd, max_val=False):
        cdef list intervals = other_ibd.to_list()
        for interval in intervals:
            self.add_interval(interval[0], interval[1], interval[2])
        self.merge_intervals(max_val=max_val)
    
    cpdef to_string(PairIBD self):
        s = []
        intervals = self.to_list()
        for inter in intervals:
            s += [str(inter[0]),",",str(inter[1]),",",str(inter[2]),";"]
        if len(intervals)>0:
            s.pop()
        return "".join(s)
    
    cpdef list to_list(PairIBD self):
        return [(x.start,x.end,x.value) for x in self._tree.to_list()]
    
    @classmethod
    def from_list(cls, list l):
        p = PairIBD()
        for interval in l:
            if len(interval) > 2:
                p.add_interval(interval[0],interval[1],interval[2])
            else:
                p.add_interval(interval[0],interval[1])
        return p

    cpdef merge_intervals(PairIBD self, int ovelap = 1, max_val=False, merge_diff_vals=False):

        # get a list of all intervals in the IntervalTree
        intervals = self._tree.to_list()
        if len(intervals) == 0:
            return

        # clear the IntervalTree
        self._tree = IntervalTree()

        # sort the intervals according to their value
        decorated = [(tup.value, tup) for tup in intervals]
        if max_val:
            decorated.sort()
        else:
            decorated.sort(reverse=True)
        intervals = [tup for second, tup in decorated]

        for interval in intervals:
            intersections = self._tree.find(interval.start-ovelap,interval.end+ovelap)
            if len(intersections) > 0:
                non_intersecting = self._tree.to_list()
                for curr_intersection in intersections:
                    non_intersecting.remove(curr_intersection)

                self._tree = IntervalTree()
                for non in non_intersecting:
                    self._tree.insert_interval(non)

                new_intervals = []
                for curr_inter in intersections:
                    if curr_inter.value == interval.value or merge_diff_vals:
                        interval=interval.merge(curr_inter,max_val)
                    else:
                        if (not max_val and curr_inter.value > interval.value) or (max_val and curr_inter.value < interval.value):
                            new_intervals += curr_inter.substract(interval)

                new_intervals.append(interval)
                for new_int in new_intervals:
                    self._tree.insert_interval(new_int)
            else:
                self._tree.insert_interval(interval)

    cpdef merge_intervals_fast(PairIBD self, int overlap = 1, max_val=False, merge_diff_vals=False):

        # get a list of all intervals in the IntervalTree
        intervals = self._tree.to_list()
        if len(intervals) == 0:
            return

        # find the maximal end position of all intervals
        max_end = max([x.end for x in intervals])

        # sort the intervals according to their start position
        decorated = [(tup.start, tup) for tup in intervals]
        decorated.sort()
        intervals = [tup for second, tup in decorated]

        # clear the IntervalTree
        self._tree = IntervalTree()

        # these will hold the start,end and value of the current interval that is being created
        prev_score = -1e100
        prev_start = intervals[0].start-overlap-1
        prev_end = intervals[0].start-overlap-1

        # go over all intervals, plus one dummy Interval that just makes the merged intervals before it get inserted to the tree
        for interval in intervals+[Interval(max_end+overlap+1,max_end+overlap+1,0)]:
            # if there is an overlap, update the end and score (start position remains the same)
            if interval.start - prev_end <= overlap and (interval.value == prev_score or merge_diff_vals):
                prev_end = max(prev_end,interval.end)
                if max_val:
                    prev_score = <float>max(prev_score,interval.value)
                else:
                    prev_score = <float>min(prev_score,interval.value)
            # if there is not overlap, insert the previous merge results to the tree and update the start,end and value to the current interval
            else:
                if prev_start != intervals[0].start-overlap-1:
                    self._tree.insert_interval(Interval(prev_start,prev_end,prev_score))
                prev_start = interval.start
                prev_end = interval.end
                prev_score = interval.value
    
    cpdef get_num_windows(self, int window_size = 1):
        l = self.to_list()
        total_length = 0
        for segment in l:
            total_length += round((segment[1]-segment[0])/float(window_size))
        return total_length
    
    cpdef get_IBD_percent(PairIBD self, GeneticMap gm):
        return  100 * float(sum([i[1]-i[0] for i in self.to_list()])) / gm._snp_num
    
    cpdef stats_win(PairIBD self, PairIBD true_ibd, GeneticMap gm, int window_size = 1):
        TP = 0
        if true_ibd is not None:
            intersections = self._tree.intersect(true_ibd._tree)
            if len(intersections) > 0:
                for inter in intersections:
                    TP += round((inter[1] - inter[0])/float(window_size))
        FP = self.get_num_windows(window_size) - TP
        FN = 0
        if true_ibd is not None:
            FN = true_ibd.get_num_windows(window_size) - TP
        TN = gm._snp_num/window_size - (TP + FP + FN)
        power = TP/(TP+FN) if TP+FN>0 else 0
        specificity = TN/(TN+FP) if TN+FP>0 else 1
        FPR = 1 - specificity
        if TP+FP == 0:
            FDR = 0
        else:
            FDR = FP/(TP+FP)
        return {'TP': TP, 'FP': FP, 'TN': TN, 'power': power, 'FDR': FDR, 'FPR': FPR}
    
    cpdef get_detected_segments(PairIBD self, PairIBD true_ibd):
        true_segments = true_ibd.to_list()
        cdef list detected = []
        tree = self._tree
        for segment in true_segments:
            intersections = tree.find(segment[0],segment[1])
            if len(intersections) > 0:
                detected.append(segment)
        return detected
    
    cpdef add_ibd_from_containers(PairIBD self, f1, f2, min_ibd_length = 0):
        for founder in f1._founder_to_tree.keys():
            if f2._founder_to_tree.has_key(founder):
                intersections =\
                f1._founder_to_tree[founder].intersect(f2._founder_to_tree[founder])
                for inter in intersections:
                    if inter[1]-inter[0] > min_ibd_length:
                        self.add_interval(inter[0], inter[1])
        self.merge_intervals_fast()
        
    cpdef add_ibd_from_containers_only_admixed(PairIBD self, f1, f2, min_ibd_length, f1_other, f2_other, ceu_inds, yri_inds):
        for founder in f1._founder_to_tree.keys():
            if f2._founder_to_tree.has_key(founder):
                intersections =\
                f1._founder_to_tree[founder].intersect(f2._founder_to_tree[founder])
                for inter in intersections:
                    if inter[1]-inter[0] > min_ibd_length:
                        f1_other_founders = f1.get_founders_in_interval(inter[0], inter[1])
                        f2_other_founders = f2.get_founders_in_interval(inter[0], inter[1])
                        anc = int(f1 in yri_inds)
                        for other in f1_other_founders:
                            print other, int(other in yri_inds), int(other in yri_inds)
                            if int(other in yri_inds) != int(founder in yri_inds):
                                self.add_interval(inter[0], inter[1])
                        for other in f2_other_founders:
                            print other, int(other in yri_inds), int(other in yri_inds)
                            if int(other in yri_inds) != int(founder in yri_inds):
                                self.add_interval(inter[0], inter[1])
        self.merge_intervals_fast()
        
    cpdef filter_by_length(self, float min_length, float max_length, GeneticMap gm):
        cdef IntervalTree new_tree = IntervalTree()
        cdef float length
        cdef list intervals = self.to_list()
        for interval in intervals:
            length = gm.get_length(interval[0],interval[1])
            if length >= min_length and length <= max_length:
                new_tree.insert_interval(Interval(interval[0],interval[1],interval[2]))
        self._tree = new_tree
        
    cpdef filter_by_score(self, min_score, max_score):
        cdef IntervalTree new_tree = IntervalTree()
        cdef list intervals = self.to_list()
        for interval in intervals:
            if interval[2] >= min_score and interval[2] <= max_score:
                new_tree.insert_interval(Interval(interval[0],interval[1],interval[2]))
        self._tree = new_tree
        
    cpdef filter_by_other_ibd(self, PairIBD other_ibd):
        cdef IntervalTree new_tree = IntervalTree()
        cdef list intervals = self.to_list()
        for interval in intervals:
            if len(other_ibd.find(interval[0],interval[1])) > 0:
                new_tree.insert_interval(Interval(interval[0],interval[1],interval[2]))
        self._tree = new_tree

    cpdef filter_out_other_ibd(self, PairIBD other_ibd):
        cdef IntervalTree new_tree = IntervalTree()
        cdef list intervals = self.to_list()
        for interval in intervals:
            if len(other_ibd.find(interval[0],interval[1])) == 0:
                new_tree.insert_interval(Interval(interval[0],interval[1],interval[2]))
        self._tree = new_tree

cdef class PopIBD:
    
    cdef public dict _ibd_dic
    
    def __cinit__(self):
        self._ibd_dic = {}

    cpdef add_human_pair(self, pair, pairIBD):
        if self.has_key(pair):
            self._ibd_dic[pair].update(pairIBD)
        else:
            self._ibd_dic[pair] = pairIBD
              
    cpdef add_interval_to_pair(self, pair, start, end, value=0):
        self._ibd_dic[pair].add_interval(start,end,value)
        
    cpdef clear(self):
        for pair in self._ibd_dic.keys():
            self._ibd_dic[pair].clear()
            del self._ibd_dic[pair]
        self._ibd_dic.clear()
              
    cpdef has_key(self, pair):
        return self._ibd_dic.has_key(pair)
    
    cpdef get_value(self, pair):
        if self._ibd_dic.has_key(pair):
            return self._ibd_dic[pair]
        else:
            return None
    
    cpdef update(self, other_ibd):
        assert isinstance(other_ibd, PopIBD)
        for pair in other_ibd.keys():
            if self.has_key(pair):
                self._ibd_dic[pair].update(other_ibd.get_value(pair))
            else:
                self.add_human_pair(pair, other_ibd.get_value(pair))
    
    cpdef keys(self):
        return self._ibd_dic.keys()
        
    cpdef merge_all(self, int overlap = 1, max_val = False, merge_diff_vals=False):
        for pair in self.keys():
            self.get_value(pair).merge_intervals(overlap, max_val, merge_diff_vals)
            
    cpdef merge_all_fast(self, int overlap = 1, max_val = False, merge_diff_vals=False):
        for pair in self.keys():
            self.get_value(pair).merge_intervals_fast(overlap, max_val, merge_diff_vals)
    
    cpdef stats_win(self, true_ibd, GeneticMap gm, window_size = 1):
        stats = []
        for pair in true_ibd.keys():
            if self.has_key(pair):
                stats.append(self.get_value(pair).stats_win(true_ibd.get_value(pair),gm,window_size))
            else:
                #TN = gm._snp_num/window_size - true_ibd.get_value(pair).get_num_windows(window_size)
                #FN = true_ibd.get_value(pair).get_num_windows(window_size)
                TP = 0
                FP = 0
                FN = true_ibd.get_value(pair).get_num_windows(window_size) - TP
                TN = gm._snp_num/window_size - (TP + FP + FN)
                power = TP/(TP+FN) if TP+FN>0 else 0
                specificity = TN/(TN+FP) if TN+FP>0 else 1
                FPR = 1 - specificity
                if TP+FP == 0:
                    FDR = 0
                else:
                    FDR = FP/(TP+FP)
                stats.append({'TP': TP, 'FP': FP, 'TN': TN, 'power': power, 'FDR': FDR, 'FPR': FPR})
        # calc stats
        TP = numpy.sum([x['TP'] for x in stats])
        FP = numpy.sum([x['FP'] for x in stats])
        TN = numpy.sum([x['TN'] for x in stats])
        power = numpy.mean([x['power'] for x in stats])
        FDR = numpy.mean([x['FDR'] for x in stats])
        FPR = numpy.mean([x['FPR'] for x in stats])
        #with open("stats.txt","w") as f:
        #    f.writelines(string.join([str(x) for x in stats],"\n"))
        return {'TP': TP, 'FP': FP, 'TN': TN, 'power': power, 'FDR': FDR, 'FPR': FPR}
    
    cpdef get_all_segments(self):
        segments = []
        for pair in self._ibd_dic.keys():
            segments.append(self.get_value(pair).to_list())
        return segments 

    cpdef calc_power(self, true_ibd):
        assert isinstance(true_ibd, PopIBD)
        detected_all = []
        detected_dict = {}
        for pair in true_ibd.keys():
            if self.has_key(pair):
                detected = self._ibd_dic[pair].get_detected_segments(true_ibd.get_value(pair))
                detected_all += detected
                if len(detected) > 0:
                    detected_dict[pair] = detected
          
        return (float(len(detected_all)) / len(true_ibd.get_all_segments()), detected_dict)
         
    cpdef filter_by_human_ids(self, human_ids):
        for pair in self.keys():
            if pair[0] not in human_ids or pair[1] not in human_ids:
                self._ibd_dic.pop(pair)
    
    cpdef filter_by_human_ids_partial(self, human_ids):
        for pair in self.keys():
            if pair[0] not in human_ids and pair[1] not in human_ids:
                self._ibd_dic.pop(pair)
    
    cpdef filter_by_human_pairs(self, pairs):
        for pair in self.keys():
            if pair not in pairs:
                self._ibd_dic.pop(pair)
                
    cpdef filter_by_length(self, float min_length, float max_length, GeneticMap gm):
        for pair in self.keys():
            self._ibd_dic[pair].filter_by_length(min_length, max_length, gm)
            if len(self._ibd_dic[pair].to_list()) == 0:
                self._ibd_dic.pop(pair)
                
    cpdef filter_by_score(self, min_score, max_score):
        for pair in self.keys():
            self._ibd_dic[pair].filter_by_score(min_score,max_score)
            if len(self._ibd_dic[pair].to_list()) == 0:
                self._ibd_dic.pop(pair)
                
    cpdef filter_by_other_ibd(self, PopIBD other_ibd):
        for pair in self.keys():
            if pair in other_ibd.keys():
                self._ibd_dic[pair].filter_by_other_ibd(other_ibd.get_value(pair))
                if len(self._ibd_dic[pair].to_list()) == 0:
                    self._ibd_dic.pop(pair)
            else:
                self._ibd_dic.pop(pair)

    cpdef filter_out_other_ibd(self, PopIBD other_ibd):
        for pair in self.keys():
            if pair in other_ibd.keys():
                self._ibd_dic[pair].filter_out_other_ibd(other_ibd.get_value(pair))
                if len(self._ibd_dic[pair].to_list()) == 0:
                    self._ibd_dic.pop(pair)

    cpdef to_dict(self):
        d = {}
        for pair in self._ibd_dic.keys():
            d[pair] = self._ibd_dic[pair].to_list()
        return d
    
    @classmethod
    def from_dict(cls, d):
        p = PopIBD()
        for pair in d.keys():
            p.add_human_pair(pair,PairIBD.from_list(d[pair]))
        return p
    
    cpdef to_list(self):
        l = []
        for pair in self._ibd_dic.keys():
            l += self._ibd_dic[pair].to_list()
        return l
    
    cpdef to_string(self):
        s = []
        for pair in self._ibd_dic.keys():
            s += [str(pair[0]),",",str(pair[1]),":"]
            intervals = self._ibd_dic[pair].to_list()
            for inter in intervals:
                s += [str(inter[0]),",",str(inter[1]),",",str(inter[2]),";"]
            if len(intervals)>0:
                s.pop()
            s.append("\n")
        return "".join(s)
    
    @classmethod
    def from_string(cls, s):
        assert isinstance(s, (str,list))
        p = PopIBD()
        if type(s) is str:
            s = string.split(s, '\n')
            s = [x for x in s if x != '']
        s = [string.strip(x, '\n') for x in s]
        for pair_string in s:
            temp = pair_string.split(":")
            pair = temp[0].split(",")
            pair = (int(pair[0]),int(pair[1]))
            pairibd = PairIBD()
            for inter in temp[1].split(";"):
                points = inter.split(",")
                score = 0
                if len(points) > 2:
                    score = float(points[2])
                pairibd.add_interval(int(points[0]),int(points[1]),score)
            p.add_human_pair(pair,pairibd)
        return p
    
    @classmethod
    def fast_deserialize(cls, file_name):
        f = open(file_name)
        s = f.readlines()
        ibd = PopIBD.from_string(s)
        f.close()
        return ibd
    
    