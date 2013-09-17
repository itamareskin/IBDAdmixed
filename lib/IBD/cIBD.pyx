'''
Created on Jan 11, 2012

@author: itamares
'''
from __future__ import division
import cPickle
from bx.intervals import Interval, IntervalTree
import IntervalUtils as iu 
import os
import copy
import numpy
import string
from libcpp.pair cimport pair as cpp_pair
from libcpp.vector cimport vector as cpp_vector
#from libcpp.map cimport map as cpp_map 

cdef extern from "intervalItem.h":
    
    cdef struct struct_interval_item:
        int start
        int end
        float value
        float length

    ctypedef struct_interval_item interval_item
    
    interval_item create_interval_item(int start, int end, float value, float length)

cdef class cPairIBD:
    '''
    Represents segments of IBD shared by a pair of individuals
    Implemented using a IntervalTree to allow fast intersection queries
    '''
    
    cdef cpp_vector[interval_item] *_intervals
    
    def __cinit__(self): 
        self._intervals = new cpp_vector[interval_item]()
        #self._tree = IntervalTree()
        #if intervals != None:
        #    for interval in intervals: 
        #        self.add_interval(interval[0],interval[1])

    property _tree:
        def __get__(cPairIBD self):
            return self.__get_interval_tree()
    
    cdef __get_interval_tree(cPairIBD self):
        tree = IntervalTree()
        #for interval in self._intervals:
        #    tree.add_interval(Interval(interval.first,interval.second))   
        #cdef cpp_list[cpp_pair[int,int]].iterator i = self._intervals.begin()
        cdef interval_item p 
        #while i != self._intervals.end():
        for i in range(self._intervals.size()):
            p = self._intervals.at(i)
            tree.add_interval(Interval(<int>p.start,<int>p.end)) 
            #tree.insert(p.start,p.end,p.value)
        return tree
    
    cdef __get_values_interval_tree(cPairIBD self):
        tree = IntervalTree()
        #for interval in self._intervals:
        #    tree.add_interval(Interval(interval.first,interval.second))   
        #cdef cpp_list[cpp_pair[int,int]].iterator i = self._intervals.begin()
        cdef interval_item p 
        #while i != self._intervals.end():
        for i in range(self._intervals.size()):
            p = self._intervals.at(i)
            #tree.add_interval(Interval(<int>p.start,<int>p.end)) 
            tree.insert(p.start,p.end,p.value)
        return tree
        
    cpdef add_interval(cPairIBD self, int start, int end, float value=0, float length=0):
        cdef interval_item p 
        p.start = start 
        p.end = end
        p.value = value
        p.length = length
        self._intervals.push_back(p)
        #self._tree.add_interval(Interval(start,end))
    
    cpdef clear(cPairIBD self):
        del self._intervals
    
    cpdef find(cPairIBD self,int start,int end):
        intervals = self._tree.find(start,end)
        return self._tree.find(start,end)
    
    cpdef find_with_values(cPairIBD self,int start,int end):
        intervals = self._tree.find(start,end)
        values = self.__get_values_interval_tree().find(start,end)
        return zip(intervals,values)
    
    cpdef update(cPairIBD self, cPairIBD other_ibd):
        cdef list intervals = other_ibd.to_list()
        for interval in intervals:
            self.add_interval(interval[0], interval[1])
        self.merge_intervals()
    
    cpdef to_string(cPairIBD self):
        s = []
        intervals = self.to_list()
        for inter in intervals:
            s += [str(inter[0]),",",str(inter[1]),";"]
        if len(intervals)>0:
            s.pop()
        return "".join(s)
    
    cpdef list to_list(cPairIBD self):
    #    return self._intervals
        #l = []; self._tree.traverse(lambda x: l.append((x.start,x.end)))
        #return l
        cdef list result = [] 
        cdef interval_item p
        #while i != self._intervals.end():
        for i in range(self._intervals.size()):
            p = self._intervals.at(i)
            result.append((<int>p.start,<int>p.end,<float>p.value,<float>p.length)) 
        return result
    
    @classmethod
    def from_list(cls, list l):
        p = cPairIBD()
        for interval in l:
            p.add_interval(interval[0],interval[1])
        return p
    
    cpdef merge_intervals(cPairIBD self):
        tree = self.__get_interval_tree()
        new_tree = IntervalTree()
        intervals = cPairIBD._get_tree_list(tree)
        for interval in intervals:
            intersections = [(x.start,x.end) for x in new_tree.find(interval[0]-1,interval[1]+1)]
            if len(intersections) > 0:
                #non_intersecting = []; 
                #new_tree.traverse(lambda x: non_intersecting.append((x.start,x.end)))
                non_intersecting = cPairIBD._get_tree_list(new_tree)
                for intersection in intersections:
                    non_intersecting.remove(intersection)
                    
                new_tree = IntervalTree()
                for non in non_intersecting:
                    new_tree.add_interval(Interval(non[0],non[1]))
                
                intersections.append(interval)
                new_tree.add_interval(Interval(min([x[0] for x in intersections]),max([x[1] for x in intersections])))
            else:
                new_tree.add_interval(Interval(interval[0],interval[1]))
        self._intervals.erase(self._intervals.begin(),self._intervals.end())
        l = cPairIBD._get_tree_list(new_tree)
        cdef interval_item p
        for interval in l:
            p.start = <int>interval[0]
            p.end = <int>interval[1]
            self._intervals.push_back(p)
    
    @classmethod
    def _get_tree_list(cls,tree):
        l = []
        tree.traverse(lambda x: l.append((x.start,x.end)))
        return l
    
    cpdef get_total_segments_length(self):
        l = self.to_list()
        total_length = 0
        for segment in l:
            total_length += (segment[1]-segment[0])
        return total_length
    
    cpdef get_num_windows(self, int window_size):
        l = self.to_list()
        total_length = 0
        for segment in l:
            total_length += round((segment[1]-segment[0])/float(window_size))
        return total_length
    
    cpdef get_IBD_percent(cPairIBD self, int snp_num):
        return  100 * float(sum([i[1]-i[0] for i in self.to_list()])) / snp_num
    
#     cpdef calc_accuracy(cPairIBD self, cPairIBD true_ibd, int snp_num):
#         if true_ibd == None:
#             return 0
#         intersections = iu.IntersectIntervalTrees(self._tree, true_ibd.to_list())
#         TP = sum([i[1]-i[0] for i in intersections])
#         FP = sum([i[1]-i[0] for i in self.to_list()]) - TP
#         FN = sum([i[1]-i[0] for i in true_ibd.to_list()]) - TP
#         TN = snp_num - (TP + FP + FN)
#         accuracy = 100 * float(TP + TN) / snp_num
#         return accuracy
   
    cpdef stats(cPairIBD self, cPairIBD true_ibd, int snp_num):
        if true_ibd == None:
            return 0
        self.merge_intervals()
        true_ibd.merge_intervals()
        true_segments = true_ibd.to_list()
        tree = self._tree
        TP = 0
        FP = 0
        TN = 0
        FN = 0
        intersections = iu.IntersectIntervalTrees(self._tree, true_ibd.to_list()) 
        if len(intersections) > 0:
            for inter in intersections:
                TP += inter[1] - inter[0]
        FP = self.get_total_segments_length() - TP
        FN = true_ibd.get_total_segments_length() - TP
        TN = snp_num - (TP + FP + FN)
        power = TP/(TP+FN)
        specificity = TN/(TN+FP)
        FPR = 1 - specificity
        FDR = FP/(TP+FP)
        return (TP,FP,TN,FN,power,FDR,FPR)
    
    cpdef stats_win(cPairIBD self, cPairIBD true_ibd, int snp_num, int window_size):
        if true_ibd == None:
            return 0
        self.merge_intervals()
        true_ibd.merge_intervals()
        true_segments = true_ibd.to_list()
        tree = self._tree
        TP = 0
        FP = 0
        TN = 0
        FN = 0
        intersections = iu.IntersectIntervalTrees(self._tree, true_ibd.to_list()) 
        if len(intersections) > 0:
            for inter in intersections:
                TP += round((inter[1] - inter[0])/float(window_size))
        FP = self.get_num_windows(window_size) - TP
        FN = true_ibd.get_num_windows(window_size) - TP
        TN = snp_num/window_size - (TP + FP + FN)
        power = TP/(TP+FN) if TP+FN>0 else 0
        specificity = TN/(TN+FP) if TN+FP>0 else 1
        FPR = 1 - specificity
        FDR = FP/(TP+FP)
        return (TP,FP,TN,FN,power,FDR,FPR)
   
#     cpdef get_detected_segments(cPairIBD self, cPairIBD true_ibd):
#         true_segments = true_ibd.to_list()
#         cdef list detected = []
#         tree = self._tree
#         for segment in true_segments:
#             intersections = tree.find(segment[0],segment[1])
#             if len(intersections) > 0:
#                 detected.append(segment)
#         return detected
#     
#     cpdef get_detected_segments_inv(cPairIBD self, cPairIBD true_ibd):
#         self_segments = self.to_list()
#         cdef list detected = []
#         for segment in self_segments:
#             intersections = true_ibd.find(segment[0],segment[1])
#             if len(intersections) > 0:
#                 detected.append(segment)
#         return detected
    
    cpdef add_ibd_from_containers(cPairIBD self, f1, f2, min_ibd_length = 0): 
        for founder in f1._founder_to_tree.keys():
            if f2._founder_to_tree.has_key(founder):
                intersection =\
                iu.IntersectIntervalTrees(f1._founder_to_tree[founder], iu.get_interval_list(f2._founder_to_tree[founder]))
                for inter in intersection:
                    if inter[1]-inter[0] > min_ibd_length:
                        self.add_interval(inter[0], inter[1])
        self.merge_intervals()
        
    cpdef add_ibd_from_containers_only_admixed(cPairIBD self, f1, f2, min_ibd_length, f1_other, f2_other, ceu_inds, yri_inds): 
        for founder in f1._founder_to_tree.keys():
            if f2._founder_to_tree.has_key(founder):
                intersection =\
                iu.IntersectIntervalTrees(f1._founder_to_tree[founder], iu.get_interval_list(f2._founder_to_tree[founder]))
                for inter in intersection:
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
        self.merge_intervals()
        
    cpdef filter_by_length(self, min_length, max_length):
    
        cdef cpp_vector[interval_item] *new_intervals = new cpp_vector[interval_item]()
        cdef interval_item p
        #while i != self._intervals.end():
        for i in range(self._intervals.size()):
            p = self._intervals.at(i)
            if p.length >= min_length and p.length <= max_length:
                new_intervals.push_back(p) 
        self._intervals = new_intervals
        
    cpdef filter_by_score(self, min_score, max_score):
    
        cdef cpp_vector[interval_item] *new_intervals = new cpp_vector[interval_item]()
        cdef interval_item p
        #while i != self._intervals.end():
        for i in range(self._intervals.size()):
            p = self._intervals.at(i)
            if p.value >= min_score and p.value <= max_score:
                new_intervals.push_back(p) 
        self._intervals = new_intervals
    
    cpdef calc_dists(self, dists):
        cdef interval_item p
        #while i != self._intervals.end():
        for i in range(self._intervals.size()):
            self._intervals.at(i).length = dists[self._intervals.at(i).end] - dists[self._intervals.at(i).start] 

cdef class cPopulationIBD:
    
    cdef dict _ibd_dic
    
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
        assert isinstance(other_ibd, cPopulationIBD)
        for pair in other_ibd.keys():
            if self.has_key(pair):
                self._ibd_dic[pair].update(other_ibd.get_value(pair))
            else:
                self.add_human_pair(pair, other_ibd.get_value(pair))
    
    cpdef keys(self):
        return self._ibd_dic.keys()
        
    cpdef merge_all(self):
        for pair in self.keys():
            self.get_value(pair).merge_intervals()
    
    cpdef stats(self, true_ibd, snp_num):
        stats = []
        for pair in true_ibd.keys():
            if self.has_key(pair):
                stats.append(self.get_value(pair).stats(true_ibd.get_value(pair),snp_num))
            else:
                TN = snp_num - true_ibd.get_value(pair).get_total_segments_length()
                FN = true_ibd.get_value(pair).get_total_segments_length()
                stats.append((0,0,TN,FN,0,0,0))
        #print stats 
        power = numpy.mean([x[4] for x in stats])
        FDR = numpy.mean([x[5] for x in stats])
        FPR = numpy.mean([x[6] for x in stats])
        f = open("stats.txt","w")
        f.writelines(string.join([str(x) for x in stats],"\n"))
        f.close()
        return (power,FDR,FPR)
    
    cpdef stats_win(self, true_ibd, snp_num, window_size):
        stats = []
        for pair in true_ibd.keys():
            if self.has_key(pair):
                stats.append(self.get_value(pair).stats_win(true_ibd.get_value(pair),snp_num,window_size))
            else:
                TN = snp_num/window_size - true_ibd.get_value(pair).get_num_windows(window_size)
                FN = true_ibd.get_value(pair).get_num_windows(window_size)
                stats.append((0,0,TN,FN,0,0,0))
        #print stats 
        power = numpy.mean([x[4] for x in stats])
        FDR = numpy.mean([x[5] for x in stats])
        FPR = numpy.mean([x[6] for x in stats])
        f = open("stats.txt","w")
        f.writelines(string.join([str(x) for x in stats],"\n"))
        f.close()
        return (power,FDR,FPR)
    
    cpdef get_all_segments(self):
        segments = []
        for pair in self._ibd_dic.keys():
            segments.append(self.get_value(pair).to_list())
        return segments 
    
#     cpdef calc_power_intervals(self, true_ibd, bins=[0,500,1000,1500,2000,2500,10000]):
#         '''
#         bins defines the bin edges, inluding the leftmost and rightmost edges.
#         Returns an array of size len(bins)-1 with the power values:
#         power[i] is the power to detect a segment of length >= bins[i-1] and <bins[i]  
#         '''
#         assert isinstance(true_ibd, cPopulationIBD)
#         detected_all = []
#         detected_segments_lengths = []
#         all_segments_lengths = []
#         for pair in true_ibd.keys():
#             if self.has_key(pair):
#                 detected = self._ibd_dic[pair].get_detected_segments(true_ibd.get_value(pair))
#                 detected_all += detected
#                 detected_segments_lengths += [x[1]-x[0] for x in detected]
#                 all_segments_lengths += [x[1]-x[0] for x in true_ibd.get_value(pair).to_list()]
#          
#         (detected_hist, dummy) = numpy.histogram(detected_segments_lengths, bins)
#         (all_hist, dummy) = numpy.histogram(all_segments_lengths, bins)
#         return (numpy.nan_to_num(detected_hist/all_hist)).tolist()
#     
#     cpdef calc_power(self, true_ibd):
#         assert isinstance(true_ibd, cPopulationIBD)
#         detected_all = []
#         for pair in true_ibd.keys():
#             if self.has_key(pair):
#                 detected = self._ibd_dic[pair].get_detected_segments(true_ibd.get_value(pair))
#                 detected_all += detected
#          
#         if len(detected_all) == 0:
#             return 1
#         return float(len(detected_all)) / len(true_ibd.get_all_segments())
#     
#     cpdef calc_false_positives(self, true_ibd):
#         total_segments_num = 0
#         false_positives = 0
#         for pair in self.keys():
#             segments_num = len(self._ibd_dic[pair].to_list())
#             total_segments_num += segments_num
#             if true_ibd.has_key(pair):
#                 detected = self._ibd_dic[pair].get_detected_segments(true_ibd.get_value(pair))
#                 false_positives += segments_num - len(detected)
#             else:
#                 false_positives += segments_num
#          
#         if total_segments_num == 0:
#             return 1
#         return float(false_positives) / total_segments_num
#     
#     cpdef calc_false_positives_intervals(self, true_ibd, bins=[0,500,1000,1500,2000,2500,10000]):
#         all_segments = []
#         false_positives = []
#         for pair in self.keys():
#             segments = self._ibd_dic[pair].to_list()
#             all_segments += segments
#             if true_ibd.has_key(pair):
#                 detected = self._ibd_dic[pair].get_detected_segments_inv(true_ibd.get_value(pair))
#                 for seg in detected:
#                     if seg in segments:
#                         segments.remove(seg)
#             false_positives += segments
#          
#         all_segments_lengths = [x[1] - x[0] for x in all_segments]
#         false_positives_lengths = [x[1] - x[0] for x in false_positives]
#         (false_positives_hist, dummy) = numpy.histogram(false_positives_lengths, bins)
#         (all_hist, dummy) = numpy.histogram(all_segments_lengths, bins)
#         return (numpy.nan_to_num(false_positives_hist / all_hist)).tolist()
         
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
                
    cpdef filter_by_length(self, min_length, max_length):
        for pair in self.keys():
            self._ibd_dic[pair].filter_by_length(min_length, max_length)
            if len(self._ibd_dic[pair].to_list()) == 0:
                self._ibd_dic.pop(pair)
                
    cpdef filter_by_score(self, min_score, max_score):
        for pair in self.keys():
            self._ibd_dic[pair].filter_by_score(min_score,max_score)
            if len(self._ibd_dic[pair].to_list()) == 0:
                self._ibd_dic.pop(pair)
            
    cpdef calc_dists(self, dists):
        for pair in self.keys():
            self._ibd_dic[pair].calc_dists(dists)
    
    cpdef to_dict(self):
        d = {}
        for pair in self._ibd_dic.keys():
            d[pair] = self._ibd_dic[pair].to_list()
        return d
    
    @classmethod
    def from_dict(cls, d):
        p = cPopulationIBD()
        for pair in d.keys():
            p.add_human_pair(pair,cPairIBD.from_list(d[pair]))
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
                s += [str(inter[0]),",",str(inter[1]),";"]
            if len(intervals)>0:
                s.pop()
            s.append("\n")
        return "".join(s)
    
    @classmethod
    def from_string(cls, s, dists=[]):
        assert isinstance(s, (str,list))
        p = cPopulationIBD()
        if type(s) is str:
            s = string.split(s, '\n')
            s = [x for x in s if x != '']
        s = [string.strip(x, '\n') for x in s]
        for pair_string in s:
            temp = pair_string.split(":")
            pair = temp[0].split(",")
            pair = (int(pair[0]),int(pair[1]))
            pairIBD = cPairIBD()
            for inter in temp[1].split(";"):
                points = inter.split(",")
                score = 0
                if len(points) > 2:
                    score = float(points[2])
                dist = 0
                if len(dists) > 0:
                    if int(points[1]) >= len(dists) or int(points[0]) >= len(dists):
                        print pair, points, len(dists)
                    dist = dists[int(points[1])] - dists[int(points[0])]
                pairIBD.add_interval(int(points[0]),int(points[1]),score,dist)
            p.add_human_pair(pair,pairIBD)
            p.get_value(pair).calc_dists(dists)
        return p
    
    @classmethod
    def fast_deserialize(cls, file_name, gm_file_name=None):
        dists = []
        if gm_file_name != None:
            gm_f = open(gm_file_name)
            data = gm_f.readlines()
            dists = [float(x.split(" ")[2]) for x in data]
        f = open(file_name)
        s = f.readlines()
        ibd = cPopulationIBD.from_string(s,dists)
        #ibd.merge_all()
        f.close()
        return ibd
    
    