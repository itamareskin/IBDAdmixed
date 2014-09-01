from __future__ import division
import os, sys
from IBD.GeneticMap import GeneticMap

gm = GeneticMap(sys.argv[2])
with open(sys.argv[1]) as input_file:
    print "anc0 anc1 anc2 tot_length avg_lengths\n"
    while True:
        prev_anc1 = 0
        prev_anc2 = 0
        prev_breakpoint1=0
        prev_breakpoint2=0
        anc_lens = [0]*3
        segment_lengths = []
        prev_loc = 0
        line = input_file.readline()
        if len(line) == 0 or line is None:
            break
        line = line.split()
        for segment in line:
            segment = segment.split(":")
            anc1 = int(segment[0][0])
            anc2 = int(segment[0][1])
            anc_lens[anc1] += gm.get_length_bp(prev_loc, int(segment[1]))
            anc_lens[anc2] += gm.get_length_bp(prev_loc, int(segment[1]))
            if anc1 != prev_anc1:
                prev_anc1 = anc1
                segment_lengths.append(gm.get_length_bp(prev_breakpoint1, int(segment[1])))
                prev_breakpoint1 = int(segment[1])
            if anc2 != prev_anc2:
                prev_anc2 = anc2
                segment_lengths.append(gm.get_length_bp(prev_breakpoint2, int(segment[1])))
                prev_breakpoint2 = int(segment[1])
            prev_loc = int(segment[1])
        tot = sum(anc_lens)
        print str(anc_lens[0]/tot) + " " + str(anc_lens[1]/tot) + " " + str(anc_lens[2]/tot) + " " + str(sum(segment_lengths)) + " " + str(sum(segment_lengths)/len(segment_lengths))




