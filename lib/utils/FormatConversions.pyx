from libcpp cimport bool
import os
from itertools import islice

cdef extern from "../utils/ped_to_bgl.h":
    int ped2bgl(char *ped_file, char *map_file, char *bgl_file, char *markers_file)

cdef extern from "../utils/ped_to_lamp.h":
    int ped2lamp(char *ped_file, char *map_file, char *lamp_file, char *markers_file, bool phased)

def convert_ped_to_bgl(ped_file, map_file, bgl_file, markers_file):
    return ped2bgl(ped_file, map_file, bgl_file, markers_file)

def convert_ped_to_lamp(ped_file, map_file, lamp_file, markers_file, phased):
    return ped2lamp(ped_file, map_file, lamp_file, markers_file, phased)

def convert_map_to_markers(map_file, markers_file):
    with open(map_file) as file_map, open(markers_file, "w") as file_markers:
        done = False
        buffer_size = 100
        snp_idx = 0
        while True:
            if done:
                break
            # read next buffer_size lines from the file
            lines = list(islice(file_map, buffer_size))
            if len(lines) == 0:
                done = True
            for line in lines:
                line = line.split()
                file_markers.write(line[1] + " " + line[3] + " 1 2")
                snp_idx += 1