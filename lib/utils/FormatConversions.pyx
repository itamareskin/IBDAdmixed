from libcpp cimport bool

cdef extern from "../utils/ped_to_bgl.h":
    int ped2bgl(char *ped_file, char *map_file, char *bgl_file, char *markers_file)

cdef extern from "../utils/ped_to_lamp.h":
    int ped2lamp(char *ped_file, char *map_file, char *lamp_file, char *markers_file, bool phased)

def convert_ped_to_bgl(ped_file, map_file, bgl_file, markers_file):
    return ped2bgl(ped_file, map_file, bgl_file, markers_file)

def convert_ped_to_lamp(ped_file, map_file, lamp_file, markers_file, phased):
    return ped2lamp(ped_file, map_file, lamp_file, markers_file, phased)