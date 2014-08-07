
cdef extern from "../utils/ped_to_bgl.h":
    int ped2bgl(char *ped_file, char *map_file, char *bgl_file, char *markers_file)

def convert_ped_to_bgl(ped_file, map_file, bgl_file, markers_file):
    return ped2bgl(ped_file, map_file, bgl_file, markers_file)