from libc.stdlib cimport malloc, free
from libc.string cimport strcmp
from cpython.string cimport PyString_AsString
import subprocess
import os

# from libcpp cimport bool
# cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.cpp":
#     double MIN_MATCH_LEN
# cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.cpp":
#     int MARKER_SET_SIZE
# cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.cpp":
#     bool PRINT_MATCH_HAPS
# cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.cpp":
#     bool ROI
# cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.cpp":
#     bool HAP_EXT
# cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.cpp":
#     bool WIN_EXT
# cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.cpp":
#     bool ALLOW_HOM
# cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.cpp":
#     bool HOM_ONLY
# cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.cpp":
#     bool HAPLOID
# cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.cpp":
#     bool SILENT
# cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.cpp":
#     bool DEBUG
# cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.cpp":
#     bool BINARY_OUT
# cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.cpp":
#     int MAX_ERR_HOM
# cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.cpp":
#     int MAX_ERR_HET

cdef extern from "../../external/germline-1-5-1/GERMLINE_0001.h":
    int germline_main(int argc, char* argv[])

cdef char ** to_cstring_array(list_str):
    cdef char **ret = <char **>malloc(len(list_str) * sizeof(char *))
    for i in xrange(len(list_str)):
        ret[i] = PyString_AsString(list_str[i])
    return ret

def germline(args):
    #germline_path = os.path.join(os.path.split(os.path.split(os.path.split(__file__)[0])[0])[0], "external", "germline-1-5-1","germline.exe")
    #print germline_path
    #command = [germline_path]+args
    #print command
    #retcode = subprocess.call(args, executable=germline_path)
    return germline_main(len(args), to_cstring_array(args))