'''
Created on Mar 5, 2012

@author: Eskin3
'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext 
import numpy as np
import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("utils.StringUtils", ["utils/StringUtils.pyx"], language="c++", ), 
                   Extension("utils.FoundersContainer", ["utils/FoundersContainer.pyx"], language="c++"), 
                   Extension("IBD.intersection", ["IBD/intersection.pyx"], language="c++", include_dirs=['.'], extra_compile_args=["-g"], extra_link_args=["-g"]),
                   Extension("IBD.cIBD", ["IBD/cIBD.pyx", "IBD/intervalItem.cpp"], language="c++", include_dirs=['.', 'IBD'], extra_compile_args=["-g"], extra_link_args=["-g"]), 
                   Extension("IBD.GeneticMap", ['IBD/GeneticMap.pyx', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-g"], extra_link_args=["-g"]),
                   Extension("IBD.InnerModel", ['IBD/InnerModel.pyx', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-g"], extra_link_args=["-g"]),
                   Extension("IBD.GenotypePairModel", ['IBD/GenotypePairModel.pyx', 'IBD/InnerModel.pxd', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-g"], extra_link_args=["-g"]),
                   Extension("IBD.LDModel", ['IBD/LDModel.pyx', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-g"], extra_link_args=["-g"]),
                   Extension("IBD.TestSet", ['IBD/TestSet.pyx', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-g"], extra_link_args=["-g"])
                   ]
)

#Extension("IBDGenoHMM", ['IBDGenoHMM.pyx', 'LDModel.pyx', 'structs.cpp'], language="c++", include_dirs=['.', np.get_include()], extra_compile_args=["-g"], extra_link_args=["-g"])