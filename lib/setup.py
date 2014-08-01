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
    ext_modules = [Extension("utils.StringUtils", ["utils/StringUtils.pyx"], language="c++", extra_compile_args=["-O3"]), 
                   Extension("utils.FoundersContainer", ["utils/FoundersContainer.pyx"], language="c++", extra_compile_args=["-O3"]), 
                   Extension("IBD.IntervalTree", ["IBD/IntervalTree.pyx"], language="c++", include_dirs=['.'], extra_compile_args=["-O3"]),
                   Extension("IBD.GeneticMap", ['IBD/GeneticMap.pyx', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-O3"]),
                   Extension("IBD.cIBD", ["IBD/cIBD.pyx", 'IBD/IntervalTree.pxd', 'IBD/GeneticMap.pxd'], language="c++", include_dirs=['.', 'IBD'], extra_compile_args=["-O3"]),
                   Extension("IBD.LDModel", ['IBD/LDModel.pyx', 'IBD/GeneticMap.pxd', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-O3"]),
                   Extension("IBD.TestSet", ['IBD/TestSet.pyx', 'IBD/GeneticMap.pxd', 'IBD/LDModel.pxd', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-O3"]),
                   Extension("IBD.InnerModel", ['IBD/InnerModel.pyx', 'IBD/GeneticMap.pxd', 'IBD/TestSet.pxd', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-O3"]),
                   Extension("IBD.GenotypePairModel", ['IBD/GenotypePairModel.pyx', 'IBD/GeneticMap.pxd',  'IBD/TestSet.pxd', 'IBD/LDModel.pxd', 'IBD/InnerModel.pxd', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-O3"]),
                   Extension("IBD.GenotypeModel", ['IBD/GenotypeModel.pyx', 'IBD/GeneticMap.pxd',  'IBD/TestSet.pxd', 'IBD/LDModel.pxd', 'IBD/InnerModel.pxd', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-O3"]),
                   Extension("IBD.WindowedModel", ['IBD/WindowedModel.pyx', 'IBD/GeneticMap.pxd', 'IBD/TestSet.pxd', 'IBD/InnerModel.pxd', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-O3"])
                   ]
)

#Extension("IBDGenoHMM", ['IBDGenoHMM.pyx', 'LDModel.pyx', 'structs.cpp'], language="c++", include_dirs=['.', np.get_include()], extra_compile_args=["-g"], extra_link_args=["-g"])