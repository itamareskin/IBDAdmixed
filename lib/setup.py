'''
Created on Mar 5, 2012

@author: Eskin3
'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext 
import numpy as np
import Cython.Compiler.Options
import os
Cython.Compiler.Options.annotate = True

os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("utils.StringUtils", ["utils/StringUtils.pyx"], language="c++", extra_compile_args=["-O3"]), 
                   Extension("utils.FoundersContainer", ["utils/FoundersContainer.pyx"], language="c++", extra_compile_args=["-O3"]),
                   Extension("utils.FormatConversions", ['utils/FormatConversions.pyx', 'utils/ped_to_bgl.cpp'], language="c++", include_dirs=['.'], extra_compile_args=["-O3"]),
                   #Extension("utils.GermlineWrapper", ['utils/GermlineWrapper.pyx'], language="c++", include_dirs=['.'], extra_compile_args=["-O3"]),
                   Extension("utils.GermlineWrapper", ['utils/GermlineWrapper.pyx', '../external/germline-1-5-1/GERMLINE_0001.cpp', '../external/germline-1-5-1/GERMLINE.cpp', '../external/germline-1-5-1/Share.cpp', '../external/germline-1-5-1/Chromosome.cpp', '../external/germline-1-5-1/ChromosomePair.cpp', '../external/germline-1-5-1/HMIndividualsExtractor.cpp', '../external/germline-1-5-1/MarkerSet.cpp', '../external/germline-1-5-1/Individual.cpp', '../external/germline-1-5-1/Individuals.cpp', '../external/germline-1-5-1/InputManager.cpp', '../external/germline-1-5-1/MatchFactory.cpp', '../external/germline-1-5-1/MatchesBuilder.cpp', '../external/germline-1-5-1/NucleotideMap.cpp', '../external/germline-1-5-1/PEDIndividualsExtractor.cpp', '../external/germline-1-5-1/Match.cpp', '../external/germline-1-5-1/PolymorphicIndividualsExtractor.cpp', '../external/germline-1-5-1/SNP.cpp', '../external/germline-1-5-1/SNPPositionMap.cpp', '../external/germline-1-5-1/SNPs.cpp'], language="c++", include_dirs=['.','../external/germline-1-5-1/include'], extra_compile_args=["-g", "-O3"]),
                   Extension("IBD.IntervalTree", ["IBD/IntervalTree.pyx"], language="c++", include_dirs=['.'], extra_compile_args=["-O3"]),
                   Extension("IBD.GeneticMap", ['IBD/GeneticMap.pyx', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-O3"]),
                   Extension("IBD.IBDSegments", ["IBD/IBDSegments.pyx", 'IBD/IntervalTree.pxd', 'IBD/GeneticMap.pxd'], language="c++", include_dirs=['.', 'IBD'], extra_compile_args=["-O3"]),
                   Extension("IBD.LDModel", ['IBD/LDModel.pyx', 'IBD/GeneticMap.pxd', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-O3"]),
                   Extension("IBD.TestSet", ['IBD/TestSet.pyx', 'IBD/GeneticMap.pxd', 'IBD/LDModel.pxd', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-O3"]),
                   Extension("IBD.InnerModel", ['IBD/InnerModel.pyx', 'IBD/GeneticMap.pxd', 'IBD/TestSet.pxd', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-O3"]),
                   Extension("IBD.GenotypePairModel", ['IBD/GenotypePairModel.pyx', 'IBD/GeneticMap.pxd',  'IBD/TestSet.pxd', 'IBD/LDModel.pxd', 'IBD/InnerModel.pxd', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-O3"]),
                   Extension("IBD.GenotypeModel", ['IBD/GenotypeModel.pyx', 'IBD/GeneticMap.pxd',  'IBD/TestSet.pxd', 'IBD/LDModel.pxd', 'IBD/InnerModel.pxd', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-O3"]),
                   Extension("IBD.WindowedModel", ['IBD/WindowedModel.pyx', 'IBD/GeneticMap.pxd', 'IBD/TestSet.pxd', 'IBD/InnerModel.pxd', 'IBD/structs.cpp'], language="c++", include_dirs=['.', 'IBD', np.get_include()], extra_compile_args=["-O3"])
                   ]
)

#Extension("IBDGenoHMM", ['IBDGenoHMM.pyx', 'LDModel.pyx', 'structs.cpp'], language="c++", include_dirs=['.', np.get_include()], extra_compile_args=["-g"], extra_link_args=["-g"])