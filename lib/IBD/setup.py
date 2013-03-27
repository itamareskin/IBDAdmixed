'''
Created on Mar 5, 2012

@author: Eskin3
'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext 
import numpy as np

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("StringUtils", ["StringUtils.pyx"], language="c++", ), 
                   Extension("FoundersContainer", ["FoundersContainer.pyx"], language="c++"), 
                   Extension("cIBD", ["cIBD.pyx"], language="c++"), 
                   Extension("LDModel", ['LDModel.pyx', 'structs.cpp'], language="c++", include_dirs=['.', np.get_include()], extra_compile_args=["-g"], extra_link_args=["-g"])
                   ]
)

#Extension("IBDGenoHMM", ['IBDGenoHMM.pyx', 'LDModel.pyx', 'structs.cpp'], language="c++", include_dirs=['.', np.get_include()], extra_compile_args=["-g"], extra_link_args=["-g"])