# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 10:51:05 2015

@author: John
"""


from distutils.core import setup
from Cython.Build import cythonize



setup(
    ext_modules=cythonize('def_distance.pyx'),
)