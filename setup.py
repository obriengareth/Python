# -*- coding: utf-8 -*-
"""
LSM 1 - gareth 
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy as np

setup(ext_modules = cythonize('kirchhoff_migration.pyx'))

#ext = Extension("mod_cyMigrate", ["mod_cyMigrate.pyx"], include_dirs = [np.get_include()])
#setup(ext_modules=[ext], cmdclass = {'build_ext': build_ext})

